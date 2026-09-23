"""Production-schema storage contract and optional C++ cross-language checks."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
import state_archive as archive
import h5py
import numpy as np


class StateArchiveTest(unittest.TestCase):
    def test_native_bjt_seed_converts_density_units_and_binds_mesh(self):
        import run_genius_bjt_mesh_sensitivity as bjt
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);(root/'fields').mkdir()
            mesh=dict(nodes=[dict(id=0,x=0.,y=0.),dict(id=1,x=1.,y=0.)],
                      triangles=[],regions=[],contacts=[])
            meshpath=root/'mesh.json';meshpath.write_text(json.dumps(mesh))
            for name,value in [('ElectrostaticPotential',.1),('eQuasiFermiPotential',.2),
                               ('hQuasiFermiPotential',.3),('eDensity',1e18),('hDensity',1e8)]:
                (root/'fields'/f'{name}_region0.csv').write_text(
                    f'node_id,component0\n0,{value}\n1,{value}\n')
            output=root/'seed.h5'
            self.assertEqual(bjt.make_seed_state(root,output,meshpath),2)
            fields,meta=archive.read(output,2,archive.mesh_identity(mesh,1e-6))
            self.assertEqual(float(fields['electrons_m3'][0]),1e24)
            self.assertEqual(float(fields['holes_m3'][0]),1e14)
            self.assertEqual(len(meta['source_fields_sha256']),5)
            with self.assertRaises(ValueError):archive.read(output,2,archive.mesh_identity(mesh,1.))

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)/'state.h5'
        self.fields = dict(psi=[0., .125, -.25], phin=[28.,28.,28.], phip=[0.,-0.,1e-200],
            electrons_m3=[1e26,1e12,0.], holes_m3=[1e-120,1e23,1e-30],
            electron_qf_increment_V=[1e-18,-1e-18,0.], electron_qf_reference_V=[28.,28.,28.],
            hole_qf_increment_V=[0.,-0.,1e-200], hole_qf_reference_V=[0.,0.,0.])
        self.meta = dict(mode='dd', mesh_sha256='a'*64, potential_origin_V=28., bias_V=40.)

    def assertFields(self, got):
        self.assertEqual(set(got), set(self.fields))
        for name in self.fields:
            self.assertEqual(np.asarray(got[name],dtype='<f8').tobytes(),
                             np.asarray(self.fields[name],dtype='<f8').tobytes(), name)

    def test_roundtrip_thermal_and_split(self):
        self.meta['mode'] = 'electrothermal'
        self.fields['temperature_K'] = [300.,327.25,401.5]
        archive.write(self.path, self.fields, self.meta)
        got, metadata = archive.read(self.path, 3, 'a'*64)
        self.assertFields(got)
        self.assertEqual(metadata, self.meta)

    def test_row_bridge_preserves_signed_zero_and_split_low_bits(self):
        rows=archive.fields_to_rows(self.fields)
        self.assertIsInstance(rows[0]['psi'],float)
        self.assertFields(archive.rows_to_fields(rows))

    def test_reject_wrong_mesh_even_same_node_count(self):
        archive.write(self.path, self.fields, self.meta)
        with self.assertRaises(ValueError): archive.read(self.path, 3, 'b'*64)
        with self.assertRaises(ValueError): archive.read(self.path, 4, 'a'*64)

    def test_failed_commit_preserves_old_checkpoint(self):
        archive.write(self.path, self.fields, self.meta)
        before = self.path.read_bytes()
        with patch.object(archive.os, 'replace', side_effect=OSError('injected commit failure')):
            with self.assertRaises(OSError): archive.write(self.path, self.fields, self.meta)
        self.assertEqual(before, self.path.read_bytes())
        self.assertEqual(list(self.path.parent.iterdir()), [self.path])

    def test_reject_invalid_fields(self):
        self.fields.pop('electron_qf_increment_V')
        with self.assertRaises(ValueError): archive.write(self.path, self.fields, self.meta)
        self.assertFalse(self.path.exists())

    def test_translation_preserves_increment_rounding_remainder(self):
        archive.write(self.path,self.fields,self.meta)
        target=self.path.parent/'translated.h5'
        archive.translate(self.path,target,1e-18,{0,1},3,'a'*64)
        fields,meta=archive.read(target,3,'a'*64)
        self.assertEqual(meta['potential_origin_V'],28.)
        self.assertEqual(meta['bias_V']+meta['potential_origin_V'],self.meta['bias_V']+self.meta['potential_origin_V'])
        self.assertEqual(fields['electron_qf_increment_V'][0],0.)
        self.assertEqual(fields['electron_qf_increment_V'][1],-2e-18)
        self.assertEqual(fields['phin'][2],0.)
        np.testing.assert_array_equal(fields['electrons_m3'],self.fields['electrons_m3'])

    def test_reject_corrupt_units_dtype_and_temperature(self):
        for mutation in ('unit', 'type', 'temperature'):
            archive.write(self.path, self.fields, self.meta)
            with h5py.File(self.path, 'r+') as f:
                if mutation == 'unit': f['fields/psi'].attrs['unit'] = 'mV'
                if mutation == 'type':
                    del f['fields/psi']
                    f['fields'].create_dataset('psi', data=np.zeros(3,dtype='<f4')).attrs['unit'] = 'V'
                if mutation == 'temperature':
                    m = dict(self.meta, mode='electrothermal')
                    f.attrs['metadata_json'] = json.dumps(m)
            with self.assertRaises(ValueError): archive.read(self.path, 3, 'a'*64)

    @unittest.skipUnless(os.environ.get('VELA_STATE_ARCHIVE_TOOL'), 'C++ tool supplied by CTest')
    def test_cpp_python_interop(self):
        tool = os.environ['VELA_STATE_ARCHIVE_TOOL']
        self.meta['mode'] = 'electrothermal'
        self.fields['temperature_K'] = [300.,327.25,401.5]
        archive.write(self.path, self.fields, self.meta)
        output = subprocess.check_output([tool,'decode',str(self.path),'3','a'*64], text=True)
        got = json.loads(output)
        self.assertFields(got['fields'])
        self.assertEqual(got['metadata'], self.meta)
        input_path = self.path.with_suffix('.json')
        input_path.write_text(json.dumps(dict(fields=self.fields, metadata=self.meta)))
        subprocess.run([tool,'encode',str(input_path),str(self.path)], check=True)
        self.assertFields(archive.read(self.path,3,'a'*64)[0])
        with h5py.File(self.path,'r+') as f: f['fields/psi'].attrs['unit'] = 'mV'
        result = subprocess.run([tool,'decode',str(self.path),'3','a'*64], capture_output=True)
        self.assertNotEqual(result.returncode, 0)


if __name__ == '__main__': unittest.main()
