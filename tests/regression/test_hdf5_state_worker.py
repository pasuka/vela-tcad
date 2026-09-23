"""Check production-schema mesh binding and failure isolation in the DC worker."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import state_archive as archive
import numpy as np
RUNNER = None


class Hdf5WorkerTest(unittest.TestCase):
    def test_mesh_units_frame_and_failed_request_isolation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mesh = dict(nodes=[dict(id=i,x=x*1e-6,y=y*1e-6) for i,(x,y) in enumerate(((0,0),(1,0),(1,1),(0,1)))],
                triangles=[dict(id=0,region_id=0,node_ids=[0,1,2]),dict(id=1,region_id=1,node_ids=[0,2,3])],
                regions=[dict(id=0,name='n',material='Si',cell_ids=[0]),dict(id=1,name='p',material='Si',cell_ids=[1])],
                contacts=[dict(id=0,name='anode',region_id=1,node_ids=[0,3]),dict(id=1,name='cathode',region_id=0,node_ids=[1,2])])
            (root/'mesh.json').write_text(json.dumps(mesh))
            fields = dict(psi=[.11*(i+1) for i in range(4)], phin=[1e-17]*4,phip=[-1e-17]*4,
                electrons_m3=[1e10]*4,holes_m3=[1e10]*4,electron_quantum_potential_V=[0.]*4,
                electron_qf_increment_V=[1e-17]*4,hole_qf_increment_V=[-1e-17]*4,
                electron_qf_reference_V=[0.]*4,hole_qf_reference_V=[0.]*4)
            identity=archive.mesh_identity(mesh,1.)
            metadata=dict(mode='dd',mesh_sha256=identity,potential_origin_V=0.)
            archive.write(root/'seed.h5',fields,metadata)
            archive.write(root/'wrong_mesh.h5',fields,dict(metadata,mesh_sha256='b'*64))
            archive.write(root/'wrong_origin.h5',fields,dict(metadata,potential_origin_V=28.))
            requests=[]
            for i,seed in enumerate(('seed.h5','wrong_mesh.h5','seed.h5','wrong_origin.h5','seed.h5')):
                cfg=dict(simulation_type='dc_sweep',state_format='hdf5',mesh_file='mesh.json',output_csv=f'curve{i}.csv',
                    doping=[dict(region='n',donors=1e23,acceptors=0),dict(region='p',donors=0,acceptors=1e23)],
                    contacts=[dict(name='anode',bias=0),dict(name='cathode',bias=0)],
                    solver=dict(method='frozen_state',impact_ionization='none'),
                    sweep=dict(contact='anode',current_contact='anode',start=0.,stop=0.,step=.1,write_vtk=False,
                        initial_state_file=seed,write_state_file=f'out{i}.h5',
                        write_state_every_point_prefix=f'point{i}'))
                path=root/f'config{i}.json';path.write_text(json.dumps(cfg))
                requests.append(json.dumps(dict(id=i,config=str(path))))
            requests.append('{"shutdown":true}')
            result=subprocess.run([str(RUNNER),'--dc-worker'],input='\n'.join(requests)+'\n',
                                  capture_output=True,text=True,timeout=60)
            self.assertEqual(result.returncode,0,result.stderr)
            responses=[json.loads(line) for line in result.stdout.splitlines()]
            for i in (0,2,4):
                self.assertEqual(responses[i]['returncode'],0,responses[i])
                got,meta=archive.read(root/f'out{i}.h5',4,identity)
                self.assertEqual(meta['potential_origin_V'],0.)
                self.assertEqual(meta['bias_V'],0.)
                self.assertEqual(meta['contact_biases_V']['anode'],0.)
                self.assertEqual(meta['state_role'],'accepted')
                self.assertIn('mesh_file',meta['input_file_sha256'])
                for name in fields:
                    self.assertEqual(np.asarray(got[name],dtype='<f8').tobytes(),np.asarray(fields[name],dtype='<f8').tobytes(),name)
                self.assertTrue(list(root.glob(f'point{i}*.h5')))
                self.assertFalse(list(root.glob(f'point{i}*.csv')))
            for i in (1,3):
                self.assertIn('error',responses[i])
                self.assertFalse((root/f'out{i}.h5').exists())


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--runner',type=Path,required=True)
    args,rest=parser.parse_known_args();RUNNER=args.runner.resolve()
    unittest.main(argv=[__file__]+rest)
