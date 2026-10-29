"""P2 real launcher checks, source identity and fail-closed interface controls."""
from __future__ import annotations
import copy
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'scripts'))
import check_sde_sdevice_cases as checks
import sentaurus_sde_cases as p1
from vela_sdevice.mesh_import import import_grid, MeshImportError
from tests.regression.gmsh_test_runtime import require_gmsh_runtime
from tests.regression.test_sdevice_mesh_import import OBTUSE_GRD, OBTUSE_DAT


class InputContractTest(unittest.TestCase):
    def test_eight_named_inputs_are_pinned(self):
        data=checks.load_cases()
        self.assertEqual(8,len(data['cases']))
        self.assertEqual({'pn2d','schottky','gate_reference'},{c['mesh_case'] for c in data['cases']})

    def test_only_grid_parameter_file_semantics_are_changed(self):
        for case in checks.load_cases()['cases']:
            original=(REPO/case['command']['path']).read_text(encoding='utf-8')
            effective,_=checks.effective_text(original)
            a=checks.parse_cmd_text(original,path_name='device.cmd')
            b=checks.parse_cmd_text(effective,path_name='device.cmd')
            for key in checks.SEMANTIC_KEYS:self.assertEqual(a[key],b[key],(case['id'],key))
            self.assertEqual('device.grd',b['files']['grid']['resolved'])
            self.assertEqual('models.par',b['files']['parameter']['resolved'])

    def test_original_missing_parameter_is_not_relabelled_as_supported(self):
        case=checks.load_cases()['cases'][0]
        text=(REPO/case['command']['path']).read_text()
        self.assertIn('E301',checks.parse_cmd_text(text,path_name='device.cmd')['verdict']['fatal_ids'])

    def test_ambiguous_multiple_grid_declarations_are_rejected(self):
        with self.assertRaisesRegex(ValueError,'one explicit Grid'):
            checks.effective_text('File { Grid="a" Grid="b" }')


class InterfaceImportTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(dir=REPO/'build-release',prefix='p2 interface ')
        self.root=Path(self.temp.name);self.addCleanup(self.temp.cleanup)

    def grid(self,record='1 2 0'):
        text=OBTUSE_GRD.replace('nb_elements = 5','nb_elements = 6').replace('nb_regions = 3','nb_regions = 4')
        text=text.replace('"R.Si" "R.Ox" "gate"','"R.Si" "R.Ox" "gate" "R.Si+R.Ox"')
        text=text.replace('Silicon Oxide Contact','Silicon Oxide Contact Interface')
        text=text.replace('Elements (5) {','Elements (6) {').replace(' 1 3 4\n  }',' 1 3 4\n '+record+'\n  }')
        i=text.rfind('}')
        text=text[:i]+'Region ("R.Si+R.Ox") {\n material = Interface\n Elements (1) {\n 5\n }\n}\n'+text[i:]
        path=self.root/'grid.grd';path.write_text(text)
        path.with_suffix('.dat').write_text(OBTUSE_DAT)
        return path

    def test_interface_is_preserved_without_becoming_a_material(self):
        report=import_grid(self.grid(),self.root/'out')
        self.assertEqual(['Si','SiO2'],report['materials'])
        self.assertEqual(2,len(report['regions']))
        mesh=json.loads((self.root/'out/mesh.json').read_text())
        self.assertEqual([{'name':'R.Si+R.Ox','node_ids':[0,2],'edge_node_ids':[[2,0]],
                          'adjacent_regions':['R.Ox','R.Si']}],mesh['interfaces'])
        self.assertEqual(4,len(mesh['triangles']))

    def test_non_line_interface_is_not_silently_discarded(self):
        with self.assertRaisesRegex(MeshImportError,'non-line'):
            import_grid(self.grid('2 0 1 2'),self.root/'bad')

    def test_unattached_interface_is_rejected(self):
        with self.assertRaisesRegex(MeshImportError,'two volume regions'):
            import_grid(self.grid('1 0 1'),self.root/'bad')


class NativeCheckTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        require_gmsh_runtime()
        cls.launcher=Path(os.environ.get('VELA_SD_SIM',REPO/'build-release'/('sd-sim.exe' if os.name=='nt' else 'sd-sim')))
        if not cls.launcher.is_file():raise unittest.SkipTest('build sd-sim or set VELA_SD_SIM')
        cls.temp=tempfile.TemporaryDirectory(dir=REPO/'build-release',prefix='p2 native paths with spaces ')
        cls.root=Path(cls.temp.name);cls.addClassCleanup(cls.temp.cleanup)
        cls.meshes=Path(os.environ.get('VELA_REFERENCE_SDE_CASES_DIR',cls.root/'generated meshes'))
        cls.runtime=cls.root/'gmsh runtime'
        result=(json.loads((cls.meshes/'case_summary.json').read_text()) if os.environ.get('VELA_REFERENCE_SDE_CASES_DIR') else
                p1.run_suite(p1.DEFAULT_MANIFEST,cls.meshes,cls.runtime,repeat=1))
        if result['status']!='passed':raise AssertionError(result)
        cls.output=Path(os.environ.get('VELA_REFERENCE_SDEVICE_TEST_OUTPUT',cls.root/'native checks'))
        cls.specs=checks.load_cases()['cases']
        cls.results=checks.run_suite(cls.meshes,cls.output,cls.launcher)

    def test_all_eight_real_commands_pass_check_only(self):
        self.assertEqual('passed',self.results['status'],self.results)
        self.assertEqual(8,len(self.results['cases']))
        for row in self.results['cases']:
            self.assertEqual('checked',row['status'])
            self.assertEqual('not_run',row['numerical_qualification'])

    def reject(self,label,mutate,identifier,exit_code):
        directory=self.root/label
        case=self.specs[0]
        checks.prepare_case(case,self.meshes,directory)
        mutate(directory)
        result=checks.invoke_check(self.launcher,directory)
        self.assertEqual(exit_code,result.returncode,result.stderr)
        work=directory/'device_vela'
        failure=json.loads((work/'fail_report.json').read_text())
        self.assertEqual(identifier,failure['message_id'])
        run=json.loads((work/'run_manifest.json').read_text())
        self.assertEqual([],run['stages'])
        self.assertFalse((directory/'forbidden_tool.called').exists())

    def test_missing_parameter_declaration_E301(self):
        def change(d):
            p=d/'device.cmd';p.write_text(re.sub(r'Parameter="models.par"','',p.read_text()))
        self.reject('missing declaration',change,'E301',3)

    def test_missing_parameter_file_E303(self):
        self.reject('missing parameter',lambda d:(d/'models.par').unlink(),'E303',3)

    def test_wrong_contact_name_E401(self):
        def change(d):
            p=d/'device.cmd';p.write_text(p.read_text().replace('Name="Anode"','Name="missing_anode"'))
        self.reject('wrong contact',change,'E401',4)

    def test_unknown_volume_material_E403(self):
        def change(d):
            p=d/'device.grd';s,n=re.subn(r'material\s*=\s*Silicon','material = UnknownMaterial',p.read_text())
            self.assertGreater(n,0);p.write_text(s)
        self.reject('unknown material',change,'E403',4)

    def test_unsupported_physics_E106(self):
        def change(d):
            p=d/'device.cmd';p.write_text(p.read_text().replace('Physics {','Physics { Thermodynamic',1))
        self.reject('unsupported model',change,'E106',1)

    def test_sidecar_cannot_replace_physics_E409(self):
        def change(d):
            (d/'device.vela.json').write_text(json.dumps({'schema':'vela.sdevice_sidecar.v1',
                                                       'solver':{'mobility':{'model':'constant'}}}))
        self.reject('forbidden sidecar',change,'E409',4)

    def test_emitted_goal_corruption_is_detected(self):
        case=next(c for c in self.specs if c['id']=='pn2d_iv')
        directory=self.root/'changed stage goal'
        shutil.copytree(self.output/case['id'],directory)
        row=next(r for r in self.results['cases'] if r['id']==case['id'])
        packet=copy.deepcopy(row['inputs'])
        packet['original_ir']=checks.parse_cmd_text((REPO/case['command']['path']).read_text(),path_name='device.cmd')
        p=directory/'device_vela/simulation_stage01.json';data=json.loads(p.read_text())
        data['sweep']['stop']=9.0;p.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError,'sweep goal changed'):
            checks.verify_check(case,packet,self.meshes,directory)

    def test_modified_p1_grid_is_rejected_before_native_run(self):
        altered=self.root/'altered meshes';shutil.copytree(self.meshes,altered)
        with (altered/'pn2d/device.grd').open('a') as f:f.write('\nchanged\n')
        output=self.root/'should not start'
        with self.assertRaisesRegex(ValueError,'generated file modified'):
            checks.run_suite(altered,output,self.launcher)
        self.assertFalse(output.exists())

    def test_existing_output_is_preserved(self):
        old=(self.output/'check_summary.json').read_bytes()
        with self.assertRaisesRegex(ValueError,'fresh empty output'):
            checks.run_suite(self.meshes,self.output,self.launcher)
        self.assertEqual(old,(self.output/'check_summary.json').read_bytes())


if __name__=='__main__':unittest.main()
