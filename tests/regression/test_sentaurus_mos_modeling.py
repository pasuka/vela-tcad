"""P4 geometry, shortcut semantics, source safety, and complete-device gates."""
from pathlib import Path
import sys,tempfile,json,math,unittest,copy,os
REPO=Path(__file__).resolve().parents[2];sys.path.insert(0,str(REPO/'scripts'))
from sentaurus_device_ir import parse_sde_device_ir
from sentaurus_ir_v2 import migrate_device_ir_v1
from sentaurus_curved import polygon,area,contains
from jsonschema import validate

def parse(text,parameters=None):
    with tempfile.TemporaryDirectory() as t:
        p=Path(t)/'input.cmd';p.write_text(text);return parse_sde_device_ir(p,parameters=parameters)

class GeometryTests(unittest.TestCase):
    def test_all_convex_corners_meet_circle_chord_and_area_bounds(self):
        for corner in [[0,0],[1,0],[1,1],[0,1]]:
            b={'lower_left':[0,0],'upper_right':[1,1],'fillets':[{'corner':corner,'radius_um':.05}]}
            exact=1-.05**2*(1-math.pi/4)
            self.assertLessEqual(abs(area(b)-exact),math.pi*.05/2*1e-5)
            self.assertFalse(contains(b,corner))
            self.assertTrue(contains(b,[.5,.5]))

    def test_originals_preserve_fillet_and_mirror_in_migration(self):
        for family,name in [('singledevice_sentaurus2018','singledevice_sde.cmd'),('transportmodels_sentaurus2022','nMOS_dvs.cmd')]:
            d=parse_sde_device_ir(REPO/'reference_tcad'/family/'source'/name,parameters={'node':'7'})
            v=migrate_device_ir_v1(d)
            self.assertEqual(.05,next(b for b in v['geometry']['bodies'] if b['region']=='R.Spacer')['geometry']['fillets'][0]['radius_um'])
            self.assertEqual('n7_half_msh',v['mesh_control']['operations'][0]['input'])
            for ver,ir in [(1,d),(2,v)]:validate(ir,json.loads((REPO/f'schemas/vela.sentaurus_device_ir.v{ver}.schema.json').read_text()))

    def test_unknown_shell_and_unsafe_parameter_are_rejected(self):
        for command in ['system:command','system']:
            with self.assertRaises(ValueError):parse(f'({command} "touch should_not_exist")')
        with self.assertRaisesRegex(ValueError,'unresolved'):parse('(sde:build-mesh "n@node@")')
        with self.assertRaisesRegex(ValueError,'unsafe'):parse('(sde:build-mesh "n@node@")',{'node':'0; rm something'})
        with self.assertRaises(ValueError):parse('(sde:build-mesh "half") (system:command "tdx -mtt -y -M 0 -S 0 -ren drain=source other full")')

    def test_oversized_or_ambiguous_corner_rejects(self):
        text='(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R")'
        for p,r in [('0 0 0',2),('.5 .5 0',.1)]:
            with self.assertRaises(ValueError):parse(text+f'(sdegeo:fillet-2d (find-vertex-id (position {p})) {r})')

class ShortcutTests(unittest.TestCase):
    def test_documented_defaults_and_directional_coordinates(self):
        d=parse('''(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R")
        (sdedr:refine-box "W" (position 0 0 0) (position 1 1 0) .1 .2 .3)
        (sdedr:refine-doping .001)
        (sdedr:refine-interface "material" "Silicon" "Oxide" "hlocal" .001)''')
        sizes=d['mesh_control']['refinement_sizes']
        box=next(s for s in sizes if s['name']=='Shortcut.W')
        self.assertEqual((.1,.2),(box['max_x'],box['max_y']))
        functions=d['mesh_control']['refinement_functions']
        self.assertEqual([1.0],functions[0]['arguments'])
        self.assertEqual(2.,functions[1]['factor']);self.assertFalse(functions[1]['double_side'])
        place=d['mesh_control']['refinement_placements'][1]
        self.assertEqual(('material','Silicon'),(place['target_kind'],place['target']))

class AcceptedMeshTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        value=os.environ.get('VELA_P4_MESHES')
        if not value:raise unittest.SkipTest('set VELA_P4_MESHES to independently validated generated corpus')
        cls.root=Path(value)

    def test_full_originals_and_corruption_detection(self):
        from sentaurus_mos_contract import check
        for case in ['singledevice','transportmodels']:
            b=json.loads((self.root/case/'mesh_field_bundle.json').read_text())
            self.assertEqual('passed',check(b)['status'])
            bad=copy.deepcopy(b);bad['fields']['DonorConcentration']['values'][-1]+=1e20
            self.assertEqual('failed',check(bad)['status'])
            bad=copy.deepcopy(b);bad['nodes'][-1]['y']+=.001
            self.assertEqual('failed',check(bad)['status'])

    def test_half_and_full_native_imports_match_authoritative_fields(self):
        from vela_sdevice.mesh_import import import_grid
        from sentaurus_sde_cases import compare_native
        from sentaurus_mesh_exporters import write_doping_csv,bundle_to_vela_mesh
        import shutil
        for case in ['singledevice','transportmodels']:
            out=self.root/case
            with tempfile.TemporaryDirectory() as t:
                import_grid(out/'device.grd',Path(t));compare_native(out,Path(t))
            with tempfile.TemporaryDirectory() as t:
                base=Path(t);stage=base/'half';stage.mkdir()
                b=json.loads((out/'half_mesh_field_bundle.json').read_text())
                (stage/'mesh_field_bundle.json').write_text(json.dumps(b))
                (stage/'mesh.json').write_text(json.dumps(bundle_to_vela_mesh(b)))
                for ext in ['grd','dat']:shutil.copy2(out/('half_device.'+ext),stage/('device.'+ext))
                write_doping_csv(stage/'doping.csv',b)
                import_grid(stage/'device.grd',base/'native');compare_native(stage,base/'native')

class ConsumerBoundaryTests(unittest.TestCase):
    def test_original_physics_commands_reject_without_solver_execution(self):
        import shutil,re,hashlib
        from check_sde_sdevice_cases import invoke_check
        meshes=os.environ.get('VELA_P4_MESHES')
        launcher=Path(os.environ.get('VELA_SD_SIM',REPO/'build-release/sd-sim.exe'))
        if not meshes or not launcher.is_file():self.skipTest('requires generated P4 meshes and sd-sim')
        config=json.loads((REPO/'configs/reference_mos_acceptance.json').read_text())
        with tempfile.TemporaryDirectory() as temp:
            out=Path(os.environ.get('VELA_P4_MOS_CHECK_OUTPUT',temp))
            out.mkdir(parents=True,exist_ok=True)
            self.assertFalse(any(out.iterdir()),'use fresh consumer output directory')
            reports=[]
            for case in config['consumer_boundaries']:
                source=REPO/case['source'];directory=out/case['id'];directory.mkdir()
                self.assertEqual(case['source_sha256_lf'],hashlib.sha256(source.read_bytes().replace(b'\r\n',b'\n')).hexdigest())
                text,n=re.subn(r'(\bGrid\s*=\s*)"[^"]+"',r'\1"device.grd"',source.read_text());self.assertEqual(1,n)
                (directory/'device.cmd').write_text(text)
                for p in source.parent.glob('*.par'):shutil.copy2(p,directory/p.name)
                for extension in ['grd','dat']:shutil.copy2(Path(meshes)/case['id']/('device.'+extension),directory/('device.'+extension))
                result=invoke_check(launcher,directory)
                (directory/'stdout.log').write_text(result.stdout);(directory/'stderr.log').write_text(result.stderr)
                self.assertEqual(case['exit_code'],result.returncode,result.stderr)
                report=json.loads((directory/'device_vela/fail_report.json').read_text())
                ids={report['message_id']}|{f['id'] for f in report.get('details',{}).get('findings',[])}
                self.assertTrue(set(case['diagnostic_ids'])<=ids)
                reports.append({'id':case['id'],'status':'expected_rejection','exit_code':result.returncode,'diagnostic_ids':sorted(ids),
                                'grid_sha256':hashlib.sha256((directory/'device.grd').read_bytes()).hexdigest(),
                                'data_sha256':hashlib.sha256((directory/'device.dat').read_bytes()).hexdigest(),
                                'source_sha256_lf':case['source_sha256_lf'],'effective_command_sha256':hashlib.sha256((directory/'device.cmd').read_bytes()).hexdigest()})
            (out/'consumer_boundary_summary.json').write_text(json.dumps({'status':'passed','mode':'check_only','solver_invoked':False,'cases':reports},indent=2))

if __name__=='__main__':unittest.main()
