"""P6 offline packaging boundaries and extracted-runtime behavior."""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys,tempfile,unittest,zipfile
from unittest.mock import patch

REPO=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'scripts'))
from package_sde_modeling import install_wheels,module_closure,verify_inputs
from prepare_sde_offline_assets import download


class InputBoundaryTests(unittest.TestCase):
    def test_corrupted_existing_download_is_preserved(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);target=root/'asset.whl';target.write_bytes(b'owned data')
            with self.assertRaisesRegex(ValueError,'hash mismatch'):
                download({'filename':target.name,'sha256':'0'*64,'url':'https://invalid.example'},root)
            self.assertEqual(b'owned data',target.read_bytes())

    def test_wheel_path_traversal_fails_before_writing_outside_vendor(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);wheel=root/'bad.whl'
            with zipfile.ZipFile(wheel,'w') as z:z.writestr('../escape.py','bad')
            with self.assertRaisesRegex(ValueError,'Unsafe archive'):
                install_wheels([wheel],root/'vendor')
            self.assertFalse((root/'escape.py').exists())

    def test_wheel_collisions_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);wheels=[]
            for i in (0,1):
                path=root/f'{i}.whl';wheels.append(path)
                with zipfile.ZipFile(path,'w') as z:z.writestr('same.py',str(i))
            with self.assertRaisesRegex(ValueError,'collision'):install_wheels(wheels,root/'vendor')
            self.assertEqual('0',(root/'vendor/same.py').read_text())

    def test_complete_offline_inputs_are_required_before_unpacking(self):
        lock=json.loads((REPO/'configs/sde_offline_inputs.json').read_text())
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(ValueError,'offline input missing'):
                verify_inputs(lock,Path(temp),'win64')

    def test_deferred_local_imports_are_included_but_unrelated_solver_tools_are_not(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);(root/'scripts').mkdir()
            for name,text in [('sentaurus_gmsh_pipeline','def f():\n import late\n'),('late','from helper import f\n'),('helper','f=1\n'),('unrelated','raise RuntimeError()\n')]:
                (root/'scripts'/(name+'.py')).write_text(text)
            self.assertEqual(['helper','late','sentaurus_gmsh_pipeline'],module_closure(root))

    def test_ogawa_boundary_is_explicit_before_external_import(self):
        from vela_sdevice.mesh_import import import_grid,MeshImportError
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);tdr=root/'input.tdr';tdr.write_bytes(b'Ogawa'+bytes(20))
            with patch('vela_sdevice.mesh_import.find_importer') as external:
                with self.assertRaisesRegex(MeshImportError,'Ogawa TDR is unsupported'):
                    import_grid(tdr,root/'out')
                external.assert_not_called()

    def test_both_platforms_pin_offline_dependencies_and_external_gmsh(self):
        lock=json.loads((REPO/'configs/sde_offline_inputs.json').read_text())
        self.assertFalse(lock['gmsh']['bundled'])
        for platform in ('win64','ubuntu26.04-x86_64'):
            self.assertEqual({'jsonschema','attrs','jsonschema-specifications','referencing','rpds-py','typing-extensions'},
                             {p['name'] for p in lock['dependencies'][platform]})
            for item in lock['dependencies'][platform]:self.assertRegex(item['sha256'],r'^[0-9a-f]{64}$')


class ExtractedRuntimeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        value=os.environ.get('VELA_P6_PACKAGE')
        if not value:raise unittest.SkipTest('requires an extracted P6 package')
        cls.source=Path(value)
        cls.temp=tempfile.TemporaryDirectory(prefix='p6-package-')
        cls.addClassCleanup(cls.temp.cleanup)
        cls.base=Path(cls.temp.name)/'moved space \u4e2d\u6587';cls.base.mkdir()
        cls.root=cls.base/'vela-sde';shutil.copytree(cls.source,cls.root,symlinks=True)
        cls.python=cls.root/('bin/python/python.exe' if os.name=='nt' else 'bin/python/bin/python3')

    def run_entry(self,*args):
        env={k:v for k,v in os.environ.items() if not k.startswith(('PYTHON','VELA_GMSH'))}
        env['PYTHONPATH']=str(self.base/'poison')
        return subprocess.run([str(self.python),'-B','-I','-S',str(self.root/'run_sde.py'),*args],
                              cwd=self.base,env=env,capture_output=True,text=True)

    def test_relocation_and_explicit_external_dependency_profile(self):
        result=self.run_entry('--version');self.assertEqual(0,result.returncode,result.stderr)
        manifest=json.loads(result.stdout);self.assertEqual('explicit external 4.15.2 (not bundled)',manifest['gmsh_profile'])
        self.assertTrue(manifest['python_version'].startswith('3.12.'))
        self.assertRegex(manifest['source_commit'],r'^[0-9a-f]{40}$')

    def test_absent_external_gmsh_fails_closed_without_output(self):
        output=self.base/'no-runtime-output'
        result=self.run_entry('--sde',str(self.base/'input.cmd'),'--output-dir',str(output))
        self.assertEqual(2,result.returncode,result.stderr)
        self.assertIn('explicit Gmsh',result.stderr);self.assertFalse(output.exists())

    def test_package_corruption_is_detected_before_execution(self):
        path=self.root/'scripts/sentaurus_device_ir.py';original=path.read_bytes()
        try:
            path.write_bytes(original+b'\n# changed\n')
            result=self.run_entry('--version');self.assertEqual(2,result.returncode)
            self.assertIn('hash mismatch',result.stderr)
        finally:path.write_bytes(original)

    def test_no_reference_outputs_or_gmsh_binary_are_distributed(self):
        manifest=json.loads((self.root/'PACKAGE_MANIFEST.json').read_text())
        for item in manifest['files']:
            path=Path(item['path'])
            self.assertNotIn(path.suffix.lower(),{'.tdr','.grd','.dat','.msh','.log','.cmd'})
            self.assertNotIn('reference_tcad',path.parts)
            self.assertNotIn(path.name.lower(),{'gmsh.py','gmsh.exe','libgmsh.dll','gmsh-4.15.dll','libgmsh.so.4.15'})


class QualificationFailureTests(unittest.TestCase):
    def run_failed_child(self, payload, expected_message):
        import qualify_sde_package as qualification
        with tempfile.TemporaryDirectory(prefix='sde-qualification-failure-') as temp:
            base=Path(temp);package=base/'package';output=base/'output';output.mkdir()
            interpreter=package/('bin/python/python.exe' if os.name=='nt' else 'bin/python/bin/python3')
            interpreter.parent.mkdir(parents=True);interpreter.touch()
            # A failed rerun must not leave a previous success receipt behind.
            (output/'qualification.json').write_text(json.dumps({'status':'passed'}))
            case={'id':'bvmethods','source':'unused.cmd','parameters':{},'contract':{}}
            def fail(command, **kwargs):
                folder=Path(command[command.index('--output-dir')+1]);folder.mkdir(parents=True)
                if payload is not None:
                    text=payload if isinstance(payload,str) else json.dumps(payload)
                    (folder/'fail_report.json').write_text(text,encoding='utf-8')
                return subprocess.CompletedProcess(command,2,'child stdout\n','Gmsh preprocessing failed; inspect report\n')
            with patch('sentaurus_sde_cases.load_manifest',return_value={'cases':[case]}), \
                 patch('qualify_sde_package.subprocess.run',side_effect=fail), \
                 patch('sentaurus_sde_cases.check_contract') as contract:
                with self.assertRaisesRegex(RuntimeError,expected_message):
                    qualification.qualify(package,output,[base/'gmsh.py',base/'libgmsh',base/'gmsh'],1)
                contract.assert_not_called()
            report=json.loads((output/'qualification.json').read_text())
            self.assertEqual('failed',report['status']);self.assertEqual([],report['cases'])
            self.assertEqual(2,report['failure']['exit_code'])
            self.assertEqual('bvmethods',report['failure']['id'])
            self.assertEqual(0,report['failure']['iteration'])
            self.assertEqual([],report['failure']['completed_runs'])
            self.assertIn('child stdout',(output/'bvmethods-0.log').read_text())
            return report['failure']

    def test_modeling_failure_surfaces_inner_stage_and_preserves_failure_status(self):
        failure=self.run_failed_child({
            'schema':'vela.sentaurus_gmsh_preprocess.fail_report.v1','status':'failed',
            'stage':'gmsh_backend','message':'MaxTransDiff iteration budget exhausted',
            'details':{'max_iterations':12}},r'\[gmsh_backend\].*iteration budget exhausted')
        self.assertEqual('gmsh_backend',failure['stage'])
        self.assertEqual('fail_report.json',failure['report'])

    def test_missing_or_malformed_report_keeps_original_child_failure(self):
        for payload in (None,'{not JSON',{'status':'failed','stage':'unowned','message':'wrong schema'}):
            with self.subTest(payload=payload):
                failure=self.run_failed_child(payload,'Gmsh preprocessing failed')
                self.assertEqual('subprocess',failure['stage'])
                self.assertNotIn('report',failure)


if __name__=='__main__':unittest.main()
