"""P5 lifecycle, signed grading, Erf endpoints, and original BV contracts."""
from pathlib import Path
import copy, json, math, os, subprocess, sys, tempfile, unittest
R=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(R/'scripts'))
from sentaurus_device_ir import parse_sde_device_ir
from sentaurus_ir_v2 import migrate_device_ir_v1
from sentaurus_doping import evaluate_m2_doping_at
from sentaurus_multibox import spacing, validate as validate_multibox
from sentaurus_gmsh_backend import GmshBackendError
from jsonschema import validate

RECT='(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R.Base")'
GATE='(sdegeo:create-rectangle (position .2 .7 0) (position .8 1.2 0) "PolySi" "R.Gate")'
REMOVE='(sdegeo:set-contact (find-body-id (position .5 .9 0)) "gate" "remove")'
DOPE='''(sdedr:define-constant-profile "C" "ArsenicActiveConcentration" 1e20)
(sdedr:define-constant-profile-region "P" "C" "R.Gate")'''
BV=R/'reference_tcad/bvmethods_sentaurus2018/source/bvmethods_nmos_sde.cmd'


def parse(text):
    with tempfile.TemporaryDirectory() as t:
        p=Path(t)/'test.cmd';p.write_text(text)
        return migrate_device_ir_v1(parse_sde_device_ir(p))


class LifecycleTests(unittest.TestCase):
    def test_alias_identity_and_deleted_profile_stay_in_ir(self):
        ir=parse(RECT+GATE+REMOVE+DOPE)
        body=ir['geometry']['bodies'][1]
        self.assertEqual({'sde':'PolySi','vela':'PolySilicon'},body['material'])
        self.assertIn('removed',body)
        self.assertEqual('body_boundary',ir['contacts'][0]['assignments'][0]['kind'])
        self.assertEqual('R.Gate',ir['doping']['placements'][0]['target'])
        self.assertEqual((0.,0.,{}),evaluate_m2_doping_at(ir,.5,.9))
        validate(ir,json.loads((R/'schemas/vela.sentaurus_device_ir.v2.schema.json').read_text()))

    def test_removed_overlap_does_not_resurrect_base(self):
        from sentaurus_gmsh_backend import _point_owner
        ir=parse(RECT+GATE+REMOVE)
        self.assertIsNone(_point_owner(ir,.5,.9))
        self.assertEqual('R.Base',_point_owner(ir,.1,.9)['region'])

    def test_bab_created_after_remove_fills_vacancy(self):
        from sentaurus_gmsh_backend import _point_owner
        ir=parse(RECT+GATE+REMOVE+'(sdegeo:set-default-boolean "BAB")'+
                 '(sdegeo:create-rectangle (position .3 .8 0) (position .7 1.1 0) "Oxide" "R.Fill")')
        self.assertEqual('R.Fill',_point_owner(ir,.5,.9)['region'])

    def test_partial_or_prior_hidden_bab_removal_fails_closed(self):
        for boolean in ['ABA','BAB']:
            with self.assertRaisesRegex(ValueError,'later overlapping'):
                parse(RECT+GATE+f'(sdegeo:set-default-boolean "{boolean}")'+
                      '(sdegeo:create-rectangle (position .3 .8 0) (position .7 1.1 0) "Oxide" "R.Other")'+REMOVE.replace('.5 .9','.25 .9'))
        with self.assertRaisesRegex(ValueError,'interior'):
            parse(RECT+GATE+REMOVE.replace('.5 .9','.2 .9'))

    def test_bv_live_fillet_and_x_mirror_survive_both_schemas(self):
        v1=parse_sde_device_ir(BV,parameters={'node':'0'});v2=migrate_device_ir_v1(v1)
        spacer=next(b for b in v2['geometry']['bodies'] if b['region']=='R.Spacer')
        self.assertEqual(.08,spacer['geometry']['fillets'][0]['radius_um'])
        self.assertEqual('x',v2['mesh_control']['operations'][0]['axis'])
        for v,d in [(1,v1),(2,v2)]:
            validate(d,json.loads((R/f'schemas/vela.sentaurus_device_ir.v{v}.schema.json').read_text()))


class GradingTests(unittest.TestCase):
    def setUp(self):
        self.size=dict(min_x=.01,max_x=.08,min_y=.01,max_y=.08,ratio_x=1.4,ratio_y=-1.4)
        self.window=dict(lower_left=[0.,0.],upper_right=[1.,1.])

    def test_geometric_layers_have_opposite_signed_directions(self):
        for n in range(6):
            d=.01*(1.4**n-1)/.4
            hx,hy=spacing(self.size,self.window,d,1-d,1.)
            self.assertAlmostEqual(min(.08,.01*1.4**n),hx)
            self.assertAlmostEqual(hx,hy)

    def test_ratio_one_and_zero_are_max_and_disabled(self):
        self.size.update(ratio_x=1.,ratio_y=0.)
        self.assertEqual([.08,1.],spacing(self.size,self.window,.5,.5,1.))
        self.size['ratio_x']=-1.
        self.assertEqual(.08,spacing(self.size,self.window,.5,.5,1.)[0])

    def test_malformed_ir_does_not_bypass_parser(self):
        ir=migrate_device_ir_v1(parse_sde_device_ir(BV,parameters={'node':'0'}))
        for key,value in [('ratio_x',.5),('min_y',0.),('max_y',float('nan')),('min_x',2.)]:
            bad=copy.deepcopy(ir);bad['mesh_control']['multibox_sizes'][0][key]=value
            with self.assertRaises(GmshBackendError):validate_multibox(bad)

    def test_overlapping_windows_enforce_finest_and_clip_outside(self):
        from sentaurus_multibox import marked_edges
        s=dict(self.size,name='S',ratio_x=1.,ratio_y=0.)
        w=dict(self.window,name='W',shape='rectangle')
        ir={'mesh_control':{'multibox_sizes':[s],'windows':[w],
                          'multibox_placements':[dict(name='P',multibox='S',window='W')]}}
        points={0:(.2,.2),1:(.4,.2),2:(.2,.4),3:(2.,2.),4:(3.,2.),5:(2.,3.)}
        marks=marked_edges(ir,points,[[0,1,2],[3,4,5]],1.)
        self.assertEqual({(0,1),(1,2)},marks)
        smaller=dict(s,name='T',max_x=.02)
        ir['mesh_control']['multibox_sizes'].append(smaller)
        ir['mesh_control']['multibox_placements'].append(dict(name='Q',multibox='T',window='W'))
        self.assertEqual(marks,marked_edges(ir,points,[[0,1,2],[3,4,5]],1.))


class ErfTests(unittest.TestCase):
    def test_continuous_finite_baseline_endpoint_and_symm(self):
        source='''(sdegeo:create-rectangle (position -.1 -.1 0) (position 1.1 .1 0) "Silicon" "R")
        (sdedr:define-refeval-window "B" "Line" (position 0 0 0) (position 1 0 0))
        (sdedr:define-gaussian-profile "G" "ArsenicActiveConcentration" "PeakPos" 0 "PeakVal" 1e18 "ValueAtDepth" 1e16 "Depth" .1 "Erf" "Factor" .4)
        (sdedr:define-analytical-profile-placement "P" "G" "B" "Symm" "NoReplace" "Eval")'''
        ir=parse(source)
        self.assertEqual('Both',ir['doping']['placements'][0]['direction'])
        self.assertAlmostEqual(.5,evaluate_m2_doping_at(ir,0.,0.)[0]/1e18)
        left=evaluate_m2_doping_at(ir,-1e-10,0.)[0]
        right=evaluate_m2_doping_at(ir,1e-10,0.)[0]
        self.assertLess(abs(right-left)/1e18,1e-7)
        for x,y in [(-.025,.02),(.025,.05),(.5,.1),(.99,.05)]:
            self.assertEqual(evaluate_m2_doping_at(ir,x,y),evaluate_m2_doping_at(ir,x,-y))


class MultiboxPipelineTests(unittest.TestCase):
    def test_live_polysi_alias_survives_native_export_import(self):
        from tests.regression.gmsh_test_runtime import require_gmsh_runtime
        require_gmsh_runtime()
        from sentaurus_sde_cases import compare_native
        from vela_sdevice.mesh_import import import_grid
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);source=root/'alias.cmd'
            source.write_text('''(sdegeo:create-rectangle (position 0 0 0) (position .32 .32 0) "Silicon" "R.Si")
            (sdegeo:create-rectangle (position 0 .32 0) (position .32 .4 0) "PolySi" "R.Gate")
            (sdedr:define-constant-profile "C" "ArsenicActiveConcentration" 1e20)
            (sdedr:define-constant-profile-region "P" "C" "R.Gate")
            (sde:build-mesh "alias")''')
            out=root/'out'
            result=subprocess.run([sys.executable,str(R/'scripts/sentaurus_gmsh_pipeline.py'),
                '--sde',str(source),'--output-dir',str(out),'--runtime-dir',str(root/'runtime')],
                cwd=R,capture_output=True,text=True)
            self.assertEqual(0,result.returncode,result.stderr)
            import_grid(out/'device.grd',root/'native');compare_native(out,root/'native')
            mesh=json.loads((root/'native/mesh.json').read_text())
            gate=next(r for r in mesh['regions'] if r['name']=='R.Gate')
            self.assertEqual('PolySilicon',gate['material'])
            ir=json.loads((out/'device_ir.json').read_text())
            authored=next(b for b in ir['geometry']['bodies'] if b['region']=='R.Gate')
            self.assertEqual('PolySi',authored['material']['sde'])

    def test_multibox_only_positive_negative_and_overlap_export(self):
        from tests.regression.gmsh_test_runtime import require_gmsh_runtime
        require_gmsh_runtime()
        from sentaurus_sde_cases import compare_native
        from vela_sdevice.mesh_import import import_grid
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for sign in [1,-1]:
                src=root/f'case{sign}.cmd'
                src.write_text('''(sdegeo:create-rectangle (position 0 0 0) (position .32 .32 0) "Silicon" "R")
                (sdedr:define-refeval-window "W" "Rectangle" (position -.1 0 0) (position .4 .32 0))
                (sdedr:define-multibox-size "M" .08 .08 .01 .01 1 RATIO)
                (sdedr:define-multibox-placement "P" "M" "W")
                (sdedr:define-multibox-size "N" .04 .06 .02 .01 1 0)
                (sdedr:define-multibox-placement "Q" "N" "W")
                (sde:build-mesh "mb")'''.replace('RATIO',str(sign*1.4)))
                out=root/f'out{sign}'
                p=subprocess.run([sys.executable,str(R/'scripts/sentaurus_gmsh_pipeline.py'),'--sde',str(src),'--output-dir',str(out),'--runtime-dir',str(root/'runtime')],capture_output=True,text=True,cwd=R)
                self.assertEqual(0,p.returncode,p.stderr)
                q=json.loads((out/'mesh_qualification.json').read_text())
                self.assertEqual('passed',q['refinement']['status'])
                self.assertEqual(2,len(q['refinement']['multibox']))
                import_grid(out/'device.grd',root/f'native{sign}');compare_native(out,root/f'native{sign}')


class CompleteBVTests(unittest.TestCase):
    def test_complete_device_and_mutation_detection(self):
        from sentaurus_bv_contract import check
        value=os.environ.get('VELA_P5_MESHES')
        if not value:self.skipTest('requires generated corpus')
        out=Path(value)/'bvmethods'
        b=json.loads((out/'mesh_field_bundle.json').read_text())
        self.assertEqual('passed',check(b)['status'])
        summary=json.loads((Path(value)/'case_summary.json').read_text())
        row=next(c for c in summary['cases'] if c['id']=='bvmethods')
        self.assertEqual(2,row['repeats'])
        self.assertTrue(all(len(run['half_sha256'])==4 for run in row['runs']))
        for mutation in ['doping','axis','geometry','adjacency','inactive']:
            bad=copy.deepcopy(b)
            if mutation=='doping':bad['fields']['DonorConcentration']['values'][-1]+=1e20
            if mutation=='axis':bad['provenance']['mirror']['axis']='y'
            if mutation=='geometry':bad['nodes'][-1]['x']+=.001
            if mutation=='adjacency':bad['provenance']['contact_owners']['gate']=['R.Gateox']
            if mutation=='inactive':next(n for n in bad['provenance']['refinement_field_dag']['nodes'] if n['id']=='multibox:MBPlace.Gate')['kind']='MathEvalAniso'
            self.assertEqual('failed',check(bad)['status'],mutation)


if __name__=='__main__':unittest.main()
