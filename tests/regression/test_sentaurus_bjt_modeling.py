from __future__ import annotations
import copy
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

REPO=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'scripts'))
from jsonschema import validate, ValidationError
from sentaurus_device_ir import parse_sde_device_ir, SdeParseError
from sentaurus_ir_v2 import migrate_device_ir_v1
from sentaurus_general import expression_tree, initialise, evaluate, validate_profile
from sentaurus_doping import evaluate_m2_doping_at
from sentaurus_bjt_contract import concentration, integrals, exact_roots
from sentaurus_adaptive import analyse, validate as validate_adaptive
from sentaurus_gmsh_backend import GmshBackendError

BJT=REPO/'reference_tcad/genius_bjt_sentaurus2022/source/bjt_sde.cmd'

def parse(text):
    with tempfile.TemporaryDirectory() as tmp:
        p=Path(tmp)/'case.cmd';p.write_text(text)
        return parse_sde_device_ir(p)

def general(expression,init=''):
    params=initialise(init,{})
    return dict(expression=expression,parameters=params,expression_tree=expression_tree(expression,set(params)|set('xyz')),default_value_cm3=0.0)

class TypedDefineTests(unittest.TestCase):
    def test_body_handles_aliases_and_creation_order_survive_migration(self):
        v1=parse('''(define a (sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "A"))
        (define alias a)
        (define b (sdegeo:create-rectangle (position 1 0 0) (position 2 1 0) "Oxide" "B"))''')
        v2=migrate_device_ir_v1(v1)
        self.assertEqual(v1['defines'],v2['defines'])
        self.assertEqual(v1['defines']['a'],v1['defines']['alias'])
        self.assertEqual(['body:0','body:1'],[b['id'] for b in v2['geometry']['bodies']])
        for v,d in [(1,v1),(2,v2)]:
            validate(d,json.loads((REPO/f'schemas/vela.sentaurus_device_ir.v{v}.schema.json').read_text()))

    def test_nested_bounded_numeric_and_position_helpers(self):
        d=parse('''(define scale 2) (define (twice x) (* scale x))
        (define (point x y) (position (twice x) y 0))
        (sdegeo:create-rectangle (point 0 0) (point 1 1) "Silicon" "R")''')
        self.assertEqual([2,1],d['geometry']['regions'][0]['upper_right'])
        self.assertEqual('pure_function',d['defines']['point']['kind'])
        self.assertEqual(2,d['defines']['scale'])

    def test_impure_body_and_unbounded_recursion_rejected(self):
        with self.assertRaisesRegex(SdeParseError,'impure/unsupported'):
            parse('(define (bad x) (system x))')
        with self.assertRaisesRegex(ValueError,'depth exceeds 16'):
            parse('(define (loop x) (loop x)) (define n (loop 1))')

    def test_total_call_budget_bounds_nonrecursive_exponential_expansion(self):
        forms=['(define (f0 x) (+ x 1))']
        for i in range(1,14):
            forms.append(f'(define (f{i} x) (+ (f{i-1} x) (f{i-1} x)))')
        forms.append('(define n (f13 0))')
        with self.assertRaisesRegex(ValueError,'call budget exceeds 4096'):
            parse(' '.join(forms))

    def test_cross_kind_rebinding_is_rejected(self):
        for text in ['(define a 1) (define (a x) (+ x 1))',
                     '(define (a x) (+ x 1)) (define a 1)',
                     '(define (exp x) (+ x 1))']:
            with self.assertRaises(SdeParseError):parse(text)

    def test_handle_is_not_a_scalar(self):
        with self.assertRaisesRegex(ValueError,'not numeric'):
            parse('''(define a (sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R"))
            (define b (+ a 1))''')

class GeneralTests(unittest.TestCase):
    def test_global_initialization_first_assignment_wins(self):
        namespace={}
        self.assertEqual({'a':2.0},initialise('a=2',namespace))
        self.assertEqual({'a':2.0,'b':4.0},initialise('a=3;b=a^2',namespace))
        self.assertEqual((16.0,False),evaluate(general('b*(x+y)','b=4'),1,3))

    def test_ast_rejects_execution_and_resource_abuse(self):
        for text in ["__import__('os')",'x.real','x[0]','[x for x in y]','lambda: x','unknown+1','exp(x,y)','1e999', '(' * 1000+'x'+')'*1000]:
            with self.subTest(text=text[:30]),self.assertRaises(ValueError):
                expression_tree(text,{'x','y'})
        with self.assertRaises(ValueError):expression_tree('+'.join(['x']*200),{'x'})

    def test_defaults_are_explicit_and_negative_concentrations_reject(self):
        for expression in ['sqrt(-1)','1/0','exp(10000)','2^65']:
            p=general(expression);p['default_value_cm3']=7
            self.assertEqual((7.0,True),evaluate(p,0,0))
        with self.assertRaisesRegex(ValueError,'negative'):evaluate(general('-1'),0,0)

    def test_ast_tampering_is_detected(self):
        p=general('x+1');p['expression_tree']['op']='mul'
        with self.assertRaisesRegex(ValueError,'does not match'):validate_profile(p)

    def test_original_bjt_profiles_match_independent_max_distance_oracle(self):
        v1=parse_sde_device_ir(BJT);ir=migrate_device_ir_v1(v1)
        profiles=json.loads((REPO/'reference_tcad/genius_bjt_sentaurus2022/contracts/geometry_contract.json').read_text())['profiles']
        for x in [0,1.13,1.25,1.625,2.6,2.75,3.5,4.25,4.75,6]:
            for y in [0,.1,.17,.35,.42,.55,.7,1.45,1.73,2]:
                d,a,_=evaluate_m2_doping_at(ir,x,y);ed,ea=concentration(profiles,x,y)
                self.assertLessEqual(abs(d-ed)/max(ed,1),1e-12)
                self.assertLessEqual(abs(a-ea)/max(ea,1),1e-12)
        for v,d in [(1,v1),(2,ir)]:validate(d,json.loads((REPO/f'schemas/vela.sentaurus_device_ir.v{v}.schema.json').read_text()))

    def test_general_unsupported_placement_modes_fail_closed(self):
        from sentaurus_gmsh_backend import validate_m5_device_ir
        for key,value in [('direction','Positive'),('evaluation','NoEval'),('restriction_decay_length',.01)]:
            ir=migrate_device_ir_v1(parse_sde_device_ir(BJT))
            ir['doping']['placements'][1][key]=value
            with self.assertRaisesRegex(GmshBackendError,'General requires Both/Eval'):
                validate_m5_device_ir(ir)

    def test_general_schema_requires_ast(self):
        d=parse_sde_device_ir(BJT);del d['doping']['profiles'][1]['expression_tree']
        with self.assertRaises(ValidationError):validate(d,json.loads((REPO/'schemas/vela.sentaurus_device_ir.v1.schema.json').read_text()))

    def test_oracle_integral_and_junctions_have_independent_numerical_checks(self):
        profiles=[dict(species='donor',peak_cm3=3.,window_um=[0,0,1,1],char_um=None)]
        self.assertEqual({'donor':3.,'acceptor':0.},integrals(profiles,[0,0,1,1]))
        profiles=json.loads((REPO/'reference_tcad/genius_bjt_sentaurus2022/contracts/geometry_contract.json').read_text())['profiles']
        values=integrals(profiles,[0,0,6,2])
        # Independent midpoint quadrature at sufficient resolution for these smooth Gaussians.
        measured=[0.,0.]
        for i in range(600):
            for j in range(200):
                d,a=concentration(profiles,(i+.5)/100,(j+.5)/100)
                measured[0]+=d/10000;measured[1]+=a/10000
        for k,m in zip(('donor','acceptor'),measured):self.assertLess(abs(m-values[k])/values[k],1e-5)
        self.assertEqual([1,2,2],[len(exact_roots(profiles,x)) for x in [1.625,2.6,3.5]])

class AdaptiveTests(unittest.TestCase):
    def setUp(self):
        self.ir=migrate_device_ir_v1(parse_sde_device_ir(BJT))

    def test_signed_asinh_and_minimum_stop(self):
        # Opposite signs across a junction retain a large transformed defect.
        coords={1:(3.5,.2,0),2:(3.51,.2,0),3:(3.5,.201,0)}
        values={1:1e17,2:-1e17,3:1e17}
        targets,report=analyse(self.ir,coords,[[1,2,3]],values)
        self.assertGreater(report['max_transformed_difference'],70)
        self.assertTrue(targets)
        tiny={1:(3.5,.2,0),2:(3.501,.2,0),3:(3.5,.201,0)}
        targets,report=analyse(self.ir,tiny,[[1,2,3]],values)
        self.assertFalse(targets)
        self.assertEqual(2,report['minimum_stopped_edges'])

    def test_uniform_doping_does_not_refine(self):
        targets,report=analyse(self.ir,{1:(0,0,0),2:(1,0,0),3:(0,1,0)},[[1,2,3]],{1:1e17,2:1e17,3:1e17})
        self.assertFalse(targets);self.assertEqual(0,report['max_transformed_difference'])

    def test_unsupported_scopes_and_thresholds_fail_closed(self):
        function=self.ir['mesh_control']['refinement_functions'][0]
        function['arguments']=[float('nan')]
        with self.assertRaises(GmshBackendError):validate_adaptive(self.ir)
        function['arguments']=[.8]
        function['variable']='AcceptorConcentration'
        with self.assertRaises(GmshBackendError):validate_adaptive(self.ir)

class ExportPrecisionTests(unittest.TestCase):
    def test_dfise_roundtrip_preserves_cancellation_sensitive_net(self):
        from sentaurus_mesh_exporters import write_dfise, POINT_FIELD_NAMES
        from compare_sentaurus_tdr_tdx import parse_dat
        from vela_sdevice.mesh_import import import_grid
        import csv
        donor,acceptor=1.2345678901234567e19,1.2345678901234565e19
        net=donor-acceptor
        values=[donor,acceptor,net,abs(net),math.log10(abs(net))]
        bundle={"nodes":[dict(vela_id=i,dfise_id=i,x=x,y=y,z=0) for i,(x,y) in enumerate([(0,0),(1,0),(0,1)])],
                "triangles":[dict(node_ids=[0,1,2],dfise_id=0,region_id=0,material_id=0)],
                "lines":[],"regions":[dict(id=0,name="R")],"materials":[dict(id=0,name="Silicon")],
                "contacts":[],"interfaces":[],
                "fields":{**{name:{"values":[v]*3} for name,v in zip(POINT_FIELD_NAMES,values)},"species":{}}}
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);write_dfise(root/'device.grd',root/'device.dat',bundle)
            data=parse_dat(root/'device.dat')
            self.assertEqual(donor,data['DonorConcentration']['values'][0])
            self.assertEqual(acceptor,data['AcceptorConcentration']['values'][0])
            import_grid(root/'device.grd',root/'native')
            with (root/'native/doping.csv').open() as f:row=next(csv.DictReader(f))
            self.assertEqual(net,float(row['donors_cm3'])-float(row['acceptors_cm3']))

class LocalSplitTests(unittest.TestCase):
    def test_all_split_masks_preserve_area_and_shared_edge_conformity(self):
        from sentaurus_adaptive import split_triangles
        from collections import Counter
        coordinates={1:(0.,0.,0.),2:(1.,0.,0.),3:(0.,1.,0.),4:(1.,1.,0.)}
        triangles=[[1,2,3],[2,4,3]]
        edges=[(1,2),(2,3),(1,3),(2,4),(3,4)]
        for mask in range(32):
            marked={edge for i,edge in enumerate(edges) if mask & (1<<i)}
            points,result,mid=split_triangles(coordinates,triangles,marked)
            area=0.;counts=Counter()
            for tri in result:
                a,b,c=[points[t] for t in tri]
                signed=((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))/2
                self.assertGreater(signed,0)
                area+=signed
                counts.update(tuple(sorted((tri[i],tri[(i+1)%3]))) for i in range(3))
            self.assertAlmostEqual(1.,area)
            self.assertEqual(len(marked),len(mid))
            for (a,b),count in counts.items():
                pa,pb=points[a],points[b]
                boundary=(pa[0]==pb[0] and pa[0] in (0,1)) or (pa[1]==pb[1] and pa[1] in (0,1))
                self.assertEqual(1 if boundary else 2,count,(mask,a,b))

    def test_longest_edge_closure_preserves_locality_on_separated_mesh(self):
        from sentaurus_adaptive import longest_edge_closure
        points={1:(0,0,0),2:(2,0,0),3:(0,1,0),4:(10,0,0),5:(12,0,0),6:(10,1,0)}
        marked=longest_edge_closure(points,[[1,2,3],[4,5,6]],{(1,3)})
        self.assertEqual({(1,3),(2,3)},marked)

class BudgetTests(unittest.TestCase):
    def test_node_and_iteration_exhaustion_are_explicit_failures(self):
        from types import SimpleNamespace
        from unittest.mock import patch
        from sentaurus_adaptive import remesh
        ir=migrate_device_ir_v1(parse_sde_device_ir(BJT))
        coords=[1.625,0,0,1.625,1,0,2.0,0,0]
        gmsh=SimpleNamespace(logger=SimpleNamespace(write=lambda *args:None),model=SimpleNamespace(
            mesh=SimpleNamespace(getNodes=lambda:([1,2,3],coords,[])),
            getEntities=lambda dimension:[(2,1)]))
        with patch('sentaurus_gmsh_backend._triangles_for_surface',return_value=[(1,[1,2,3])]):
            with self.assertRaises(GmshBackendError) as context:
                remesh(gmsh,ir,{'root_field_tag':1},.2,max_iterations=0)
            self.assertEqual('iteration_budget',context.exception.details['reason'])
            validate(context.exception.details,json.loads((REPO/'schemas/vela.gmsh_maxtransdiff.v1.schema.json').read_text()))
            self.assertEqual(1,len(context.exception.details['rounds']))
        gmsh.model.mesh.getNodes=lambda:([1,2,3,4],coords+[2,1,0],[])
        with self.assertRaises(GmshBackendError) as context:
            remesh(gmsh,ir,{'root_field_tag':1},.2,max_nodes=3)
        self.assertEqual('node_budget',context.exception.details['reason'])

    def test_budget_failure_invalidates_previous_success_outputs(self):
        from tests.regression.gmsh_test_runtime import require_gmsh_runtime
        import sentaurus_gmsh_pipeline as pipeline
        require_gmsh_runtime()
        source_text='\n'.join([
            '(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R")',
            '(sdedr:define-refeval-window "W" "Line" (position 0 0 0) (position 1 0 0))',
            '(sdedr:define-analytical-profile "D" "PhosphorusActiveConcentration" "" "1e16+1e17*x" 0 "General")',
            '(sdedr:define-analytical-profile-placement "P" "D" "W" "Both" "NoReplace" "Eval")',
            '(sdedr:define-refinement-size "S" .2 .2 .02 .02)',
            '(sdedr:define-refinement-region "RP" "S" "R")'])
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);src=root/'case.cmd';out=root/'outputs'
            src.write_text(source_text)
            args=pipeline.build_parser().parse_args(['--sde',str(src),'--output-dir',str(out),'--runtime-dir',str(root/'runtime')])
            import subprocess
            def run():
                return subprocess.run([sys.executable,str(REPO/'scripts/sentaurus_gmsh_pipeline.py'),
                    '--sde',str(src),'--output-dir',str(out),'--runtime-dir',str(root/'runtime'),
                    '--maxtransdiff-max-iterations','0'],capture_output=True,text=True)
            self.assertEqual(0,run().returncode)
            self.assertTrue((out/'device.grd').exists())
            src.write_text(source_text+'\n(sdedr:define-refinement-function "S" "DopingConcentration" "MaxTransDiff" .000001)')
            args.maxtransdiff_max_iterations=0
            self.assertEqual(2,run().returncode)
            self.assertIn('iteration budget',json.loads((out/'fail_report.json').read_text())['message'])
            self.assertFalse((out/'run_manifest.json').exists())
            self.assertFalse(any((out/name).exists() for name in pipeline.GENERATED_FILES))

if __name__=='__main__':unittest.main()
