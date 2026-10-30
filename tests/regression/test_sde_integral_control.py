"""Physical integral convergence and resource bounds, independent of mesh topology."""
from pathlib import Path
import json,math,subprocess,sys,tempfile,unittest
from unittest.mock import patch
REPO=Path(__file__).resolve().parents[2];sys.path.insert(0,str(REPO/'scripts'))
from sentaurus_device_ir import parse_sde_device_ir
from sentaurus_ir_v2 import migrate_device_ir_v1
from sentaurus_integral import analyse
from tests.regression.gmsh_test_runtime import require_gmsh_runtime


SOURCE='''(sdegeo:create-rectangle (position 0 0 0) (position .1 .2 0) "Silicon" "R")
(sdedr:define-refeval-window "B" "Line" (position 0 0 0) (position .1 0 0))
(sdedr:define-gaussian-profile "G" "ArsenicActiveConcentration" "PeakPos" 0 "PeakVal" 1e18 "ValueAtDepth" 1e16 "Depth" .1 "Gauss" "Factor" 0)
(sdedr:define-analytical-profile-placement "P" "G" "B" "Both" "NoReplace" "Eval")
(sde:build-mesh "test")'''


class IntegralTests(unittest.TestCase):
    def test_triangle_rule_integrates_a_quadratic_density(self):
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)/'input.cmd';source.write_text('(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) \"Silicon\" \"R\")')
            ir=migrate_device_ir_v1(parse_sde_device_ir(source));owner=ir['geometry']['bodies'][0]
            with patch('sentaurus_doping._evaluate_m2_doping_for_owner',
                       side_effect=lambda ir,o,r,w,x,y:(x*x+y*y,0.,{})):
                marks,report=analyse(ir,{0:(0.,0.,0.),1:(1.,0.,0.),2:(0.,1.,0.)},
                                     {1:[[0,1,2]]},{1:owner},.01,{})
            self.assertAlmostEqual(1/6,report['quadrature_integrals_cm3_um2'][0],places=13)
            self.assertAlmostEqual(1/3,report['linear_integrals_cm3_um2'][0],places=13)
            self.assertTrue(marks)

    def test_gaussian_dose_converges_to_analytic_half_space_integral(self):
        require_gmsh_runtime()
        with tempfile.TemporaryDirectory(prefix='p6-integral-') as temp:
            root=Path(temp);source=root/'input.cmd';source.write_text(SOURCE)
            errors=[]
            char=.1/math.sqrt(math.log(100));exact=.1*1e18*char*math.sqrt(math.pi)/2*math.erf(.2/char)
            # This deliberately coarse mesh gets an explicit 20-round budget;
            # the product default remains 12 and exhaustion is tested below.
            for target in [0.,.00025]:
                out=root/str(target)
                cmd=[sys.executable,str(REPO/'scripts/sentaurus_gmsh_pipeline.py'),'--sde',str(source),
                     '--capability-level','M3','--output-dir',str(out),'--runtime-dir',str(root/'runtime'),
                     '--integral-relative-target',str(target),'--maxtransdiff-max-iterations','20']
                result=subprocess.run(cmd,cwd=REPO,capture_output=True,text=True)
                detail=(out/'fail_report.json').read_text() if (out/'fail_report.json').exists() else ''
                self.assertEqual(0,result.returncode,result.stderr+detail)
                bundle=json.loads((out/'mesh_field_bundle.json').read_text());points=bundle['nodes'];dose=0.
                for t in bundle['triangles']:
                    ids=t['node_ids'];a,b,c=[points[i] for i in ids]
                    area=abs((b['x']-a['x'])*(c['y']-a['y'])-(b['y']-a['y'])*(c['x']-a['x']))/2
                    dose+=area*sum(bundle['fields']['DonorConcentration']['values'][i] for i in ids)/3
                errors.append(abs(dose-exact)/exact)
                if target:
                    report=bundle['provenance']['maxtransdiff']['integral_control']
                    self.assertLessEqual(max(report['relative_errors']),target)
                    self.assertLessEqual(max(report['relative_error_bounds']),target)
            self.assertLess(errors[1],.001)
            self.assertLess(errors[1],errors[0]/10)

    def test_integral_refinement_cannot_bypass_iteration_budget(self):
        require_gmsh_runtime()
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);source=root/'input.cmd';source.write_text(SOURCE);out=root/'out'
            result=subprocess.run([sys.executable,str(REPO/'scripts/sentaurus_gmsh_pipeline.py'),
                '--sde',str(source),'--capability-level','M3','--output-dir',str(out),
                '--runtime-dir',str(root/'runtime'),'--integral-relative-target','.00025',
                '--maxtransdiff-max-iterations','0'],cwd=REPO,capture_output=True,text=True)
            self.assertEqual(2,result.returncode)
            self.assertFalse((out/'run_manifest.json').exists())
            self.assertEqual('iteration_budget',json.loads((out/'fail_report.json').read_text())['details']['reason'])


if __name__=='__main__':unittest.main()
