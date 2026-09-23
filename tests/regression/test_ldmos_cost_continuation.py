import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from ldmos_cost_continuation import work_units,next_proposal,Timings
from analyze_templates_ldmos_d5_cost_controls import analyze

class CostContinuationTest(unittest.TestCase):
    def test_more_rejected_work_cannot_increase_step(self):
        a=work_units({'linear.factorize_calls':4,'dd.jacobian_calls':4,'dd.residual_calls':4})
        b=work_units({'linear.factorize_calls':4,'dd.jacobian_calls':4,'dd.residual_calls':20})
        self.assertGreater(b,a)
        self.assertLess(next_proposal(.2,.2,b,1,False)[0],next_proposal(.2,.2,a,1,False)[0])
    def test_expensive_recovery_and_damping_shrink(self):
        for work,alpha,recovery in ((20,1,False),(2,.001,False),(2,1,True)):
            self.assertLess(next_proposal(.2,.2,work,alpha,recovery)[0],.2)
    def test_exact_output_fragment_does_not_grow_or_collapse_trusted_step(self):
        self.assertEqual(next_proposal(.4,.01,1,1,False)[0],.4)
        self.assertLess(next_proposal(.4,.01,20,1,False)[0],.01)
    def test_limits_and_invalid_metrics(self):
        self.assertLessEqual(next_proposal(.8,.8,1,1,False)[0],.8)
        self.assertGreaterEqual(next_proposal(.0025,.0025,30,1,False)[0],.0025)
        with self.assertRaises(ValueError):next_proposal(.2,.2,float('nan'),1,False)
        with self.assertRaises(ValueError):work_units({'dd.residual_calls':-1})
    def test_nested_timer_has_no_double_counting_even_on_exception(self):
        ticks=iter([0,1,3,5]);t=Timings(lambda:next(ticks))
        with self.assertRaises(ValueError):
            with t.stage('all'):
                with t.stage('child'):raise ValueError('expected')
        self.assertEqual(t.data['all']['inclusive'],5)
        self.assertEqual(t.data['all']['exclusive']+t.data['child']['exclusive'],5)
        self.assertEqual(t.stack,[])

class EvidenceTest(unittest.TestCase):
    def batch(self):
        return dict(status='completed',planned_curves=1,points=8,cases=[dict(
            status='completed',mode='baseline',gate=4,repeat=0,wall_seconds=5.,
            exact_points=8,trajectory_exact=True,max_potential_difference_V=0.,
            audit=dict(integrity_pass=True,Newton_updates=5,services=2,advances=1,rollbacks=1),
            performance=dict(blas_threads=[1],counters={'linear.factorize_calls':7,
                'dd.jacobian_calls':7,'dd.residual_calls':10},child_wall_seconds=3.,child_cpu_seconds=3.),
            parent_timing=dict(controller=dict(inclusive=5.,exclusive=2.),
                worker_roundtrip=dict(inclusive=3.,exclusive=3.)))])
    def test_attempt_cost_and_parent_exclusive_survive_analysis(self):
        row=analyze(self.batch())['rows'][0]
        self.assertEqual(row['rollbacks'],1)
        self.assertGreater(row['factorization_calls'],row['updates'])
        self.assertEqual(row['parent_exclusive_seconds'],2.)
    def test_incomplete_and_failed_states_are_not_qualified(self):
        b=self.batch();b['status']='running'
        with self.assertRaises(ValueError):analyze(b)
        b=self.batch();b['cases'][0]['max_potential_difference_V']=2e-8
        with self.assertRaises(ValueError):analyze(b)
    def test_overlapping_and_nonfinite_timing_are_rejected(self):
        for value in (3.,float('nan')):
            b=self.batch();b['cases'][0]['parent_timing']['controller']['exclusive']=value
            with self.assertRaises(ValueError):analyze(b)
    def test_cache_cannot_replace_final_cold_audit_or_change_trajectory(self):
        b=self.batch();b['cases'][0]['mode']='files'
        with self.assertRaisesRegex(ValueError,'uncached'):analyze(b)
        b['cases'][0]['cold_final_audit']=True
        self.assertEqual(analyze(b)['status'],'pass')
        b['cases'][0]['trajectory_exact']=False
        with self.assertRaisesRegex(ValueError,'baseline'):analyze(b)

if __name__=='__main__':unittest.main()
