import os
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import tempfile

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from run_templates_ldmos_thread_pilot import configurations, environment, audit_threads
from analyze_templates_ldmos_thread_pilot import analyze
from run_templates_ldmos_linked_d5 import write


class ThreadPilotTests(unittest.TestCase):
    def test_analysis_keeps_rejected_cases_without_speed_claim(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            write(root/'summary.json',dict(status='completed_with_rejections',
                configurations=[dict(key='mumps_t2',backend='mumps',level=2)],
                cases=[dict(key='mumps_t2',status='not_run',reason='fixed_matrix_screen_rejected')],
                screen=[dict(key='mumps_t2',round=0,status='rejected')],runner_sha256='frozen'))
            result=analyze(root)
            self.assertEqual(result['cases'][0]['status'],'not_run')
            self.assertNotIn('speed_ratio',result['cases'][0])

    def test_analysis_rejects_incomplete_batch_and_incomplete_screening(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            report=dict(status='running',configurations=[],cases=[],screen=[])
            write(root/'summary.json',report)
            with self.assertRaisesRegex(ValueError,'not complete'):analyze(root)
            report.update(status='pass',configurations=[dict(key='umfpack_t1',backend='umfpack')],
                          cases=[dict(key='umfpack_t1',status='pass')])
            write(root/'summary.json',report)
            with self.assertRaisesRegex(ValueError,'three qualified'):analyze(root)
            report['cases']=[]
            write(root/'summary.json',report)
            with self.assertRaisesRegex(ValueError,'Missing'):analyze(root)

    def test_serial_backends_are_not_mislabeled_as_parallel_solvers(self):
        rows=configurations()
        self.assertEqual(len(rows),13)
        self.assertEqual(len({r['key'] for r in rows}),13)
        self.assertEqual([r['level'] for r in rows if r['backend']=='sparselu'],[1])
        for r in rows:
            if r['backend']=='umfpack':
                self.assertEqual(r['solver_threads'],1)
                self.assertEqual(r['blas_threads'],r['level'])
                self.assertEqual(r['omp_threads'],r['level'])
            else:self.assertEqual(r['blas_threads'],1)

    def test_inherited_settings_cannot_override_experiment_or_enable_capture(self):
        cfg=next(c for c in configurations() if c['key']=='umfpack_t4')
        with patch.dict(os.environ,{'OMP_NUM_THREADS':'1','VELA_BLAS_THREADS':'1','VELA_LINEAR_CAPTURE_DIR':'unexpected'}):
            env=environment(cfg)
        self.assertEqual(env['OMP_NUM_THREADS'],'4')
        self.assertEqual(env['VELA_BLAS_THREADS'],'4')
        self.assertEqual(env['VELA_LINEAR_THREADS'],'1')
        self.assertEqual(env['OMP_MAX_ACTIVE_LEVELS'],'1')
        self.assertNotIn('VELA_LINEAR_CAPTURE_DIR',env)

    def test_thread_audit_requires_runtime_observations_and_rejects_drift(self):
        cfg=next(c for c in configurations() if c['key']=='mumps_t2')
        obs={k:dict(min=v,max=v,count=2) for k,v in {
            'linear.openblas_threads':1,'linear.mumps_omp_threads':2,'linear.mumps_omp_max_active_levels':1}.items()}
        profile=dict(counters={'linear.solve_calls':2},observations=obs)
        self.assertTrue(audit_threads([profile],cfg)['verified'])
        obs['linear.mumps_omp_threads']['max']=4
        with self.assertRaises(ValueError):audit_threads([profile],cfg)
        with self.assertRaises(ValueError):audit_threads([],cfg)
        with self.assertRaises(ValueError):audit_threads([dict(counters={'linear.solve_calls':1})],cfg)


if __name__=='__main__':unittest.main()
