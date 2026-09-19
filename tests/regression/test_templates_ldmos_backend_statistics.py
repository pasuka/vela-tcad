"""Aggregation checks for heterogeneous point-service sample counts."""
import sys
import unittest
import tempfile
from types import SimpleNamespace
from unittest.mock import patch
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from analyze_templates_ldmos_isothermal_backends import factor_statistics
from windows_system_load import cpu_interval
from analyze_templates_ldmos_multi_backend import matrix_summary, curve_summary
from run_templates_ldmos_isothermal_backends import archive_interrupted_artifacts, validate_resume_settings, backend_threads


class StatisticsTests(unittest.TestCase):
    def test_threaded_candidates_keep_serial_control_fixed(self):
        for requested in [1,2,4]:
            self.assertEqual(backend_threads('strumpack',requested),requested)
            self.assertEqual(backend_threads('superlu_mt',requested),requested)
            self.assertEqual(backend_threads('umfpack',requested),1)
            self.assertEqual(backend_threads('sparselu',requested),1)

    def test_resume_does_not_change_experiment_settings(self):
        record=dict(backends=['sparselu','umfpack'],threads=1,points=8,rounds=3,
                    profiles=['D5'],factor_statistics='off',idle_cpu_percent=None)
        validate_resume_settings(record,SimpleNamespace(**record))
        for key,value in dict(threads=2,points=31,profiles=['D4'],factor_statistics='on').items():
            with self.subTest(key=key),self.assertRaises(ValueError):
                validate_resume_settings(record,SimpleNamespace(**(record|{key:value})))

    def test_interrupted_archive_preserves_bytes_and_rejects_collisions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);case=root/'case';case.mkdir();(case/'state.csv').write_bytes(b'original state')
            (root/'case.log').write_bytes(b'partial log')
            (root/'case.log.interrupted_t').write_bytes(b'prior evidence')
            with self.assertRaises(ValueError):archive_interrupted_artifacts(root,'case','t')
            self.assertTrue(case.exists())
            moves=archive_interrupted_artifacts(root,'case','next')
            self.assertEqual(len(moves),2)
            self.assertEqual((root/'case.interrupted_next/state.csv').read_bytes(),b'original state')
            self.assertEqual((root/'case.log.interrupted_next').read_bytes(),b'partial log')
            self.assertFalse(case.exists())

    def test_load_is_cpu_capacity_not_single_thread_percent(self):
        a=dict(total=0,idle=0,monotonic=10.)
        b=dict(total=80000000,idle=60000000,monotonic=12.)
        self.assertEqual(cpu_interval(a,b),dict(interval_seconds=2.,busy_percent=25.,busy_cpu_seconds=2.))
        for invalid in [dict(total=0,idle=0,monotonic=12.),dict(total=1,idle=2,monotonic=12.)]:
            with self.assertRaises(ValueError):cpu_interval(a,invalid)
    def test_weighted_factors_and_process_peak(self):
        def profile(count, mean, peak):
            return dict(observations={'linear.numeric_factor_fill_ratio':
                dict(count=count, average=mean, min=mean, max=mean)},
                resources={'process_peak_working_set_bytes':peak})
        value=factor_statistics([profile(1,2.,100),profile(9,4.,200)])
        obs=value['observations']['linear.numeric_factor_fill_ratio']
        self.assertEqual(obs['count'],10)
        self.assertEqual(obs['average'],3.8)
        self.assertEqual(obs['total'],38.)
        self.assertEqual(obs['min'],2.)
        self.assertEqual(obs['max'],4.)
        self.assertEqual(value['max_solver_process_peak_working_set_bytes'],200)

    def test_missing_is_not_zero(self):
        value=factor_statistics([dict(observations={})])
        self.assertIsNone(value['max_solver_process_peak_working_set_bytes'])
        self.assertEqual(value['observations'],{})

    def test_matrix_rejections_and_missing_memory_remain_visible(self):
        summary=dict(status='completed_with_rejections',threads=2,cases=[
            dict(backend='superlu_mt',status='pass',round=0,wall_seconds=3.),
            dict(backend='mumps_metis',status='rejected',round=0,wall_seconds=1.)])
        systems=dict(systems=[dict(profiling=dict(stages=[dict(name='linear.total',total_ns=2e9)]),
            cumulative_analyses=1,quality=dict(raw_available=False,scaled=dict(normwise_backward_error=1e-17)))])
        with patch('analyze_templates_ldmos_multi_backend.read',side_effect=[summary,systems]):
            value=matrix_summary(Path('unused'))
        row=value['groups']['superlu_mt']['samples'][0]
        self.assertIsNone(row['peak_process_bytes'])
        self.assertFalse(row['raw_available'])
        self.assertEqual(len(value['groups']['mumps_metis']['failures']),1)

    def test_curve_cpu_includes_children_and_preserves_backend_name(self):
        summary=dict(status='pass',comparisons=[],cases=[dict(name='d5_vg4_r0_superlu_mt',
            backend='superlu_mt',status='pass',wall_seconds=20.,cpu_seconds=dict(total=2.),
            audit=dict(child_cpu_seconds=13.,Newton_updates=7))])
        ledger=dict(runs=[],total_Newton_updates=7,rollbacks=[])
        with patch('analyze_templates_ldmos_multi_backend.read',side_effect=[summary,ledger]), \
             patch('analyze_templates_ldmos_multi_backend.repeat_statistics',return_value=[]):
            value=curve_summary(Path('unused'))
        row=value['cases'][0]
        self.assertEqual(row['backend'],'superlu_mt')
        self.assertEqual(row['process_tree_cpu_seconds'],15.)
        self.assertEqual(row['parent_cpu_seconds'],2.)


if __name__=='__main__':
    unittest.main()
