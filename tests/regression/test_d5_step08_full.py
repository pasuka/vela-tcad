import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from run_templates_ldmos_d5_step08_full import select_history,qf_limit_for


class HistoryTest(unittest.TestCase):
    def test_high_only_qf_preserves_low_voltage_baseline(self):
        self.assertEqual(qf_limit_for(1.433,.1,'high_only_qf'),.1)
        self.assertEqual(qf_limit_for(20.8,.8,'high_only_qf'),.2)
        self.assertEqual(qf_limit_for(21.333,.133,'high_only_qf'),.2)
        self.assertEqual(qf_limit_for(1.433,.1,'uniform_qf'),.2)

    def test_unpredicted_recovery_keeps_legacy_small_step_limit(self):
        self.assertEqual(qf_limit_for(12.9,.1,'protected_qf',False),.1)
        self.assertEqual(qf_limit_for(12.9,.1,'protected_qf',True),.2)
        self.assertEqual(qf_limit_for(20.8,.8,'protected_qf',True),.2)
        with self.assertRaises(ValueError):qf_limit_for(20,.8,'unknown')

    def test_old_frame_future_and_pre_recovery_states_are_excluded(self):
        points=[dict(bias_V=19.5,frame_V=28,exact=True),dict(bias_V=20.2,frame_V=0,exact=True),
                dict(bias_V=19,frame_V=0,exact=True)]
        self.assertIsNone(select_history(points,20,20.8,0,19.5))

    def test_clipped_nearest_point_can_use_safe_older_history(self):
        points=[dict(bias_V=19.99,frame_V=0,exact=True),dict(bias_V=19.2,frame_V=0,exact=False)]
        result=select_history(points,20,20.8,0,0)
        self.assertEqual(result['bias_V'],19.2)

    def test_recent_exact_history_is_preferred_with_ratio_guard(self):
        points=[dict(bias_V=20-4/3,frame_V=0,exact=True),dict(bias_V=19.2,frame_V=0,exact=False)]
        self.assertEqual(select_history(points,20,20.8,0,0),points[0])
        self.assertIsNone(select_history(points,20,20.8,0,19.8))

    def test_no_prediction_for_zero_or_backward_target(self):
        points=[dict(bias_V=19,frame_V=0,exact=True)]
        self.assertIsNone(select_history(points,20,20,0,0))
        self.assertIsNone(select_history(points,20,19,0,0))


if __name__=='__main__':unittest.main()
