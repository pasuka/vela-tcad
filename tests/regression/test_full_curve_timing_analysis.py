"""Timing ranges and voltage partitioning preserve complete solver work."""
import sys
from pathlib import Path
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import analyze_templates_ldmos_full_curve_timing as analysis


class FullCurveAnalysisTests(unittest.TestCase):
    def test_statistics_keep_all_runs_instead_of_fastest(self):
        result=analysis.repeat_statistics([12.,8.,10.])
        self.assertEqual(result['median'],10.)
        self.assertEqual(result['relative_range_percent'],40.)
        self.assertEqual(result['values'],[12.,8.,10.])
        self.assertEqual(analysis.repeat_statistics([5.])['count'],1)
        for bad in ([],[0.],[-1.]):
            with self.assertRaises(ValueError):analysis.repeat_statistics(bad)

    def test_voltage_segments_retain_repeated_and_failed_attempt_work(self):
        runs=[dict(target_V=v,Newton_updates=n,case=str(i)) for i,(v,n) in
              enumerate(((0.,0),(9.33333333333333,3),(9.34,5),(9.34,2)))]
        profile=dict(counters={'linear.solve_calls':5,'linear.analyze_calls':1},
            stages=[dict(name='newton.linear_l2_row_column',calls=5,total_ns=1000000000)])
        with patch.object(analysis,'read',side_effect=[{'runs':runs},profile,profile,profile,profile]):
            result=analysis.profile_breakdown(Path('unused'))
        self.assertEqual(result['prefix_0_to_9p333333']['updates'],3)
        self.assertEqual(result['remaining_high_voltage']['updates'],7)
        self.assertEqual(sum(x['requests'] for x in result.values()),4)
        self.assertEqual(result['remaining_high_voltage']['main_only_analyses'],2)


if __name__=='__main__':unittest.main()
