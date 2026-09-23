import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from analyze_templates_ldmos_d5_newton_cost import trace_cost, native_tables, segment, recovery_counts


class CostTest(unittest.TestCase):
    def test_success_does_not_erase_internal_recovery_work(self):
        curve=[dict(converged='1',carrier_row_recovery_attempted='1',
                    carrier_row_recovery_cycles='2',carrier_row_recovery_density_passes='4'),
               dict(converged='1',carrier_row_recovery_attempted='0')]
        self.assertEqual(recovery_counts(curve),dict(internal_recovery_points=1,
                         internal_recovery_cycles=2,internal_density_passes=4))

    def test_row_tail_is_measured_before_update(self):
        gates=dict(psi_residual_ceiling=5e-8,electron_residual_ceiling=1e-11,hole_residual_ceiling=3e-10)
        initial=dict(run_id='0',segment_id='0',attempt_id='1',event='initial',residual_norm='1e-12',
            block_psi='1e-10',block_phin='1e-12',block_phip='1e-12',carrier_row_max_ratio='1e-3')
        update=dict(initial,event='accepted_iteration',residual_norm='1e-13',damping='0.5',
                    line_search_attempts='3',carrier_row_max_ratio='1e-10')
        c,_=trace_cost([initial,update],gates,1e-8)
        self.assertEqual(c['row_only_tail_updates'],1)
        self.assertEqual(c['damped_updates'],1)
        self.assertEqual(c['trace_extra_trials'],2)
        self.assertEqual(trace_cost([initial],gates,1e-8)[0]['updates'],0)
        with self.assertRaises(ValueError): trace_cost([update],gates,1e-8)

    def test_native_iteration_zero_is_not_an_update(self):
        text='Iteration   |Rhs| factor step error inner iterative time\n----\n'
        text+='0  1.0  0.01\n1 1e-4 1.0 0.1 0.1 0 1 0.2\n2 1e-6 0.5 0.1 0.01 0 1 0.3\n'
        text+='Finished, because...\nError smaller than 1\n'
        self.assertEqual(native_tables(text)[0],dict(updates=2,damped=1,stop='Error smaller than 1',success=True))

    def test_segment_boundaries_tolerate_native_print_precision(self):
        self.assertEqual(segment(9.33334),'low_0_9p333')
        self.assertEqual(segment(26.66668),'mid_9p333_26p667')
        self.assertEqual(segment(40),'high_26p667_40')


if __name__=='__main__': unittest.main()
