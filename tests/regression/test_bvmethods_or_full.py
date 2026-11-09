import copy
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts/diagnostics'))
from bvmethods_or_full import or_gate, gate_summary, ORSupervisor


def row():
    g = dict(coupled_device_residual_inf=1.01157e-8, actual_load_line_residual_V=1e-8,
             density_update_rms=.254, density_update_valid=True, density_update_mode='enforce',
             residual_stopping_scope='device_and_circuit', load_line_enforced=True,
             load_line_convergence_mode='enforce', physical_residuals_finite=True,
             terminal_kcl_satisfied=True, actual_load_line_satisfied=True,
             carrier_row_enforced=False, carrier_row_satisfied=False,
             global_continuity_enforced=False, global_continuity_satisfied=False)
    return dict(validation_diagnostics=json.dumps(g), outer_voltage_V='.56',
                inner_voltage_V='.53', current_total_A_per_um='2.5e-9', iterations='3')


def change(item, **kw):
    result=copy.deepcopy(item);g=json.loads(result['validation_diagnostics']);g.update(kw)
    result['validation_diagnostics']=json.dumps(g);return result


class ORFullTests(unittest.TestCase):
    def test_update_branch_and_independent_gates(self):
        r=row();self.assertTrue(or_gate(r))
        for kw in [dict(actual_load_line_residual_V=.01),dict(terminal_kcl_satisfied=False),
                   dict(carrier_row_enforced=True),dict(global_continuity_enforced=True),
                   dict(load_line_enforced=False),dict(coupled_device_residual_inf=float('nan')),
                   dict(density_update_valid=False)]:
            self.assertFalse(or_gate(change(r,**kw)),kw)

    def test_residual_branch_does_not_require_small_update(self):
        r=change(row(),coupled_device_residual_inf=4e-9,density_update_rms=1083.)
        self.assertTrue(or_gate(r));self.assertTrue(or_gate(change(r,density_update_valid=False,density_update_rms=None)))
        self.assertFalse(or_gate(change(r,actual_load_line_residual_V=.1)))

    def test_gate_summary_preserves_report_failures_and_branches(self):
        a=row();b=change(a,coupled_device_residual_inf=4e-9,density_update_rms=1083.)
        result=gate_summary([a,b]);self.assertEqual(result['update_only'],1)
        self.assertEqual(result['residual_only'],1);self.assertEqual(result['failed_points'],0)
        self.assertEqual(result['carrier_report_unsatisfied_points'],2)
        self.assertEqual(result['max_valid_update_rms'],1083.)

    def test_configuration_only_selects_OR_in_addition_to_existing_group_policy(self):
        from bvmethods_full_ablation import Supervisor
        template={'solver':{},'runtime_log':{},'sweep':{'external_circuit':{'coupled_update_convergence':{}}}}
        controller=object.__new__(ORSupervisor)
        for group in ['G1','G2']:
            expected=Supervisor.config(controller,template,Path('case'),group)
            expected['sweep']['external_circuit']['coupled_update_convergence']['mode']='enforce'
            self.assertEqual(controller.config(template,Path('case'),group),expected)
        self.assertEqual(template['solver'],{})


if __name__=='__main__':unittest.main()
