import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from validate_wp14_zero_qs import BASELINE, source_checks, trajectory_checks
from validate_wp14_external_recovery import RESISTANCE, STOP


def state(outer, inner, current):
    row = {'time': outer / 100000}
    for name in ('drain', 'source', 'gate', 'substrate'):
        value = current if name == 'drain' else -current if name == 'source' else 0
        row[name] = {'total_current_A_per_um': value, 'e_current_A_per_um': value,
                     'h_current_A_per_um': 0, 'DisplacementCurrent': 0,
                     'outer_voltage_V': outer if name == 'drain' else 0,
                     'inner_voltage_V': inner if name == 'drain' else 0}
    return row


class ZeroQSAcquisitionTest(unittest.TestCase):
    def setUp(self):
        self.source = (BASELINE / 'external_resistor_sdevice.cmd').read_text()
        current = STOP * 1.001
        self.rows = [state(0, 0, 0), state(6.57 + RESISTANCE * current, 6.57, current)]
        self.log = 'Contact current exceeded limit:\n Good Bye !'

    def test_preserves_model_and_resistor_and_rejects_non_qs_seeds(self):
        self.assertTrue(all(source_checks(self.source, self.source).values()))
        for addition in ('Load(FilePrefix="seed")', 'Transient(InitialTime=0 FinalTime=1){}',
                         'Set("drain" mode current)', 'Continuation(){}'):
            source = self.source.replace('Solve {', 'Solve {\n' + addition)
            self.assertFalse(source_checks(source, self.source)['voltage_qs_only'])
        source = self.source.replace('Resistor=1e7', 'Resistor=1e6')
        self.assertFalse(source_checks(source, self.source)['same_electrodes_and_fixed_resistor'])
        source = self.source.replace('Avalanche(Eparallel)', 'Avalanche(GradQuasiFermi)')
        self.assertFalse(source_checks(source, self.source)['same_physics'])

    def test_comments_and_math_transient_do_not_execute_transient(self):
        source = self.source.replace('Solve {', 'Solve {\n * Transient and Load are historical comments\n')
        self.assertTrue(source_checks(source, self.source)['voltage_qs_only'])

    def test_alternative_source_approximation_is_a_separate_profile(self):
        source = self.source.replace('Math {', 'Math {\n AvalDensGradQF')
        checks = source_checks(source, self.source)
        self.assertTrue(checks['same_physics'])
        self.assertFalse(checks['original_avalanche_current_approximation'])

    def test_auxiliary_gate_path_cannot_qualify_the_original_fixed_gate_path(self):
        source = self.source.replace('Solve {', 'Solve {\n Quasistationary(Goal { Name="gate" Voltage=0.2 }) { Coupled { Poisson Electron Hole } }')
        self.assertFalse(source_checks(source, self.source)['voltage_qs_only'])
        rows = copy.deepcopy(self.rows)
        intermediate = copy.deepcopy(rows[-1])
        intermediate['gate']['outer_voltage_V'] = 0.2
        rows.insert(1, intermediate)
        self.assertFalse(trajectory_checks(rows, self.log)['checks']['original_fixed_bias_contacts_all_states'])

    def test_oversized_cancelled_currents_cannot_pass_a_stop_marker(self):
        rows = copy.deepcopy(self.rows)
        rows[-1]['drain'].update(total_current_A_per_um=8.689e129,
                                 e_current_A_per_um=8.501e146,
                                 h_current_A_per_um=-8.501e146)
        self.assertEqual('fail', trajectory_checks(rows, self.log)['status'])

    def test_zero_qs_stop_with_load_line_and_kcl_passes(self):
        self.assertEqual('pass', trajectory_checks(self.rows, self.log)['status'])

    def test_good_bye_cannot_hide_min_step_failure_or_missing_stop(self):
        result = trajectory_checks(self.rows, self.log + '\nStep-size less than MinStep')
        self.assertFalse(result['checks']['native_success'])
        self.assertEqual('fail', result['status'])
        self.assertFalse(trajectory_checks(self.rows, 'Good Bye !')['checks']['native_current_stop_reported'])

    def test_restart_backtracking_and_capacitive_current_do_not_prove_acquisition(self):
        self.assertFalse(trajectory_checks(self.rows[1:], self.log)['checks']['starts_at_zero_voltage'])
        self.assertFalse(trajectory_checks(self.rows + self.rows[:1], self.log)['checks']['monotone_outer_voltage'])
        rows = copy.deepcopy(self.rows)
        rows[-1]['drain']['DisplacementCurrent'] = STOP
        rows[-1]['drain']['e_current_A_per_um'] -= STOP
        self.assertFalse(trajectory_checks(rows, self.log)['checks']['negligible_displacement_all_states'])

    def test_a_bad_intermediate_state_is_not_hidden_by_a_good_endpoint(self):
        rows = copy.deepcopy(self.rows)
        intermediate = copy.deepcopy(rows[-1])
        intermediate['drain']['outer_voltage_V'] -= 0.1
        rows.insert(1, intermediate)
        self.assertFalse(trajectory_checks(rows, self.log)['checks']['fixed_resistor_all_states'])
        rows = copy.deepcopy(self.rows)
        rows[-1]['source']['total_current_A_per_um'] *= 0.9
        self.assertFalse(trajectory_checks(rows, self.log)['checks']['terminal_kcl_all_states'])


if __name__ == '__main__':
    unittest.main()
