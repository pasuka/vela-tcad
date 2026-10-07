import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from validate_wp14_gate_qs import gate_source_checks, gate_trajectory_checks
from validate_wp14_zero_qs import BASELINE, source_checks
from tests.regression.test_wp14_zero_qs import state


def held(gate_voltage, inner):
    row = state(13000, inner, (13000 - inner) / 1e7)
    row['gate']['outer_voltage_V'] = gate_voltage
    row['gate']['inner_voltage_V'] = gate_voltage
    return row


class AuxiliaryGateQSTest(unittest.TestCase):
    def setUp(self):
        baseline = (BASELINE / 'external_resistor_sdevice.cmd').read_text(encoding='utf-8')
        self.baseline = baseline
        self.source = baseline[:baseline.index('Solve {')] + '''Solve {
          Coupled { Poisson Electron Hole }
          Quasistationary(Goal{Name="drain" Voltage=13000}){Coupled{Poisson Electron Hole}}
          NewCurrentPrefix="gate_seed_"
          Quasistationary(Goal{Name="gate" Voltage=0.2}){Coupled{Poisson Electron Hole}}
          Quasistationary(Goal{Name="gate" Voltage=0}){Coupled{Poisson Electron Hole}}
          NewCurrentPrefix="final_"
          Quasistationary(Goal{Name="drain" Voltage=100000}){Coupled{Poisson Electron Hole}}
        }'''
        start, restored = held(0, 7.13), held(0, 6.64)
        self.phases = [[state(0, 0, 0), start],
                       [copy.deepcopy(start), held(0.1, 6.9), held(0.2, 6.8), held(0.1, 6.7), restored],
                       [copy.deepcopy(restored), state(14506.57, 6.57, 0.00145)]]
        self.log = 'Contact current exceeded limit:\n Good Bye !'

    def test_authorized_source_is_distinct_from_original_fixed_gate_scope(self):
        self.assertTrue(all(gate_source_checks(self.source, self.baseline).values()))
        self.assertFalse(source_checks(self.source, self.baseline)['voltage_qs_only'])

    def test_seed_and_changed_model_or_goal_cannot_qualify(self):
        for token in ('Load(FilePrefix="seed")', 'Transient(){}', 'Set("drain" mode current)'):
            checks = gate_source_checks(self.source.replace('Solve {', 'Solve {' + token), self.baseline)
            self.assertFalse(checks['qs_only_no_loaded_seed'])
        self.assertFalse(gate_source_checks(self.source.replace('Voltage=0.2', 'Voltage=0.3'), self.baseline)['authorized_goal_sequence'])
        self.assertFalse(gate_source_checks(self.source.replace('Resistor=1e7', 'Resistor=5e7'), self.baseline)['same_electrodes_and_fixed_resistor'])

    def test_complete_physically_closed_auxiliary_path_passes(self):
        result = gate_trajectory_checks(self.phases, self.log)
        self.assertEqual('pass', result['status'])
        self.assertEqual([2, 5, 2], result['phase_state_counts'])

    def test_gate_must_reach_peak_and_restore_zero(self):
        for index, value in ((2, 0.15), (-1, 0.01)):
            phases = copy.deepcopy(self.phases)
            phases[1][index]['gate']['outer_voltage_V'] = value
            self.assertFalse(gate_trajectory_checks(phases, self.log)['checks']['gate_sweep_0_to_0_2_to_0'])

    def test_prefix_discontinuity_and_drain_drift_are_rejected(self):
        phases = copy.deepcopy(self.phases)
        phases[2][0]['drain']['inner_voltage_V'] -= 0.1
        self.assertFalse(gate_trajectory_checks(phases, self.log)['checks']['phase_states_continuous'])
        phases = copy.deepcopy(self.phases)
        phases[1][2]['drain']['outer_voltage_V'] += 1
        self.assertFalse(gate_trajectory_checks(phases, self.log)['checks']['gate_phase_drain_held_13000V'])

    def test_bad_middle_state_or_native_failure_is_not_hidden(self):
        phases = copy.deepcopy(self.phases)
        phases[1][2]['source']['total_current_A_per_um'] *= 0.9
        self.assertFalse(gate_trajectory_checks(phases, self.log)['checks']['terminal_kcl_all_states'])
        self.assertEqual('fail', gate_trajectory_checks(self.phases, self.log + '\nStep-size less than MinStep')['status'])
        self.assertEqual('fail', gate_trajectory_checks(self.phases, 'Good Bye !')['status'])

    def test_final_scan_cannot_stop_while_gate_is_biased(self):
        phases = copy.deepcopy(self.phases)
        phases[2][-1]['gate']['outer_voltage_V'] = 0.1
        self.assertFalse(gate_trajectory_checks(phases, self.log)['checks']['gate_zero_during_acquisition_and_final'])

    def test_continuing_after_failed_gate_ramp_does_not_hide_the_missing_goal(self):
        phases = copy.deepcopy(self.phases)
        for row in phases[1]:
            row['gate']['outer_voltage_V'] *= 0.483477
        # Native may continue later Solve statements and report a final stop.
        # An incomplete auxiliary sweep must fail even with a closed endpoint.
        result = gate_trajectory_checks(phases, self.log)
        self.assertTrue(result['checks']['original_stop_reached'])
        self.assertFalse(result['checks']['gate_sweep_0_to_0_2_to_0'])
        self.assertEqual('fail', result['status'])


if __name__ == '__main__':
    unittest.main()
