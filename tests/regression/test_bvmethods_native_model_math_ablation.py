"""Contracts for isolated native sensitivity experiments, independent of sdevice."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from bvmethods_native_model_math_ablation import SOURCE, VARIANTS, make_variant, classify_vela_control
from validate_wp14_zero_qs import block, trajectory_checks


class NativeAblationContract(unittest.TestCase):
    def test_each_control_preserves_circuit_scan_and_other_section(self):
        source = (SOURCE / 'external_resistor_sdevice.cmd').read_text()
        for name, (section, _, _) in VARIANTS.items():
            with self.subTest(name=name):
                candidate = make_variant(source, name)
                for key in ('Electrode', 'Solve', 'Plot', 'File'):
                    self.assertEqual(block(candidate, key), block(source, key))
                if section == 'Physics':
                    self.assertEqual(block(candidate, 'Math'), block(source, 'Math'))
                elif section == 'Math':
                    self.assertEqual(block(candidate, 'Physics'), block(source, 'Physics'))
                else:
                    self.assertEqual(candidate, source)

    def test_changed_source_is_rejected_instead_of_silent_noop(self):
        source = (SOURCE / 'external_resistor_sdevice.cmd').read_text()
        with self.assertRaises(ValueError):
            make_variant(source.replace('Band2Band(E2)', 'Band2Band(Hurkx)'), 'no_btbt')

    def test_math_derivative_switch_is_explicit_and_keeps_avalanche_model(self):
        source = (SOURCE / 'external_resistor_sdevice.cmd').read_text()
        candidate = make_variant(source, 'no_aval_derivatives')
        self.assertIn('-AvalDerivatives', block(candidate, 'Math'))
        self.assertIn('Avalanche(Eparallel)', block(candidate, 'Physics'))

    def test_vela_control_identity_comes_from_actual_selector(self):
        config = {'solver': {'band_to_band': {'model': 'e2'},
                            'mobility': {'model': 'masetti_field_lombardi',
                                         'jacobian_field_derivatives': True}}}
        self.assertEqual(classify_vela_control(config), 'baseline')
        config['solver']['mobility']['jacobian_field_derivatives'] = False
        self.assertEqual(classify_vela_control(config), 'frozen_hfs_jacobian')
        config['solver']['band_to_band']['model'] = 'none'
        with self.assertRaises(ValueError):
            classify_vela_control(config)

    def test_good_bye_after_minstep_cannot_pass_trajectory(self):
        def state(current):
            row = {'time': 0}
            for contact in ('drain', 'source', 'gate', 'substrate'):
                sign = 1 if contact == 'drain' else -1 if contact == 'substrate' else 0
                value = sign * current
                row[contact] = {'inner_voltage_V': 0, 'outer_voltage_V': 1e7 * current if contact == 'drain' else 0,
                                'total_current_A_per_um': value, 'e_current_A_per_um': value,
                                'h_current_A_per_um': 0, 'DisplacementCurrent': 0}
            return row
        result = trajectory_checks([state(0), state(0.00145)],
                                    'Contact current exceeded limit:\nStep-size less than MinStep\nGood Bye')
        self.assertEqual(result['status'], 'fail')
        self.assertFalse(result['checks']['native_success'])


if __name__ == '__main__':
    unittest.main()
