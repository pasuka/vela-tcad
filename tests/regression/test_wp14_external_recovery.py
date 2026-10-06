import sys
import unittest
import re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from validate_wp14_external_recovery import check_state, STOP, RESISTANCE


class FixedResistorStateTest(unittest.TestCase):
    def test_qualified_sources_preserve_models_parameters_and_valid_load_prefix(self):
        root = Path(__file__).resolve().parents[2]
        original = root / 'reference_tcad/bvmethods_sentaurus2018/source'
        qualified = root / 'reference_tcad/bvmethods_sentaurus2018/wp14_sources/v4'
        self.assertEqual((original/'models.par').read_bytes(), (qualified/'models.par').read_bytes())
        physics = lambda text: re.sub(r'\s+', '', re.search(r'Physics\s*\{([^}]*)\}', text).group(1))
        baseline = physics((original/'external_resistor_sdevice.cmd').read_text())
        for name in ('external_transient_saved','external_dc_saved','external_recovery_saved'):
            text = (qualified/(name+'_sdevice.cmd')).read_text()
            self.assertEqual(baseline, physics(text))
            self.assertIn('Resistor=1e7', text)
            if name != 'external_dc_saved': self.assertIn('AbsVal=1.443e-3', text)
            if name != 'external_transient_saved':
                self.assertIn('Load(FilePrefix="../external_transient_saved_v2/external_transient_seed")', text)
            if name == 'external_recovery_saved': self.assertIn('Voltage=100000.0', text)

    def state(self):
        current = STOP * 1.01
        row = {name: {'total_current_A_per_um': 0.0} for name in ('drain','source','gate','substrate')}
        row['drain'] = {'total_current_A_per_um': current, 'e_current_A_per_um': current,
                        'h_current_A_per_um': 0.0, 'DisplacementCurrent': 0.0,
                        'inner_voltage_V': 6.57, 'outer_voltage_V': 6.57 + RESISTANCE * current}
        row['source']['total_current_A_per_um'] = -current
        return row

    def test_fixed_resistor_conductive_state_passes_and_broken_load_line_fails(self):
        row = self.state()
        self.assertEqual('pass', check_state(row)['status'])
        row['drain']['inner_voltage_V'] = row['drain']['outer_voltage_V']
        self.assertFalse(check_state(row)['checks']['fixed_resistor_load_line'])

    def test_capacitive_threshold_crossing_and_bad_kcl_are_not_qualified(self):
        row = self.state()
        row['drain']['DisplacementCurrent'] = STOP
        row['drain']['e_current_A_per_um'] -= STOP
        self.assertFalse(check_state(row)['checks']['negligible_displacement'])
        row = self.state()
        row['source']['total_current_A_per_um'] *= 0.9
        self.assertFalse(check_state(row)['checks']['terminal_kcl'])

    def test_cancelled_carriers_in_bad_restart_do_not_prove_stop(self):
        row = self.state()
        row['drain'].update(total_current_A_per_um=0, e_current_A_per_um=5e8,
                            h_current_A_per_um=-5e8, inner_voltage_V=row['drain']['outer_voltage_V'])
        row['source']['total_current_A_per_um'] = 0
        self.assertFalse(check_state(row)['checks']['declared_stop'])
        self.assertEqual('fail', check_state(row)['status'])


if __name__ == '__main__':
    unittest.main()
