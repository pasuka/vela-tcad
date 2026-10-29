"""Physical balance checks for the bounded PN2D resistor experiment."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from run_pn2d_external_resistor_control import qualify_rows, validate_mobility_basis


class QualificationTests(unittest.TestCase):
    def test_impurity_scattering_rejects_compensated_net_doping(self):
        for basis in ['total_impurity', 'cell_reconstructed_total_impurity']:
            validate_mobility_basis({'solver': {'mobility': {'doping_concentration_basis': basis}}})
        with self.assertRaisesRegex(ValueError, 'Nd \\+ Na'):
            validate_mobility_basis({'solver': {'mobility': {'doping_concentration_basis': 'net_doping'}}})

    def rows(self):
        return [dict(converged='1', carrier_row_violations='0',
                     global_continuity_closure_satisfied='1',
                     current_total_A_per_um='-1e-10',
                     contact_Anode_total_A_per_um='-1e-10',
                     contact_Cathode_total_A_per_um='1e-10',
                     outer_voltage_V=str(v), inner_voltage_V=str(v + 0.001))
                for v in [0, -0.002, -0.004, -0.006, -0.008, -0.01]]

    def test_passive_resistor_drop_and_balanced_contacts(self):
        self.assertTrue(qualify_rows(self.rows(), 'coupled_newton', 1e7)['passed'])

    def test_wrong_current_sign_fails_load_equation(self):
        rows = self.rows()
        for row in rows:
            row['current_total_A_per_um'] = '1e-10'
        self.assertFalse(qualify_rows(rows, 'coupled_newton', 1e7)['passed'])

    def test_kcl_and_device_gates_remain_independent(self):
        for key, value in [('contact_Cathode_total_A_per_um', '2e-10'),
                           ('carrier_row_violations', '1'),
                           ('global_continuity_closure_satisfied', '0')]:
            rows = copy.deepcopy(self.rows())
            rows[-1][key] = value
            self.assertFalse(qualify_rows(rows, 'nested_scalar', 1e7)['passed'])

    def test_missing_requested_endpoint_fails(self):
        self.assertFalse(qualify_rows(self.rows()[:-1], 'coupled_newton', 1e7)['passed'])


if __name__ == '__main__':
    unittest.main()
