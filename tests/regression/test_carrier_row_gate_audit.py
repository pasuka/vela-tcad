import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from audit_carrier_row_gate import audit_row, audit_rows, normalized_row_conversion, audit_global_closure


class CarrierRowAuditTests(unittest.TestCase):
    def row(self, **kw):
        return dict(dict(node=7, carrier='hole', step=1, active=True,
                         residual=2e-15, flux_abs_sum=1e-12,
                         recombination=3e-14, impact=0.), **kw)

    def test_local_imbalance_survives_small_absolute_residual(self):
        r = audit_row(self.row())
        self.assertAlmostEqual(r['ratio'], .002)
        self.assertTrue(r['failed'])

    def test_common_unit_change_requires_floor_conversion(self):
        original = self.row()
        for factor in (1e-20, 1e20):
            changed = dict(original)
            for key in ('residual', 'flux_abs_sum', 'recombination', 'impact'):
                changed[key] *= factor
            r = audit_row(changed, min_source_scale=1e-18*factor, scale_floor=1e-30*factor)
            self.assertAlmostEqual(r['ratio'], audit_row(original)['ratio'])
            self.assertEqual(r['qualified'], audit_row(original)['qualified'])

    def test_balanced_row_and_threshold_equality(self):
        self.assertFalse(audit_row(self.row(residual=0.))['failed'])
        self.assertFalse(audit_row(self.row(residual=1e-15))['failed'])

    def test_inactive_and_unqualified_are_not_successful_equations(self):
        for row in (self.row(active=False), self.row(recombination=0., impact=0.)):
            result = audit_row(row)
            self.assertFalse(result['qualified'])
            self.assertFalse(result['failed'])
            self.assertGreater(result['ratio'], 1e-3)

    def test_source_cancellation_does_not_remove_denominator(self):
        r = audit_row(self.row(flux_abs_sum=0., recombination=1., impact=-1., residual=1e-4))
        self.assertEqual(r['scale'], 1.)
        self.assertTrue(r['qualified'])

    def test_invalid_evidence_is_rejected(self):
        for row in (self.row(residual=float('nan')), self.row(flux_abs_sum=-1.),
                    self.row(active='unknown'), self.row(flux_abs_sum=0., recombination=0.)):
            with self.assertRaises(ValueError):
                audit_row(row, scale_floor=0.)

    def test_carriers_are_separate(self):
        _, s = audit_rows([self.row(), self.row(carrier='electron', residual=0.)])
        self.assertEqual([x['failed'] for x in s], [0, 1])

    def test_tcad_particle_line_units_include_length_and_current_factors(self):
        si = normalized_row_conversion(3., 7., 1., 1.)
        tcad = normalized_row_conversion(3., 7., 1e4, 1e-6)
        self.assertAlmostEqual(tcad, 100*si)
        # The same physical target scale represented in each unit system agrees.
        self.assertAlmostEqual(normalized_row_conversion(3., 700., 1e4, 1e-6), si)

    def test_wrong_unit_scale_changes_eligibility(self):
        row = self.row(residual=2e-20, flux_abs_sum=1e-18, recombination=2e-20)
        self.assertFalse(audit_row(row)['qualified'])
        for key in ('residual', 'flux_abs_sum', 'recombination', 'impact'):
            row[key] *= 100
        self.assertTrue(audit_row(row)['failed'])

    def test_global_source_floor_is_qualification_not_passed_physics(self):
        r = audit_global_closure(9e-15, 1e-16)
        self.assertFalse(r['qualified'])
        self.assertTrue(r['satisfied'])
        self.assertGreater(r['ratio'], .01)

    def test_small_net_source_can_fail_despite_large_balanced_terminal_currents(self):
        # Individual terminal fluxes 1 and -1+2e-9 nearly cancel.
        r = audit_global_closure(2e-9, 1e-9)
        self.assertTrue(r['qualified'])
        self.assertEqual(r['ratio'], .5)
        self.assertFalse(r['satisfied'])

    def test_global_closure_units_and_sign(self):
        a = audit_global_closure(-1.001, -1.)
        b = audit_global_closure(-1.001e-20, -1e-20, source_floor=1e-34)
        self.assertAlmostEqual(a['ratio'], b['ratio'])
        self.assertTrue(a['satisfied'])
        with self.assertRaises(ValueError):
            normalized_row_conversion(1., 0., 1., 1.)


if __name__ == '__main__':
    unittest.main()
