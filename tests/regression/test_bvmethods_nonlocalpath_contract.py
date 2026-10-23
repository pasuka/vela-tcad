import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from bvmethods_nonlocalpath_contract import (Q_C, generation_current_per_um,
    state_identity, pair_balance, weighted_source_difference, indirect_parameters, temperature_identity)


class NonlocalPathContractTests(unittest.TestCase):
    def test_automatic_mass_relation_and_reduced_mass(self):
        p = indirect_parameters(4e14, 1.9e7, 1.124, .037)
        mc, mv, mr = p['m_c_m0'], p['m_v_m0'], p['m_r_m0']
        self.assertAlmostEqual(1/mc - 1/mv, 2.)
        self.assertAlmostEqual(mc*mv/(mc+mv), mr)

    def test_A_changes_coupling_without_mass_change(self):
        a = indirect_parameters(4e14, 1.9e7, 1.124, .037)
        b = indirect_parameters(8e14, 1.9e7, 1.124, .037)
        self.assertEqual(a['m_c_m0'], b['m_c_m0'])
        self.assertEqual(b['Cpath_J2_cm_kg'], 2*a['Cpath_J2_cm_kg'])

    def test_B_mass_scaling_at_fixed_ratio(self):
        a = indirect_parameters(4e14, 1.9e7, 1.124, .037, 1.)
        b = indirect_parameters(4e14, 3.8e7, 1.124, .037, 1.)
        self.assertEqual(a['m_c_m0'], a['m_v_m0'])
        self.assertAlmostEqual(b['m_r_m0']/a['m_r_m0'], 4.)
        self.assertAlmostEqual(b['Cpath_J2_cm_kg']/a['Cpath_J2_cm_kg'], 2**-3.5)

    def test_reference_gap_is_not_a_runtime_bgn_default(self):
        a = indirect_parameters(4e14, 1.9e7, 1.124, .037)
        b = indirect_parameters(4e14, 1.9e7, 1.1, .037)
        self.assertNotEqual(a['m_r_m0'], b['m_r_m0'])
        with self.assertRaises(ValueError): indirect_parameters(4e14, 1.9e7, 1.124, 0.)
        with self.assertRaises(ValueError): indirect_parameters(4e14, 1e10, 1.124, .037)

    def test_constant_sets_are_explicit_and_complete(self):
        p = indirect_parameters(4e14, 1.9e7, 1.124, .037)
        with self.assertRaises(ValueError):
            indirect_parameters(4e14, 1.9e7, 1.124, .037, constants={'q_C': Q_C})
        physical = dict(p['constants']); physical['m0_kg'] *= 2
        other = indirect_parameters(4e14, 1.9e7, 1.124, .037, 1., constants=physical)
        baseline = indirect_parameters(4e14, 1.9e7, 1.124, .037, 1.)
        self.assertAlmostEqual(other['m_r_m0']/baseline['m_r_m0'], .5)

    def test_volume_units_and_thickness(self):
        # 1 cm^-3 s^-1 in 1 cm^2 cross-section and 1 um depth.
        self.assertAlmostEqual(generation_current_per_um([1.], [1e8]) / Q_C, 1e-4)

    def test_signed_recombination_and_cancellation(self):
        self.assertEqual(generation_current_per_um([1e20, -1e20], [2., 2.]), 0.)
        self.assertLess(generation_current_per_um([-1.], [1.]), 0.)

    def test_invalid_source_units_inputs(self):
        for g, a in [([], []), ([1], []), ([1], [-1]), ([math.nan], [1])]:
            with self.assertRaises(ValueError): generation_current_per_um(g, a)

    def test_complete_state_identity(self):
        s = {(0, 1): 1., (1, 3): 0.}
        self.assertTrue(state_identity(s, s)["identical"])
        self.assertFalse(state_identity(s, {**s, (1, 3): 1e-25})["identical"])
        with self.assertRaises(ValueError): state_identity(s, {(0, 1): 1.})

    def test_nonfinite_state_rejected(self):
        with self.assertRaises(ValueError): state_identity({0: math.nan}, {0: math.nan})

    def test_loaded_temperature_is_checked(self):
        self.assertFalse(temperature_identity([300., 300.], 400.)['identical'])
        self.assertFalse(temperature_identity([400., 300.], 400.)['identical'])
        self.assertTrue(temperature_identity([400., 400.], 400.)['identical'])
        self.assertTrue(temperature_identity([400.00000000000006], 400.)['identical'])
        self.assertFalse(temperature_identity([400.000000001], 400.)['identical'])
        with self.assertRaises(ValueError): temperature_identity([], 400.)

    def test_nonlocal_endpoints_and_signed_events(self):
        # Electron and hole need not reside at the same nodes.
        self.assertEqual(pair_balance([0., 2., -1.], [2., -1., 0.], [2., 1.]), 0.)
        self.assertAlmostEqual(pair_balance([1e-30], [0.], [1e-30]), 1.)

    def test_zero_events_cannot_hide_source(self):
        self.assertEqual(pair_balance([0.], [0.], []), 0.)
        with self.assertRaises(ValueError): pair_balance([1e-100], [1e-100], [])
        with self.assertRaises(ValueError): pair_balance([1.], [1.], [-1.])

    def test_source_tail_is_included(self):
        m = weighted_source_difference([1., 0.], [1., 1e-20], [1., 2.])
        self.assertEqual(m["relative_l1"], 2e-20)

    def test_zero_reference_is_not_a_passing_relative_error(self):
        self.assertIsNone(weighted_source_difference([0.], [1e-30], [1.])["relative_l1"])
        self.assertEqual(weighted_source_difference([0.], [0.], [1.])["relative_l1"], 0.)


if __name__ == "__main__": unittest.main()
