"""Observer invariants; no simulation or optional HDF5 dependency."""
import importlib.util
from pathlib import Path
import unittest

p = Path(__file__).resolve().parents[2]/"scripts/analysis/slot_shadow_convergence.py"
spec = importlib.util.spec_from_file_location("shadow", p)
s = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)


class ReplayTests(unittest.TestCase):
    def test_legacy_gate_retains_report_and_enforce_semantics(self):
        for enforced in (False, True):
            for satisfied in (False, True):
                q = s.legacy_gate({"enforced": enforced, "satisfied": satisfied})
                self.assertEqual(q["status"], "pass" if not enforced or satisfied else "fail")
                self.assertEqual(q["satisfied"], satisfied)
        self.assertEqual(s.legacy_gate({})["status"], "not_evaluated")

    def test_unpacked_reference_increment_needs_no_packed_scale(self):
        fields = {"psi": [5.], "phin": [5.], "phip": [0.],
                  "electron_qf_reference_V": [5.],
                  "electron_qf_increment_V": [1e-29]}
        state = s.compose_state({}, fields)
        self.assertEqual(state[1]-s.exact(5.), s.exact(1e-29))
        self.assertEqual(state[2], 0)

    def test_packed_low_is_scaled_before_reference_is_added(self):
        fields = {"packed_psi": [2.], "packed_psi_low": [1e-29],
                  "packed_electron_qf_increment": [0.],
                  "packed_electron_qf_increment_low": [1e-29],
                  "packed_hole_qf_increment": [0.],
                  "electron_qf_reference_V": [5.],
                  "hole_qf_reference_V": [0.]}
        state = s.compose_state({"packed_potential_scale_V": .025}, fields)
        self.assertEqual(state[0], (s.exact(2.)+s.exact(1e-29))*s.exact(.025))
        self.assertEqual(state[1]-s.exact(5.), s.exact(1e-29)*s.exact(.025))

    def test_reference_delta_retains_below_high_ulp(self):
        hi, lo = s.exact(5.), s.exact(1.42e-29)
        self.assertEqual(float(hi+lo), 5.)
        self.assertEqual((hi+lo)-hi, lo)

    def test_common_mode_eta_is_zero_but_psi_is_not(self):
        delta = [s.exact(1e-30)]*3
        q = s.eta_components(delta, .025)
        self.assertNotEqual(q[0][0], 0)
        self.assertEqual(q[1][0], 0)
        self.assertEqual(q[2][0], 0)

    def test_large_direction_projected_to_zero_remains_visible(self):
        q = s.transition([s.exact(100.)], [s.exact(0.)], [True])
        self.assertEqual(q["zeroed_components"], 1)
        self.assertEqual(q["max_relative_reduction"], 1.)
        self.assertEqual(q["status"], "not_calibrated")

    def test_constraint_does_not_enter_update_max(self):
        q = s.summarize([1e10, 1e-9], [False, True], [{"Si"}, {"PolySilicon"}])
        self.assertEqual(q["node"], 1)
        self.assertEqual(q["value"], 1e-9)

    def test_missing_and_uncalibrated_are_not_success(self):
        self.assertEqual(s.combine(["pass","not_evaluated"]), "not_evaluated")
        self.assertEqual(s.combine(["not_evaluated","not_calibrated"]), "not_calibrated")
        self.assertEqual(s.combine(["fail","not_calibrated"]), "fail")


if __name__ == "__main__":
    unittest.main()
