from fractions import Fraction as F
from pathlib import Path
import importlib.util
import unittest

spec = importlib.util.spec_from_file_location(
    "matched", Path(__file__).resolve().parents[2]/"scripts/audit_matched_direction.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class MatchedDirection(unittest.TestCase):
    def test_asymmetric_loaded_steps_need_their_own_actions(self):
        base = {0: 16, 1: 8}
        plus = m.bind_increment(base, {0: 19, 1: 9}, "load", "load")
        minus = m.bind_increment(base, {0: 15, 1: 6}, "load", "load")
        matrix = {0: {0: 2, 1: -3}}
        self.assertEqual(m.linear_action(matrix, plus)[0], 3)
        self.assertEqual(m.linear_action(matrix, minus)[0], 4)
        chord = {k: plus[k]-minus[k] for k in plus}
        self.assertEqual(m.linear_action(matrix, chord)[0], -1)

    def test_rounded_away_coordinate_is_zero_not_intended_step(self):
        x = float(2**54)
        self.assertEqual(m.bind_increment({0: x}, {0: x+1}, "a", "a")[0], 0)

    def test_cancellation_preserves_small_response(self):
        result = m.linear_action({0: {0: 10**30, 1: 1, 2: -10**30}},
                                 {0: 1, 1: 1, 2: 1})
        self.assertEqual(result[0], 1)

    def test_reference_and_coordinate_mismatches_rejected(self):
        with self.assertRaises(ValueError):
            m.bind_increment({0: 1}, {0: 2}, "load", "qs")
        with self.assertRaises(ValueError):
            m.bind_increment({0: 1}, {1: 2}, "a", "a")
        with self.assertRaises(ValueError):
            m.linear_action({0: {2: 1}}, {0: 1})

    def test_uncertainty_does_not_manufacture_a_pass(self):
        self.assertEqual(m.compare_response(1, 1)["status"], "pass")
        self.assertEqual(m.compare_response(1, 2)["status"], "fail")
        self.assertEqual(m.compare_response(1, 1, "0.1")["status"], "unresolved")
        self.assertEqual(m.compare_response(0, 0)["status"], "zero")
        self.assertEqual(m.compare_response(1, 0)["status"], "fail")
        self.assertEqual(m.compare_response(0, 1, 2)["status"], "unresolved")

    def test_negative_response_and_exact_threshold(self):
        self.assertEqual(m.compare_response(-1, -1)["status"], "pass")
        self.assertEqual(m.compare_response(-1, 1)["relative_error"], 1)
        self.assertEqual(m.compare_response(99, 101, tolerance="0.01")["status"], "pass")
        self.assertEqual(m.compare_response(99, 101, tolerance="0.009")["status"], "fail")

    def test_nonlinear_remainder_shrinks_without_changing_direction(self):
        errors = [m.compare_response(2*h+h*h, 2*h)["relative_error"]
                  for h in [F(1, 10), F(1, 100), F(1, 1000)]]
        self.assertTrue(errors[0] > errors[1] > errors[2])

    def test_nonfinite_or_negative_bounds_rejected(self):
        for x in [True, "NaN", "Infinity"]:
            with self.assertRaises(ValueError): m.number(x)
        with self.assertRaises(ValueError): m.compare_response(1, 1, -1)

    def test_endpoint_gate_uses_actual_chord_not_two_intended_steps(self):
        self.assertTrue(m.endpoint_intersection(4, 3, 1, 1, -1)["overlap"])
        self.assertFalse(m.endpoint_intersection(4, 3, 1, 1, 6)["overlap"])

    def test_endpoint_gate_zero_action_and_exact_touch(self):
        x = m.endpoint_intersection(0, 2, 1, 1, 0,
                                     absolute_tolerance=1)
        self.assertTrue(x["overlap"])
        self.assertEqual(x["absolute_change_interval"], (0, 4))
        self.assertFalse(m.endpoint_intersection(0, 2, 0, 0, 0)["overlap"])
        with self.assertRaises(ValueError):
            m.endpoint_intersection(0, 1, -1, 1, 1)


if __name__ == "__main__":
    unittest.main()
