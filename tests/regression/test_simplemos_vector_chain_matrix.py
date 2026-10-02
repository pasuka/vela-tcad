"""Qualification must reject unresolved, nonlinear and unqualified responses."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
import validate_simplemos_vector_chain_matrix as matrix


class ResponseQualification(unittest.TestCase):
    def setUp(self):
        self.gates = {'fd_vs_adjoint_relative': .001, 'even_nonlinear_fraction': .01,
                      'two_amplitude_derivative_relative': .001,
                      'minimum_signal_to_zero_drift': 100, 'zero_drift_dex': 1e-5}
        self.values = {'zero': 1., 'plus_1e13': 1.0001, 'minus_1e13': .9999,
                       'plus_5e12': 1.00005, 'minus_5e12': .99995}

    def check(self, values=None, qualified=True):
        return matrix.score(self.values if values is None else values, .0001, 1., qualified, self.gates)

    def test_odd_linear_response_passes_both_amplitudes(self):
        self.assertTrue(all(r['passed'] for r in self.check()))

    def test_strict_failure_cannot_be_hidden_by_accurate_response(self):
        self.assertFalse(any(r['passed'] for r in self.check(qualified=False)))

    def test_even_response_rejected(self):
        v = dict(self.values, plus_1e13=1.00012, minus_1e13=.99992)
        self.assertFalse(self.check(v)[0]['passed'])

    def test_amplitude_dependence_rejected(self):
        v = dict(self.values, plus_5e12=1.00006, minus_5e12=.99994)
        self.assertFalse(any(r['passed'] for r in self.check(v)))

    def test_restart_drift_rejected(self):
        v = {k: x + 2e-6 for k, x in self.values.items()}
        self.assertFalse(any(r['passed'] for r in self.check(v)))

    def test_missing_and_zero_response_rejected(self):
        for v in ({}, {k: 1. for k in self.values}):
            with self.assertRaises(ValueError):
                self.check(v)


if __name__ == '__main__':
    unittest.main()
