"""M79 gates reject invalid observers even if paired current trends look good."""
import copy
from pathlib import Path
import sys
import unittest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m79_terminal_response_qualification as m79


class QualificationGates(unittest.TestCase):
    def setUp(self):
        self.thresholds = {"current_replay_relative": 1e-8, "adjoint_relative_residual": 1e-10,
                           "jvp_relative": 1e-3, "current_fd_relative": 1e-3,
                           "duality_relative": 1e-8, "linear_solve_relative": 1e-8,
                           "continuity_source_norm": 0, "direct_current_term": 0}
        self.status = {"current_A_per_um": 1e-6, "state_derivative_norm": 1e-3,
                       "adjoint_relative_residual": 1e-14,
                       "electron_volume_response": {
                           "finite_difference_checks": [{"current_analytic": 1e-6,
                                                          "current_relative_error": 1e-5,
                                                          "jvp_relative_error": 1e-5}],
                           "duality_relative_error": 1e-12, "linear_solve_relative_residual": 1e-12,
                           "continuity_source_norm": 0, "parameter_direct_current_A_per_um": 0}}
        self.functional = {"current_A_per_um": 1e-6, "contact_current_extractor_A_per_um": 1e-6}

    def gates(self, status=None, functional=None):
        return m79.qualify_endpoint(status or self.status, functional or self.functional,
                                   1e-6, self.thresholds)[0]

    def test_valid_observer(self):
        self.assertTrue(all(self.gates().values()))

    def test_extractor_mismatch_is_not_hidden_by_matching_functional(self):
        functional = dict(self.functional, contact_current_extractor_A_per_um=1.001e-6)
        self.assertFalse(self.gates(functional=functional)["current_identity"])

    def test_transpose_closure_does_not_qualify_wrong_jacobian(self):
        status = copy.deepcopy(self.status)
        status["electron_volume_response"]["finite_difference_checks"][0]["jvp_relative_error"] = 0.2
        gates = self.gates(status=status)
        self.assertTrue(gates["adjoint"])
        self.assertFalse(gates["jvp"])

    def test_continuity_source_rejects_non_single_axis(self):
        status = copy.deepcopy(self.status)
        status["electron_volume_response"]["continuity_source_norm"] = 1e-20
        self.assertFalse(self.gates(status=status)["single_axis"])

    def test_no_measurable_current_direction_cannot_pass(self):
        status = copy.deepcopy(self.status)
        status["electron_volume_response"]["finite_difference_checks"][0]["current_analytic"] = 1e-30
        self.assertFalse(self.gates(status=status)["current_derivative"])


if __name__ == "__main__":
    unittest.main()
