import unittest

import numpy as np

from scripts.audit_templates_ldmos_g3_residual_scaling import (
    frozen_kernel_sensitivities,
    regression_report,
    validate_state_specs,
)


class TemplatesLdmosG3ResidualScalingTests(unittest.TestCase):
    def test_requires_exact_eight_bias_contract(self) -> None:
        with self.assertRaisesRegex(ValueError, "exactly eight"):
            validate_state_specs([])

    def test_frozen_kernel_classifies_saturated_endpoint(self) -> None:
        row = {
            "edge_id": "4",
            "length_m": "1e-8",
            "couple_m": "1e-8",
            "electron_mobility_m2_V_s": "0.1",
            "electron_generalized_einstein_factor": "1",
            "electron_bernoulli_argument": "20",
            "electron_density0_m3": "1e20",
            "electron_density1_m3": "1e20",
            "electron_eta0": "-10",
            "electron_eta1": "-10",
            "electron_bernoulli_plus": str(20.0 / np.expm1(20.0)),
            "electron_bernoulli_minus": str(20.0 / (1.0 - np.exp(-20.0))),
        }
        d0, d1 = frozen_kernel_sensitivities(row, 1.0, 300.0)
        self.assertGreater(abs(d0), abs(d1) * 1.0e7)
        self.assertLess(d0, 0.0)
        self.assertGreater(d1, 0.0)

    def test_collinear_scaling_is_downgraded(self) -> None:
        rows = []
        for index in range(8):
            driver = index + 1.0
            rows.append({
                "bias_V": 0.1 * driver,
                "n_interface_geomean_m3": float(np.exp(driver)),
                "Id_sentaurus_A_per_um": float(np.exp(2.0 * driver)),
                "hotspot_grad_phin_rms_V_per_m": float(np.exp(3.0 * driver)),
                "response": float(np.exp(0.5 * driver)),
            })
        report = regression_report(rows, "response")
        self.assertTrue(report["collinearity_exceeds_preregistered_guard"])
        self.assertEqual(report["evidence_role"], "descriptive_only")
        self.assertGreater(max(report["vif"].values()), 10.0)


if __name__ == "__main__":
    unittest.main()
