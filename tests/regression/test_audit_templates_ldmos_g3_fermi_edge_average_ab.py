import math
import unittest

from scripts.audit_templates_ldmos_g3_fermi_edge_average_ab import (
    KB_J_K,
    Q_C,
    bernoulli,
    candidate_flux,
    factor_variants,
    fermi_half,
    fermi_half_derivative,
    inverse_fermi_half,
    local_generalized_factor,
)


class FermiEdgeAverageAbTest(unittest.TestCase):
    def test_bednarczyk_functions_and_inverse(self) -> None:
        self.assertAlmostEqual(fermi_half(0.0), 0.765147, delta=0.004)
        self.assertGreater(fermi_half_derivative(4.0), 0.0)
        for eta in (-12.0, -1.0, 0.0, 4.0):
            self.assertAlmostEqual(inverse_fermi_half(fermi_half(eta)), eta, places=10)

    def test_local_factor_has_boltzmann_limit(self) -> None:
        self.assertAlmostEqual(local_generalized_factor(-40.0), 1.0, places=8)
        self.assertGreater(local_generalized_factor(4.0), 1.0)

    def test_variants_are_symmetric_under_edge_reversal(self) -> None:
        row = {
            "electron_eta0": "-1.0",
            "electron_eta1": "3.0",
            "electron_density0_m3": str(fermi_half(-1.0)),
            "electron_density1_m3": str(fermi_half(3.0)),
            "electron_generalized_einstein_factor": "1.3",
        }
        reverse = dict(row)
        reverse["electron_eta0"], reverse["electron_eta1"] = (
            row["electron_eta1"], row["electron_eta0"]
        )
        reverse["electron_density0_m3"], reverse["electron_density1_m3"] = (
            row["electron_density1_m3"], row["electron_density0_m3"]
        )
        forward = factor_variants(row)
        backward = factor_variants(reverse)
        self.assertEqual(forward.keys(), backward.keys())
        for name in forward:
            self.assertAlmostEqual(forward[name], backward[name], places=14)

    def test_bernoulli_is_stable_and_obeys_reversal_identity(self) -> None:
        for value in (-100.0, -1.0e-8, 0.0, 1.0e-8, 2.0, 100.0):
            self.assertTrue(math.isfinite(bernoulli(value)))
            self.assertAlmostEqual(
                bernoulli(-value) - bernoulli(value), value,
                delta=1.0e-12 * max(1.0, abs(value)),
            )

    def test_production_factor_reconstructs_recorded_flux(self) -> None:
        factor = 1.03
        thermal_voltage = KB_J_K * 300.0 / Q_C
        row = {
            "electron_flux": "2.5",
            "electron_generalized_einstein_factor": str(factor),
            "electron_bernoulli_argument": "0.2",
            "electron_quasi_fermi_argument": "-0.03",
            "electron_drift_potential_V": str(0.2 * thermal_voltage * factor),
            "phin0_V": "0.0",
            "phin1_V": str(-0.03 * thermal_voltage * factor),
        }
        self.assertAlmostEqual(
            candidate_flux(row, factor, thermal_voltage), 2.5, places=14
        )


if __name__ == "__main__":
    unittest.main()
