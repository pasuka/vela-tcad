import math
import unittest

from scripts.audit_templates_ldmos_g3_interface_edge_state import (
    analyze_edges,
    edge_factors,
    factor_product,
    signed_at_node,
)


def edge(edge_id: int, node0: int, node1: int, flux: float, scale: float = 1.0):
    qf_argument = 0.01 * scale
    return {
        "edge_id": str(edge_id),
        "node0": str(node0),
        "node1": str(node1),
        "x0": "0", "y0": "0", "x1": "1", "y1": "0",
        "length_m": "2", "couple_m": "1",
        "psi0_V": "0.1", "psi1_V": str(0.1 + 0.02 * scale),
        "phin0_V": "0", "phin1_V": str(0.001 * scale),
        "electron_mobility_field_V_m": str(10.0 * scale),
        "electron_mobility_m2_V_s": str(2.0 * scale),
        "electron_generalized_einstein_factor": str(1.2 * scale),
        "electron_bernoulli_plus": str(0.8 * scale),
        "electron_density1_m3": str(5.0 * scale),
        "electron_bernoulli_argument": str(0.3 * scale),
        "electron_quasi_fermi_argument": str(qf_argument),
        "electron_flux": str(flux),
    }


class TemplatesLdmosG3InterfaceEdgeStateTest(unittest.TestCase):
    def test_edge_factor_product_includes_exact_qf_drive(self) -> None:
        row = edge(0, 1, 2, 3.0)
        factors = edge_factors(row)
        self.assertAlmostEqual(factors["qf_drive"], math.expm1(0.01))
        self.assertAlmostEqual(
            factor_product(factors),
            0.5 * 2.0 * 1.2 * 0.8 * 5.0 * math.expm1(0.01),
        )

    def test_signed_at_node_uses_residual_orientation(self) -> None:
        row = edge(0, 1, 2, 3.0)
        self.assertEqual(signed_at_node(row, 1, 3.0), 3.0)
        self.assertEqual(signed_at_node(row, 2, 3.0), -3.0)

    def test_analysis_reconstructs_factor_ratio_and_row_sums(self) -> None:
        low_a = edge(0, 1, 2, 0.0)
        low_b = edge(1, 3, 1, 0.0)
        for row, flux in ((low_a, 2.0), (low_b, 1.0)):
            row["electron_flux"] = str(flux)
        high_a = edge(0, 1, 2, 0.0, 2.0)
        high_b = edge(1, 3, 1, 0.0, 2.0)
        for low, high in ((low_a, high_a), (low_b, high_b)):
            ratio = factor_product(edge_factors(high)) / factor_product(
                edge_factors(low)
            )
            high["electron_flux"] = str(float(low["electron_flux"]) * ratio)
        summary, records = analyze_edges([low_a, low_b], [high_a, high_b], (1,))
        self.assertEqual(summary["unique_incident_edges"], 2)
        self.assertAlmostEqual(summary["nodes"][0]["low_residual"], 1.0)
        self.assertLess(
            summary["maximum_factor_product_flux_ratio_relative_error"], 1.0e-14
        )
        self.assertEqual(len(records), 2)


if __name__ == "__main__":
    unittest.main()
