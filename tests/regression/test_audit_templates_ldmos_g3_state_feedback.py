import unittest

from scripts.audit_templates_ldmos_g3_state_feedback import (
    attribution,
    merge_state_rows,
)


class TemplatesLdmosG3StateFeedbackAuditTest(unittest.TestCase):
    def test_merge_selects_each_state_family_independently(self) -> None:
        sentaurus = {
            0: {"psi": 1.0, "phin": 2.0, "phip": 3.0,
                "electrons_m3": 4.0, "holes_m3": 5.0}
        }
        vela = {
            0: {"psi": 10.0, "phin": 20.0, "phip": 30.0,
                "electrons_m3": 40.0, "holes_m3": 50.0}
        }
        row = merge_state_rows(sentaurus, vela, ("S", "V", "S"))[0]
        self.assertEqual(row["psi"], 1.0)
        self.assertEqual(row["phin"], 20.0)
        self.assertEqual(row["phip"], 3.0)
        self.assertEqual(row["electrons_m3"], 40.0)
        self.assertEqual(row["holes_m3"], 5.0)

    def test_attribution_separates_operator_and_feedback(self) -> None:
        currents = {
            "VVV": 2.7, "SSS": 1.05, "SVV": 2.6, "VSV": 1.2,
            "VVS": 2.7, "SSV": 1.1, "SVS": 2.6, "VSS": 1.1,
        }
        result = attribution(currents, 1.0)
        self.assertAlmostEqual(result["operator_log10_ratio_dex"], math_log10(1.05))
        self.assertEqual(result["dominant_single_family"], "phin")
        self.assertGreater(result["feedback_amplification_dex"], 0.4)


def math_log10(value: float) -> float:
    # Keep the expected expression visibly independent from implementation data.
    import math
    return math.log10(value)


if __name__ == "__main__":
    unittest.main()
