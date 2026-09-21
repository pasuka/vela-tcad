import unittest

from scripts.audit_templates_ldmos_g3_element_hfs import caughey_thomas


class TemplatesLdmosG3ElementHfsAuditTest(unittest.TestCase):
    def test_caughey_thomas_limits(self) -> None:
        low_field = 1417.0
        vsat = 1.07e7
        beta = 1.109

        self.assertAlmostEqual(
            caughey_thomas(low_field, 0.0, vsat, beta), low_field
        )
        high_field = 1.0e8
        mobility = caughey_thomas(low_field, high_field, vsat, beta)
        self.assertAlmostEqual(mobility * high_field, vsat, delta=0.01 * vsat)


if __name__ == "__main__":
    unittest.main()
