import unittest

from scripts.run_templates_ldmos_g3_contact_hfs_replay import (
    contact_hfs_variants,
)


class TemplatesLdmosG3ContactHfsReplayTest(unittest.TestCase):
    def test_variants_change_only_contact_hfs_contract(self) -> None:
        baseline = {
            "solver": {
                "mobility": {
                    "model": "constant_field",
                    "high_field_driving_force": "quasi_fermi_gradient",
                    "contact_electric_field_fallback": True,
                    "contact_electric_field_fallback_scope": "contact_node_cell",
                    "contact_electric_field_fallback_mode": "cell_gradient_magnitude",
                }
            },
            "mesh_geometry": {
                "carrier_transport_couple_profile":
                    "templates_ldmos_external_averagebox"
            },
        }

        hfs_off, hfs_on = contact_hfs_variants(baseline)

        self.assertFalse(
            hfs_off["solver"]["mobility"]["contact_electric_field_fallback"]
        )
        self.assertNotIn(
            "contact_electric_field_fallback_scope",
            hfs_off["solver"]["mobility"],
        )
        self.assertTrue(
            hfs_on["solver"]["mobility"]["contact_electric_field_fallback"]
        )
        self.assertEqual(
            hfs_off["mesh_geometry"], hfs_on["mesh_geometry"]
        )
        self.assertTrue(
            baseline["solver"]["mobility"]["contact_electric_field_fallback"]
        )


if __name__ == "__main__":
    unittest.main()
