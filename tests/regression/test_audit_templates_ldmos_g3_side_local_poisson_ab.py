import unittest

from scripts.audit_templates_ldmos_g3_side_local_poisson_ab import (
    material_permittivities,
    region_averagebox_poisson_profile,
)


class SideLocalPoissonAbTest(unittest.TestCase):
    def test_material_weighted_interface_couple_reproduces_region_sum(self) -> None:
        mesh = {
            "nodes": [
                {"id": 0, "x": 0.0, "y": 0.0},
                {"id": 1, "x": 1.0, "y": 0.0},
                {"id": 2, "x": 0.0, "y": 1.0},
                {"id": 3, "x": 0.0, "y": -1.0},
            ],
            "regions": [
                {"id": 0, "material": "Si"},
                {"id": 1, "material": "SiO2"},
            ],
            "triangles": [
                {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                {"id": 1, "region_id": 1, "node_ids": [1, 0, 3]},
            ],
        }
        # coefficient index 2 maps to the 0--1 edge for both orderings.
        coefficients = {0: [0.0, 0.0, 2.0], 1: [0.0, 0.0, 4.0]}
        rows, metadata = region_averagebox_poisson_profile(
            mesh, coefficients, {"si": 12.0, "sio2": 4.0}
        )
        row = next(row for row in rows if row["node0"] == 0 and row["node1"] == 1)
        # Equivalent couple times arithmetic-mean eps is 12*2 + 4*4.
        self.assertAlmostEqual(row["couple_m"] * 1.0e6, 5.0)
        self.assertEqual(metadata["interface_edges"], 1)

    def test_material_contract_is_case_insensitive(self) -> None:
        values = material_permittivities({
            "materials": [{"name": "SiO2", "eps_r": 3.9}]
        })
        self.assertEqual(values, {"sio2": 3.9})


if __name__ == "__main__":
    unittest.main()
