import unittest

import numpy as np

from scripts.audit_templates_ldmos_g3_interface_charge_volume import (
    build_geometry,
    poisson_flux,
)


class InterfaceChargeVolumeAuditTest(unittest.TestCase):
    def test_material_local_volume_excludes_oxide_share(self) -> None:
        mesh = {
            "nodes": [
                {"id": 0, "x": 0.0, "y": 0.0},
                {"id": 1, "x": 1.0, "y": 0.0},
                {"id": 2, "x": 0.0, "y": 1.0},
                {"id": 3, "x": 0.0, "y": -1.0},
            ],
            "regions": [
                {"id": 0, "material": "Si", "cell_ids": [0]},
                {"id": 1, "material": "SiO2", "cell_ids": [1]},
            ],
            "triangles": [
                {"id": 0, "node_ids": [0, 1, 2]},
                {"id": 1, "node_ids": [1, 0, 3]},
            ],
            "contacts": [],
        }
        geometry = build_geometry(mesh)
        self.assertAlmostEqual(geometry["global_volume"][0], 1.0e-12 / 3.0)
        self.assertAlmostEqual(
            geometry["transport_volume"][0], 0.5e-12 / 3.0
        )
        self.assertEqual(set(geometry["interface_nodes"]), {0, 1})

    def test_poisson_flux_is_conservative(self) -> None:
        mesh = {
            "nodes": [
                {"id": 0, "x": 0.0, "y": 0.0},
                {"id": 1, "x": 1.0, "y": 0.0},
                {"id": 2, "x": 0.0, "y": 1.0},
            ],
            "regions": [{"id": 0, "material": "Si", "cell_ids": [0]}],
            "triangles": [{"id": 0, "node_ids": [0, 1, 2]}],
            "contacts": [],
        }
        flux = poisson_flux(np.array([0.0, 0.1, -0.2]), build_geometry(mesh))
        self.assertAlmostEqual(float(np.sum(flux)), 0.0, places=24)


if __name__ == "__main__":
    unittest.main()
