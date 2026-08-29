import csv
import hashlib
import json
import math
import tempfile
import unittest
from pathlib import Path

from scripts.audit_templates_ldmos_g3_idvg_shift_kcl import (
    Q_C,
    cell_qf_gradient_drain_current,
    drain_density_reconstruction,
    drain_sg_cut,
    sentaurus_drain_vector_current,
)
from scripts.templates_ldmos_contracts import validate_document


class TemplatesLdmosG3ShiftKclAuditTests(unittest.TestCase):
    def test_candidate_known_difference_remains_draft_and_schema_valid(self) -> None:
        path = Path(
            "reference_tcad/templates_ldmos_sentaurus2022/"
            "known_difference_ledger.json"
        )
        document = json.loads(path.read_text(encoding="utf-8"))
        validate_document(document)
        entry = document["entries"][0]
        self.assertEqual(entry["classification"], "candidate_baseline_difference")
        self.assertEqual(entry["status"], "open")
        self.assertEqual(entry["approval"]["status"], "draft")
        for evidence in entry["evidence"]:
            digest = hashlib.sha256(Path(evidence["path"]).read_bytes()).hexdigest()
            self.assertEqual(digest, evidence["sha256"])

    def test_drain_cut_and_density_reconstruction_close_exact_state(self) -> None:
        rows = [{
            "node0": "0",
            "node1": "1",
            "electron_particle_line_flux_per_m_s": "2.0e6",
            "hole_particle_line_flux_per_m_s": "5.0e5",
            "electron_density0_m3": "1.0e20",
            "electron_density1_m3": "2.0e20",
            "phin0_V": "0.0",
            "phin1_V": "0.01",
        }]
        cut = drain_sg_cut(rows, {0})
        self.assertAlmostEqual(cut["electron_A_per_um"], -2.0 * Q_C)
        self.assertAlmostEqual(
            cut["total_A_per_um"], -1.5 * Q_C
        )
        state = {
            0: {"electrons_m3": 1.0e20},
            1: {"electrons_m3": 2.0e20},
        }
        reconstruction = drain_density_reconstruction(rows, state, {0})
        self.assertEqual(reconstruction["crossing_edges"], 1)
        self.assertEqual(
            reconstruction[
                "contact_reconstructed_minus_supplied_electron_density_dex"
            ]["maximum_abs"],
            0.0,
        )
        self.assertEqual(
            reconstruction[
                "interior_reconstructed_minus_supplied_electron_density_dex"
            ]["maximum_abs"],
            0.0,
        )

    def test_native_vector_and_cell_qf_gradient_use_same_oriented_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fields = root / "fields"
            fields.mkdir()
            mesh = {
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 0.0, "y": 1.0},
                ],
                "triangles": [
                    {"id": 0, "region_id": 0, "node_ids": [0, 1, 2]}
                ],
                "regions": [{"id": 0, "material": "Si"}],
                "contacts": [{
                    "id": 0,
                    "name": "drain",
                    "region_id": 0,
                    "node_ids": [0, 1],
                    "edge_node_ids": [[0, 1]],
                }],
            }
            mesh_path = root / "mesh.json"
            mesh_path.write_text(json.dumps(mesh), encoding="utf-8")

            density_cm3 = 1.0e10
            mobility_cm2 = 1.0e3
            jy_A_cm2 = Q_C * mobility_cm2 * density_cm3 * 1.0e4

            def write_scalar(name: str, values: list[float]) -> None:
                with (fields / f"{name}_region0.csv").open(
                    "w", newline="", encoding="utf-8"
                ) as stream:
                    writer = csv.writer(stream)
                    writer.writerow(["node_id", "component0"])
                    writer.writerows(enumerate(values))

            write_scalar("eQuasiFermiPotential", [0.0, 0.0, 1.0])
            write_scalar("eDensity", [density_cm3] * 3)
            write_scalar("eMobility", [mobility_cm2] * 3)
            with (fields / "eCurrentDensity_region0.csv").open(
                "w", newline="", encoding="utf-8"
            ) as stream:
                writer = csv.writer(stream)
                writer.writerow(["node_id", "component0", "component1"])
                writer.writerows((node, 0.0, jy_A_cm2) for node in range(3))

            expected = -jy_A_cm2 * 1.0e-8
            native = sentaurus_drain_vector_current(root, mesh_path)
            reconstructed = cell_qf_gradient_drain_current(root, mesh_path)
            self.assertEqual(native["edges"], 1)
            self.assertEqual(reconstructed["edges"], 1)
            self.assertTrue(math.isclose(
                native["oriented_electron_A_per_um"], expected, rel_tol=1.0e-14
            ))
            self.assertTrue(math.isclose(
                reconstructed["oriented_electron_A_per_um"], expected,
                rel_tol=1.0e-14,
            ))


if __name__ == "__main__":
    unittest.main()
