import csv
import json
import tempfile
import unittest
from pathlib import Path

from scripts.audit_templates_ldmos_g3_qf_single_mode_perturbation import (
    loaded_pair_validation,
    nonlinearity,
    one_ring,
    vector_metrics,
)


def write_field(root: Path, name: str, values: dict[int, float]) -> None:
    path = root / "fields" / f"{name}_region0.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["node_id", "component0"])
        writer.writeheader()
        for node, value in values.items():
            writer.writerow({"node_id": node, "component0": value})


class QfSingleModePerturbationAuditTest(unittest.TestCase):
    def test_one_ring_is_restricted_to_silicon_triangles(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            mesh = Path(root_text) / "mesh.json"
            mesh.write_text(json.dumps({"triangles": [
                {"region_id": 0, "node_ids": [1, 2, 3]},
                {"region_id": 0, "node_ids": [1, 3, 4]},
                {"region_id": 1, "node_ids": [1, 4, 5]},
            ]}), encoding="utf-8")
            self.assertEqual(one_ring(mesh, 1), [1, 2, 3, 4])

    def test_vector_metrics_recovers_opposite_sign_scalar(self) -> None:
        metrics = vector_metrics([-34.0, 68.0], [1.0, -2.0])
        self.assertAlmostEqual(metrics["cosine"], -1.0)
        self.assertAlmostEqual(
            metrics["signed_sentaurus_over_vela_least_squares"], -34.0
        )
        self.assertAlmostEqual(metrics["relative_l2_after_scalar_fit"], 0.0)

    def test_loaded_pair_has_only_requested_qf_delta(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            root = Path(root_text)
            minus = root / "minus"
            plus = root / "plus"
            for export, target_phin in ((minus, 0.2 - 1.0e-6), (plus, 0.2 + 1.0e-6)):
                write_field(export, "ElectrostaticPotential", {1: 0.1, 2: 0.4})
                write_field(export, "eQuasiFermiPotential", {1: target_phin, 2: 0.5})
                write_field(export, "hQuasiFermiPotential", {1: 0.3, 2: 0.6})
            result = loaded_pair_validation(minus, plus, 1, 1.0e-6)
            self.assertAlmostEqual(result["target_phin_delta_error_V"], 0.0)
            self.assertEqual(
                result["maximum_abs_unintended_pair_delta_V"],
                {"psi": 0.0, "phip": 0.0, "phin_excluding_target": 0.0},
            )

    def test_nonlinearity_uses_even_over_odd_central_difference(self) -> None:
        result = nonlinearity(
            {1: 0.9, 2: 2.2},
            {1: 1.0, 2: 2.0},
            {1: 1.1, 2: 1.8},
            [1, 2],
        )
        self.assertAlmostEqual(result["even_second_difference_l2"], 0.0)
        self.assertAlmostEqual(result["even_over_odd_l2"], 0.0)


if __name__ == "__main__":
    unittest.main()
