from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "run_cell_first_cross_device_validation.py"
)
SPEC = importlib.util.spec_from_file_location(
    "run_cell_first_cross_device_validation", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class CellFirstCrossDeviceValidationTest(unittest.TestCase):
    def test_legacy_restart_columns_are_normalized_without_changing_values(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "legacy.csv"
            source.write_text(
                "node_id,psi_V,phin_V,phip_V,n_m3,p_m3\n"
                "0,1.25,0.5,-0.25,2e20,3e19\n",
                encoding="utf-8",
            )
            normalized = MODULE.normalized_restart_state(source, root / "normalized.csv")
            rows = MODULE.read_rows(normalized)
            self.assertEqual(
                list(rows[0]),
                ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"],
            )
            self.assertEqual(rows[0]["psi"], "1.25")
            self.assertEqual(rows[0]["electrons_m3"], "2e20")

    def test_current_restart_columns_are_reused(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "state.csv"
            source.write_text(
                "node_id,psi,phin,phip,electrons_m3,holes_m3\n"
                "0,0,0,0,1e20,1e20\n",
                encoding="utf-8",
            )
            self.assertEqual(
                MODULE.normalized_restart_state(source, root / "unused.csv"), source
            )

    def test_vector_metrics_use_only_requested_active_nodes(self) -> None:
        reference = {0: (10.0, 0.0), 1: (1.0e-10, 0.0), 2: (5.0, 0.0)}
        actual = [(10.0, 0.0, 0.0), (1.0, 0.0, 0.0), (10.0, 0.0, 0.0)]
        metrics = MODULE.vector_metrics(reference, actual, {0, 1, 2}, 1.0e-6)
        self.assertEqual(metrics["common_node_count"], 3)
        self.assertEqual(metrics["selected_node_count"], 2)
        self.assertAlmostEqual(
            metrics["max_log10_magnitude_error_decade"], MODULE.math.log10(2.0)
        )

    def test_assessment_keeps_small_p95_tradeoff_but_rejects_large_one(self) -> None:
        def row(recovery: str, p95: float, rmse: float) -> dict[str, object]:
            return {
                "device": "fixture",
                "case": "bias",
                "carrier": "electron",
                "recovery": recovery,
                "p95_log10_magnitude_error_decade": p95,
                "normalized_vector_rmse": rmse,
            }

        small = MODULE.add_ab_assessment(
            [row("direct", 0.50, 0.40), row("cell_first", 0.519, 0.39)]
        )[0]
        self.assertTrue(small["p95_not_materially_worse"])
        large = MODULE.add_ab_assessment(
            [row("direct", 0.50, 0.40), row("cell_first", 0.521, 0.39)]
        )[0]
        self.assertFalse(large["p95_not_materially_worse"])

if __name__ == "__main__":
    unittest.main()
