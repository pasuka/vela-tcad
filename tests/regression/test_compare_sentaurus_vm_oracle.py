#!/usr/bin/env python3
"""Regression coverage for Sentaurus VM oracle comparator."""

from __future__ import annotations

import csv
import importlib.util
import math
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "compare_sentaurus_vm_oracle",
    REPO / "scripts" / "compare_sentaurus_vm_oracle.py",
)
assert SPEC is not None
module = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(module)


def write_curve(path: Path, values: list[tuple[float, float]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["voltage_V", "current_A_per_um"])
        writer.writerows(values)


class CompareSentaurusVmOracleTest(unittest.TestCase):
    def test_compare_pass_and_fail_threshold(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vm_oracle_compare_") as tmp:
            root = Path(tmp)
            vela_csv = root / "vela.csv"
            oracle_csv = root / "oracle.csv"

            write_curve(vela_csv, [(0.0, 0.0), (1.0, 1.2e-6), (2.0, 2.5e-6)])
            write_curve(oracle_csv, [(0.0, 0.0), (1.0, 1.1e-6), (2.0, 2.4e-6)])

            pass_report = module.compare_series(
                vela_csv,
                oracle_csv,
                "voltage_V",
                "current_A_per_um",
                tolerance=2.0e-7,
                x_tolerance=0.0,
            )
            self.assertTrue(pass_report["pass"])

            fail_report = module.compare_series(
                vela_csv,
                oracle_csv,
                "voltage_V",
                "current_A_per_um",
                tolerance=5.0e-8,
                x_tolerance=0.0,
            )
            self.assertFalse(fail_report["pass"])

    def test_compare_series_rejects_x_grid_mismatch(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vm_oracle_compare_") as tmp:
            root = Path(tmp)
            vela_csv = root / "vela.csv"
            oracle_csv = root / "oracle.csv"

            write_curve(vela_csv, [(0.0, 0.0), (1.0, 1.0e-6), (2.0, 2.0e-6)])
            write_curve(oracle_csv, [(0.0, 0.0), (1.2, 1.0e-6), (2.0, 2.0e-6)])

            report = module.compare_series(
                vela_csv,
                oracle_csv,
                "voltage_V",
                "current_A_per_um",
                tolerance=1.0e-6,
                x_tolerance=0.0,
            )
            self.assertFalse(report["pass"])

    def test_compare_series_x_tolerance_controls_result(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vm_oracle_compare_") as tmp:
            root = Path(tmp)
            vela_csv = root / "vela.csv"
            oracle_csv = root / "oracle.csv"
            write_curve(vela_csv, [(0.0, 0.0), (1.0, 1.0e-6)])
            write_curve(oracle_csv, [(0.0, 0.0), (1.01, 1.0e-6)])
            self.assertFalse(module.compare_series(vela_csv, oracle_csv, "voltage_V", "current_A_per_um", 1e-8, 0.005)["pass"])
            self.assertTrue(module.compare_series(vela_csv, oracle_csv, "voltage_V", "current_A_per_um", 1e-8, 0.02)["pass"])

    def test_rejects_invalid_tolerances(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vm_oracle_compare_") as tmp:
            root = Path(tmp)
            vela_csv = root / "vela.csv"
            oracle_csv = root / "oracle.csv"
            write_curve(vela_csv, [(0.0, 0.0), (1.0, 1.0e-6)])
            write_curve(oracle_csv, [(0.0, 0.0), (1.0, 1.0e-6)])
            with self.assertRaises(ValueError):
                module.compare_series(vela_csv, oracle_csv, "voltage_V", "current_A_per_um", -1.0, 0.0)
            with self.assertRaises(ValueError):
                module.compare_series(vela_csv, oracle_csv, "voltage_V", "current_A_per_um", 0.0, math.inf)

    def test_rejects_empty_nonfinite_duplicate_or_unsorted_x(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vm_oracle_compare_") as tmp:
            root = Path(tmp)
            empty = root / "empty.csv"
            with empty.open("w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh)
                writer.writerow(["voltage_V", "current_A_per_um"])
            valid = root / "valid.csv"
            write_curve(valid, [(0.0, 0.0), (1.0, 1.0e-6)])
            with self.assertRaises(ValueError):
                module.compare_series(empty, valid, "voltage_V", "current_A_per_um", 0.0, 0.0)

            nonfinite = root / "nonfinite.csv"
            with nonfinite.open("w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh)
                writer.writerow(["voltage_V", "current_A_per_um"])
                writer.writerow([0.0, "nan"])
                writer.writerow([1.0, 1.0e-6])
            with self.assertRaises(ValueError):
                module.compare_series(nonfinite, valid, "voltage_V", "current_A_per_um", 0.0, 0.0)

            duplicate_x = root / "duplicate_x.csv"
            write_curve(duplicate_x, [(0.0, 0.0), (0.0, 1.0e-6)])
            with self.assertRaises(ValueError):
                module.compare_series(duplicate_x, valid, "voltage_V", "current_A_per_um", 0.0, 0.0)

            unsorted_x = root / "unsorted_x.csv"
            write_curve(unsorted_x, [(1.0, 0.0), (0.5, 1.0e-6)])
            with self.assertRaises(ValueError):
                module.compare_series(unsorted_x, valid, "voltage_V", "current_A_per_um", 0.0, 0.0)

    def test_allows_strictly_decreasing_curves_but_rejects_direction_mismatch(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_vm_oracle_compare_") as tmp:
            root = Path(tmp)
            dec_a = root / "dec_a.csv"
            dec_b = root / "dec_b.csv"
            inc = root / "inc.csv"
            write_curve(dec_a, [(2.0, 2.0e-6), (1.0, 1.0e-6), (0.0, 0.0)])
            write_curve(dec_b, [(2.0, 2.0e-6), (1.0, 1.0e-6), (0.0, 0.0)])
            write_curve(inc, [(0.0, 0.0), (1.0, 1.0e-6), (2.0, 2.0e-6)])
            self.assertTrue(module.compare_series(dec_a, dec_b, "voltage_V", "current_A_per_um", 0.0, 0.0)["pass"])
            with self.assertRaises(ValueError):
                module.compare_series(dec_a, inc, "voltage_V", "current_A_per_um", 0.0, 0.0)


if __name__ == "__main__":
    unittest.main()
