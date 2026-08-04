#!/usr/bin/env python3
"""Regression coverage for Sentaurus VM oracle comparator."""

from __future__ import annotations

import csv
import importlib.util
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

            pass_report = module.compare(
                vela_csv,
                oracle_csv,
                "voltage_V",
                "current_A_per_um",
                tolerance=2.0e-7,
            )
            self.assertTrue(pass_report["pass"])

            fail_report = module.compare(
                vela_csv,
                oracle_csv,
                "voltage_V",
                "current_A_per_um",
                tolerance=5.0e-8,
            )
            self.assertFalse(fail_report["pass"])


if __name__ == "__main__":
    unittest.main()

