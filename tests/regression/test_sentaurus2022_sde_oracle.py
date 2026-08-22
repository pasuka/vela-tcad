#!/usr/bin/env python3
"""Integrity checks for the frozen licensed Sentaurus 2022 SDE oracle."""

from __future__ import annotations

import hashlib
import json
import sys
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from sentaurus_device_ir import parse_sde_device_ir  # noqa: E402
from run_sentaurus2022_sde_fixture_vm import portable_ir_provenance  # noqa: E402


ORACLE = REPO / "reference_tcad" / "pn2d_sentaurus2022_oracle"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Sentaurus2022SdeOracleTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads(
            (ORACLE / "manifest.json").read_text(encoding="utf-8")
        )

    def test_manifest_seals_every_generated_input_and_artifact(self) -> None:
        self.assertEqual(
            self.manifest["schema"], "vela.sentaurus2022_sde_oracle.v1"
        )
        self.assertFalse(self.manifest["dry_run"])
        self.assertEqual(self.manifest["exit_code"], 0)
        self.assertIn("Version T-2022.03-SP2", self.manifest["sde_version"])
        self.assertEqual(
            self.manifest["mesh_summary"],
            {"points": 1943, "elements": 3680, "regions": 3},
        )
        for relative, metadata in self.manifest["files"].items():
            with self.subTest(file=relative):
                path = ORACLE / relative
                self.assertTrue(path.is_file())
                self.assertEqual(path.stat().st_size, metadata["bytes"])
                self.assertEqual(sha256_file(path), metadata["sha256"])

    def test_evaluated_log_is_flattened_and_matches_golden_ir(self) -> None:
        evaluated = ORACLE / "artifacts" / self.manifest["evaluated_logs"][0]
        text = evaluated.read_text(encoding="utf-8")
        self.assertIn('(sdegeo:create-rectangle', text)
        self.assertIn('(sde:build-mesh "snmesh" "" "pn2d")', text)
        self.assertNotIn('(load "pn2d_sde.cmd")', text)
        actual = parse_sde_device_ir(
            evaluated, frontend="sentaurus2022-evaluated"
        )
        actual = portable_ir_provenance(actual, evaluated)
        expected = json.loads(
            (ORACLE / self.manifest["expected_device_ir"])
            .read_text(encoding="utf-8")
        )
        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
