#!/usr/bin/env python3
"""Integrity checks for the frozen licensed Sentaurus 2022 SDE oracle."""

from __future__ import annotations

import hashlib
import json
import os
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

    def test_manifest_seals_public_files_and_available_local_evidence(self) -> None:
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
        local_only = set(self.manifest["local_validation_files"])
        self.assertEqual(local_only, {
            "artifacts/pn2d_bnd.tdr", "artifacts/pn2d_msh.tdr",
            "artifacts/sde.exit_code.txt", "artifacts/sde.stderr.txt",
            "artifacts/sde.stdout.txt", "artifacts/sde.version.txt",
            "artifacts/vela_sde_evaluated_wrapper.log.eval",
        })
        self.assertTrue(local_only.issubset(self.manifest["files"]))
        for relative, metadata in self.manifest["files"].items():
            with self.subTest(file=relative):
                path = ORACLE / relative
                if relative in local_only and not path.is_file():
                    if os.environ.get("VELA_REQUIRE_LOCAL_SDE_ORACLE") == "1":
                        self.fail(f"required local oracle evidence missing: {relative}")
                    continue
                self.assertTrue(path.is_file())
                self.assertEqual(path.stat().st_size, metadata["bytes"])
                self.assertEqual(sha256_file(path), metadata["sha256"])

    def test_owned_source_matches_frozen_semantics_without_raw_logs(self) -> None:
        actual = parse_sde_device_ir(ORACLE / "input" / "pn2d_sde.cmd")
        expected = json.loads(
            (ORACLE / self.manifest["expected_device_ir"]).read_text(encoding="utf-8")
        )

        def semantics(value):
            if isinstance(value, dict):
                return {key: semantics(item) for key, item in value.items()
                        if key not in {"source", "location"}}
            if isinstance(value, list):
                return [semantics(item) for item in value]
            return value

        self.assertEqual(actual.pop("metadata_only"), [])
        wrapper_control = expected.pop("metadata_only")
        self.assertEqual(len(wrapper_control), 1)
        self.assertEqual(wrapper_control[0]["command"], "set!")
        self.assertEqual(wrapper_control[0]["text"], "(set! evaluate-log-file #t)")
        self.assertEqual(semantics(actual), semantics(expected))

    def test_evaluated_log_is_flattened_and_matches_golden_ir(self) -> None:
        evaluated = ORACLE / "artifacts" / self.manifest["evaluated_logs"][0]
        if not evaluated.is_file():
            if os.environ.get("VELA_REQUIRE_LOCAL_SDE_ORACLE") == "1":
                self.fail("required local evaluated log missing")
            self.skipTest("raw evaluated log is local validation evidence")
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
