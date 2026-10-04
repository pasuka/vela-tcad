from __future__ import annotations

import hashlib
import io
import importlib.util
import json
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "run_simplemos_m2_tdr_gate",
    REPO / "scripts" / "run_simplemos_m2_tdr_gate.py",
)
assert SPEC is not None and SPEC.loader is not None
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)


class SimpleMosTdrGateTest(unittest.TestCase):
    def test_hash_mismatch_fails_before_importer(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_simplemos_m2_gate_") as tmp:
            root = Path(tmp)
            tdr = root / "nominal.tdr"
            tdr.write_bytes(b"not-the-frozen-tdr")
            contract = root / "contract.json"
            contract.write_text(json.dumps({
                "nominal_tdr": {"tdr_sha256": "0" * 64},
            }), encoding="utf-8")
            output = root / "output"

            with redirect_stdout(io.StringIO()):
                returncode = gate.main([
                    "--tdr", str(tdr),
                    "--contract", str(contract),
                    "--tdr-importer", str(root / "must-not-run"),
                    "--output-dir", str(output),
                ])

            self.assertEqual(1, returncode)
            report = json.loads(
                (output / "simplemos_m2_tdr_gate_report.json").read_text(encoding="utf-8"))
            self.assertFalse(report["passed"])
            self.assertEqual("source.sha256", report["checks"][0]["code"])
            self.assertFalse(report["checks"][0]["passed"])

    def test_gate_report_combines_hash_and_importer_checks(self) -> None:
        digest = hashlib.sha256(b"same").hexdigest()
        report = gate.gate_report(
            Path("nominal.tdr"),
            Path("contract.json"),
            digest,
            digest,
            {"passed": True, "checks": [{
                "code": "geometry.contacts",
                "passed": True,
                "message": "four contacts",
            }]},
            0,
        )

        self.assertTrue(report["passed"])
        self.assertEqual(
            ["source.sha256", "geometry.contacts"],
            [check["code"] for check in report["checks"]],
        )


if __name__ == "__main__":
    unittest.main()
