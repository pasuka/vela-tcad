from __future__ import annotations

import importlib.util
import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_script(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(
        name, ROOT / "scripts" / filename)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M9 = load_script("simplemos_m9_refdens_scan", "run_simplemos_m9_refdens_scan.py")
CONTRACT_PATH = (
    ROOT / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_m9_hfs_diagnostics_contract_v1.json"
)
PARAMETER_AUDIT = (
    ROOT / "build-release" / "reference_tcad" / "simplemos_sentaurus2022"
    / "m9_hfs_diagnostics" / "parameter_audit" / "hfs_parameter_audit.json"
)
EDGE_REPORT = (
    ROOT / "build-release" / "reference_tcad" / "simplemos_sentaurus2022"
    / "m9_hfs_diagnostics" / "vela_edge_probes" / "edge_hfs_report.json"
)
EVIDENCE_PATH = (
    ROOT / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_m9_hfs_diagnostics_evidence.json"
)


class SimpleMosM9HfsDiagnosticsTest(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = M9.read_json(CONTRACT_PATH)
        M9.validate_contract(self.contract)

    def test_contract_freezes_read_only_log_scan_and_baseline(self) -> None:
        self.assertFalse(self.contract["default_model_policy"][
            "modify_vela_default_model"])
        self.assertTrue(self.contract["default_model_policy"][
            "diagnostics_are_read_only"])
        self.assertEqual(16, self.contract["default_model_policy"][
            "required_case_count"])
        self.assertEqual(816, self.contract["default_model_policy"][
            "required_direct_bias_points"])
        self.assertEqual(
            [0.0, 1.0e2, 1.0e4, 1.0e6, 1.0e8, 1.0e10],
            [float(item["refdens_cm3"]) for item in self.contract["variants"]],
        )
        guard = M9.validate_baseline_guard(self.contract)
        self.assertEqual("unchanged", guard["status"])
        self.assertEqual(16, guard["passing_cases"])
        self.assertEqual(816, guard["direct_bias_points"])

    def test_parameter_audit_matches_300K_core_formula(self) -> None:
        if not PARAMETER_AUDIT.is_file():
            self.skipTest("M9 parameter audit has not been generated")
        audit = json.loads(PARAMETER_AUDIT.read_text(encoding="utf-8"))
        self.assertEqual("complete", audit["status"])
        self.assertTrue(audit["core_formula_equivalent_at_300K"])
        self.assertEqual(107000.0, audit["vela_defaults"]["electron"][
            "saturation_velocity_m_s"])
        self.assertEqual(83700.0, audit["vela_defaults"]["hole"][
            "saturation_velocity_m_s"])
        self.assertEqual(1.109, audit["vela_defaults"]["electron"]["beta"])
        self.assertEqual(1.213, audit["vela_defaults"]["hole"]["beta"])
        self.assertFalse(audit["default_model_changed"])

    def test_edge_probe_exposes_exact_effective_limiter(self) -> None:
        if not EDGE_REPORT.is_file():
            self.skipTest("M9 edge diagnostics have not been generated")
        report = json.loads(EDGE_REPORT.read_text(encoding="utf-8"))
        self.assertEqual("complete", report["status"])
        self.assertEqual(16, report["case_count"])
        self.assertEqual(8, report["intermediate_solve_count"])
        self.assertLess(report["maximum_limiter_reconstruction_error"], 1e-12)
        self.assertGreater(report["maximum_mean_mobility_aggregation_error"], 0.0)
        self.assertFalse(report["default_model_changed"])

    def test_frozen_evidence_preserves_baseline_and_artifact_hashes(self) -> None:
        evidence = M9.read_json(EVIDENCE_PATH)
        self.assertEqual("complete", evidence["status"])
        self.assertEqual(24, evidence["execution"]["sentaurus_curves"])
        self.assertEqual(1224, evidence["execution"][
            "total_direct_bias_points"])
        self.assertFalse(evidence["execution"]["default_model_changed"])
        self.assertEqual("unchanged", evidence["baseline_guard"]["status"])
        self.assertEqual(16, evidence["baseline_guard"]["case_count"])
        self.assertEqual(816, evidence["baseline_guard"]["direct_bias_points"])
        for artifact in evidence["artifacts"].values():
            if isinstance(artifact, list):
                for item in artifact:
                    self.assertEqual(item["sha256"], sha256(
                        ROOT / item["path"]))
            else:
                self.assertEqual(artifact["sha256"], sha256(
                    ROOT / artifact["path"]))
        for relative, expected in evidence["implementation_sha256"].items():
            self.assertEqual(expected, sha256(ROOT / relative))


if __name__ == "__main__":
    unittest.main()
