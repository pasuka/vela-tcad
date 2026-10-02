from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest

from tests.regression.simplemos_evidence_chain import (
    assert_historical_source_provenance,
)


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = (ROOT / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m12_terminal_sensitivity_contract_v1.json")
EVIDENCE = (ROOT / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m12_terminal_sensitivity_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM12TerminalSensitivityTest(unittest.TestCase):
    def test_contract_is_read_only_and_preregistered(self) -> None:
        contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
        policy = contract["default_model_policy"]
        self.assertFalse(policy["modify_vela_default_model"])
        self.assertFalse(policy["modify_cpp"])
        self.assertTrue(policy["diagnostics_are_read_only"])
        self.assertEqual(16, contract["frozen_state_matrix"]["state_count"])
        self.assertEqual("forbidden",
                         contract["frozen_state_matrix"]["interpolation"])
        self.assertEqual(["source", "drain", "substrate"],
                         contract["contact_cut_definition"]["contacts"])
        self.assertEqual(0.1, contract["h1"]
                         ["contact_to_active_reduction_ratio_threshold"])

    def test_h1_and_exact_cut_identity_are_frozen(self) -> None:
        evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        findings = evidence["findings"]
        self.assertEqual("complete", evidence["status"])
        self.assertEqual(16, evidence["execution"]["state_count"])
        self.assertFalse(evidence["execution"]["cpp_changed"])
        self.assertFalse(evidence["execution"]["default_model_changed"])
        self.assertEqual(16, findings["drain_replay_identity_pass_count"])
        self.assertEqual(16, findings["h1_state_pass_count"])
        self.assertTrue(findings["h1_all_states_pass"])
        self.assertEqual("accepted", evidence["conclusions"]["h1"])
        self.assertLess(findings["maximum_contact_to_active_reduction_ratio"],
                        2.0e-4)
        self.assertLess(findings["maximum_abs_terminal_drive_shift_dex"],
                        1.1e-5)

    def test_error_spectrum_and_conditioning_are_explicit(self) -> None:
        evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        findings = evidence["findings"]
        weak = findings["high_nwell_low_drain_vg_le_0p95"]
        self.assertEqual(80, weak["point_count"])
        self.assertTrue(weak["all_vela_above_sentaurus"])
        self.assertGreater(weak["mean_signed_error_dex"], 0.04)
        self.assertEqual(0, findings[
            "deep_off_any_component_kappa_ge_100_count"])
        maxima = findings["deep_off_contact_kappa_maxima"]
        self.assertLess(maxima["drain"]["kappa_total"], 4.0)

    def test_artifact_and_implementation_hashes_match(self) -> None:
        evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        for artifact in evidence["artifacts"].values():
            self.assertEqual(artifact["sha256"], sha256(ROOT / artifact["path"]))
        for figure in evidence["figures"]:
            self.assertEqual(figure["sha256"], sha256(ROOT / figure["path"]))
        for relative, expected in evidence["implementation_sha256"].items():
            self.assertTrue((ROOT / relative).is_file())
        assert_historical_source_provenance(self, EVIDENCE)


if __name__ == "__main__":
    unittest.main()
