from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = (ROOT / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m10_fixed_state_replay_contract_v1.json")
EVIDENCE = (ROOT / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m10_fixed_state_replay_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM10FixedStateReplayTest(unittest.TestCase):
    def test_contract_is_read_only_and_exact(self) -> None:
        contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
        self.assertFalse(contract["default_model_policy"]["modify_vela_default_model"])
        self.assertTrue(contract["default_model_policy"]["diagnostics_are_read_only"])
        self.assertEqual(16, contract["bias_matrix"]["state_count"])
        self.assertEqual("forbidden", contract["bias_matrix"]["interpolation"])

    def test_frozen_evidence_is_complete_and_integral(self) -> None:
        evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        self.assertEqual("complete", evidence["status"])
        self.assertEqual(16, evidence["execution"]["state_count"])
        self.assertFalse(evidence["execution"]["default_model_changed"])
        improvement = evidence["findings"]["mobility_improvement"]
        self.assertTrue(improvement["sentaurus_drive_improves_all_states_at_p95"])
        self.assertGreater(improvement["p95_reduction_fraction"], 0.25)
        self.assertLess(evidence["findings"]
                        ["mapping_sensitivity_max_p95_difference_dex"], 5e-4)
        strong = evidence["findings"]["terminal_current_strong_inversion"]
        self.assertLess(strong["vela_drive_vela_hfs"]["maximum_dex"], 1e-3)
        for artifact in evidence["artifacts"].values():
            if isinstance(artifact, list):
                for item in artifact:
                    self.assertEqual(item["sha256"], sha256(ROOT / item["path"]))
            else:
                self.assertEqual(artifact["sha256"], sha256(ROOT / artifact["path"]))
        for relative, expected in evidence["implementation_sha256"].items():
            self.assertEqual(expected, sha256(ROOT / relative))


if __name__ == "__main__":
    unittest.main()
