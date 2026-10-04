import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def sha256(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM61PathSensitivityGateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.evidence = read_json(ROOT / "simplemos_m61_path_sensitivity_gate_evidence.json")
        cls.report = read_json(ROOT / "path_sensitivity_gate/m61_path_sensitivity_gate_report.json")

    def test_decisive_m60_skips_m61_execution(self):
        self.assertEqual(self.report["decision"], "not_required_m60_decisive")
        self.assertTrue(self.report["acceptance"]["m60_classification_decisive"])
        self.assertTrue(self.report["acceptance"]["m60_zero_remaining_bursts"])
        self.assertEqual(self.report["execution"]["deck_count"], 0)
        self.assertFalse(self.report["execution"]["new_sentaurus_execution"])

    def test_sequence_closure_is_complete(self):
        sequence = self.report["sequence_closure"]
        self.assertEqual(sequence["m63"], "threshold_like_horizontal_shift_dominant")
        self.assertEqual(sequence["m62"], "e3_stopped_no_identity_preserving_writer")
        self.assertEqual(sequence["m61"], "not_required_m60_decisive")

    def test_evidence_is_frozen(self):
        self.assertEqual(self.report["status"], "accepted")
        self.assertEqual(self.evidence["status"], "frozen")
        for relative, expected in self.evidence["source_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)
        for relative, expected in self.evidence["artifacts"].items():
            self.assertEqual(sha256(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
