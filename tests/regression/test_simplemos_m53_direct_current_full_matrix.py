import csv
import hashlib
import json
import math
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
OUTPUT = ROOT / "direct_current_full_matrix"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM53DirectCurrentFullMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract_path = (
            ROOT / "simplemos_m53_direct_current_full_matrix_contract_v1.json")
        cls.freeze = read_json(
            ROOT / "simplemos_m53_direct_current_full_matrix_contract_freeze.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m53_direct_current_full_matrix_evidence.json")
        cls.report = read_json(
            OUTPUT / "m53_direct_current_full_matrix_report.json")
        cls.cases = rows(OUTPUT / "m53_case_summary.csv")
        cls.points = rows(OUTPUT / "m53_pointwise_comparison.csv")
        cls.nwell = rows(OUTPUT / "m53_nwell_pair_summary.csv")
        cls.terminals = rows(OUTPUT / "m53_terminal_replay_ledger.csv")

    def test_contract_was_frozen_before_execution(self):
        self.assertEqual(self.freeze["status"], "frozen_before_execution")
        self.assertEqual(sha256(self.contract_path), self.freeze["contract_sha256"])
        self.assertEqual(self.evidence["contract_sha256"],
                         self.freeze["contract_sha256"])
        self.assertTrue(self.freeze["thresholds_frozen_before_execution"])

    def test_complete_exact_matrix(self):
        self.assertEqual(len(self.cases), 16)
        self.assertEqual(len(self.points), 816)
        self.assertEqual({int(row["point_count"]) for row in self.cases}, {51})
        self.assertEqual(
            {row["case"] for row in self.cases},
            {f"n{device}_vd_{drain}" for device in range(17, 25)
             for drain in ("0p05", "1")})
        for case in {row["case"] for row in self.points}:
            gates = sorted(float(row["gate_voltage_V"])
                           for row in self.points if row["case"] == case)
            self.assertEqual(len(gates), 51)
            for index, gate in enumerate(gates):
                self.assertAlmostEqual(gate, 0.05 * index, places=12)

    def test_direct_reference_curves_are_complete(self):
        references = sorted(
            (OUTPUT / "m53_direct_references").glob("*_reference.csv"))
        self.assertEqual(len(references), 16)
        self.assertTrue(all(len(rows(path)) == 51 for path in references))

    def test_nwell_and_terminal_ledgers_are_complete(self):
        self.assertEqual(len(self.nwell), 8)
        self.assertEqual(len(self.terminals), 72)
        self.assertEqual(
            {(row["low_nwell_device"], row["high_nwell_device"])
             for row in self.nwell},
            {("n17", "n21"), ("n18", "n22"),
             ("n19", "n23"), ("n20", "n24")})
        self.assertTrue(all(row["within_tolerance"] == "True"
                            for row in self.terminals))
        self.assertTrue(self.report["m52_terminal_replay"][
            "all_within_tolerance"])

    def test_classification_is_frozen_and_declared(self):
        declared = read_json(self.contract_path)["analysis"]["classifications"]
        self.assertIn(self.report["classification"], declared)
        self.assertEqual(self.evidence["classification"],
                         self.report["classification"])
        self.assertEqual(self.report["status"], "complete")
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])
        self.assertTrue(math.isfinite(float(self.report["global_findings"][
            "direct_maximum_absolute_log10_ratio_dex"])))

    def test_only_direct_current_was_intervened(self):
        self.assertEqual(self.evidence["only_intervention"],
                         "Math DirectCurrent flag")
        self.assertTrue(self.evidence["new_sentaurus_execution"])
        self.assertFalse(self.evidence["new_vela_execution"])
        self.assertFalse(self.evidence["production_defaults_changed"])
        self.assertFalse(self.evidence["historical_artifacts_rewritten"])
        self.assertFalse(self.evidence["closed_topics_reinvestigated"])
        manifest = read_json(
            REPO / "build-release/reference_tcad/simplemos_sentaurus2022/"
            "m53_direct_current_full_matrix/sentaurus_manifest.json")
        self.assertEqual(len(manifest["cases"]), 16)
        for item in manifest["cases"]:
            direct = (REPO / item["deck"]).read_text(encoding="utf-8")
            baseline = (REPO / item["baseline_deck"]).read_text(encoding="utf-8")
            self.assertEqual(direct.count(" DirectCurrent"), 1)
            self.assertNotIn("CurrentWeighting", direct)
            self.assertEqual(direct.replace(" DirectCurrent", ""), baseline)
            self.assertEqual(item["input_tdr_sha256"],
                             item["m8_input_tdr_sha256"])

    def test_frozen_artifact_and_source_hashes_match(self):
        self.assertEqual(self.evidence["status"], "frozen")
        for item in self.evidence["artifacts"]:
            self.assertEqual(sha256(REPO / item["path"]), item["sha256"])
        for relative, expected in self.evidence["source_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
