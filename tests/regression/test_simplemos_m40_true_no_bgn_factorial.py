import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/true_no_bgn_factorial"
REPORT = ROOT / "m40_true_no_bgn_factorial_report.json"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m40_true_no_bgn_factorial_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM40TrueNoBgnFactorialTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads(REPORT.read_text(encoding="utf-8"))
        cls.evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))

    def test_manual_and_all_execution_gates_pass(self) -> None:
        self.assertEqual(self.report["status"], "complete")
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])
        self.assertTrue(all(self.report["acceptance"].values()))
        manual = self.report["manual_confirmation"]
        self.assertEqual(manual["release"], "T-2022.03-SP2")
        self.assertEqual(manual["confirmed_syntax"],
                         "EffectiveIntrinsicDensity(NoBandGapNarrowing)")

    def test_no_bgn_is_confirmed_by_log_and_field(self) -> None:
        rows = {row["cell"]: row
                for row in self.report["sentaurus_model_observability"]}
        for cell in ("bgn_off_srh_on", "bgn_off_srh_off"):
            row = rows[cell]
            self.assertEqual(row["expected_model"], "none")
            self.assertIn("without bandgap narrowing", row["solver_log_line"])
            self.assertEqual(row["bandgap_narrowing_max_abs_eV"], 0.0)

    def test_donor_and_acceptor_observability_is_exact(self) -> None:
        for row in self.report["sentaurus_model_observability"]:
            self.assertEqual(row["doping_common_node_count"], 942)
            self.assertEqual(row["doping_max_relative_difference"], 0.0)

    def test_srh_on_conditional_bgn_result(self) -> None:
        result = self.report["causal_result"]
        self.assertAlmostEqual(result["production_gap_dex"],
                               0.09039871358507992, places=12)
        self.assertAlmostEqual(result["true_no_bgn_srh_on_gap_dex"],
                               0.010870601028172609, places=12)
        self.assertGreater(result["bgn_related_fraction_of_production_gap"], 0.87)
        self.assertLess(result["bgn_related_fraction_of_production_gap"], 0.89)
        self.assertFalse(result["factorial_main_effect_qualified"])

    def test_frozen_artifact_and_source_hashes_match(self) -> None:
        self.assertEqual(self.evidence["status"], "frozen")
        for artifact in self.evidence["artifacts"]:
            self.assertEqual(artifact["sha256"], sha256(REPO / artifact["path"]))
        for relative, expected in self.evidence["source_hashes"].items():
            self.assertEqual(expected, sha256(REPO / relative))


if __name__ == "__main__":
    unittest.main()
