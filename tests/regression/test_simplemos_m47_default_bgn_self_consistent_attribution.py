import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM47EvidenceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract_path = ROOT / "simplemos_m47_default_bgn_self_consistent_attribution_contract_v1.json"
        cls.freeze = read_json(
            ROOT / "simplemos_m47_default_bgn_self_consistent_attribution_contract_freeze.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m47_default_bgn_self_consistent_attribution_evidence.json")
        cls.report = read_json(
            ROOT / "default_bgn_state_attribution/m47_default_bgn_self_consistent_attribution_report.json")
        cls.summary = rows(ROOT / "default_bgn_state_attribution/m47_state_summary.csv")
        cls.flux = rows(ROOT / "default_bgn_state_attribution/m47_contact_flux_ledger.csv")
        cls.contrasts = rows(ROOT / "default_bgn_state_attribution/m47_attribution_contrasts.csv")

    def test_contract_was_frozen_before_execution(self):
        self.assertEqual(self.freeze["status"], "frozen_before_execution")
        self.assertEqual(sha256(self.contract_path), self.freeze["contract_sha256"])
        self.assertEqual(self.evidence["contract_sha256"], self.freeze["contract_sha256"])

    def test_exact_six_state_control_matrix(self):
        self.assertEqual(len(self.summary), 6)
        observed = {(row["device"], float(row["gate_voltage_V"])) for row in self.summary}
        expected = {(device, gate) for device in ("n23", "n19")
                    for gate in (0.0, 0.05, 0.1)}
        self.assertEqual(observed, expected)
        self.assertEqual(len(self.flux), 54)
        balances = [row for row in self.flux
                    if row["contact"] == "source_plus_drain"]
        self.assertEqual(len(balances), 18)
        self.assertTrue(self.report["acceptance"]["contact_flux_complete"])

    def test_m46_anchor_and_qf_coordinate_are_preserved(self):
        self.assertTrue(self.report["acceptance"]["m46_point_currents_reproduced"])
        self.assertTrue(self.report["acceptance"]["qf_reference_increment_preserved"])
        self.assertLessEqual(max(float(row["vela_anchor_shift_dex"])
                                 for row in self.summary), 1e-10)
        self.assertLessEqual(max(float(row["sentaurus_anchor_shift_dex"])
                                 for row in self.summary), 1e-10)

    def test_required_attribution_observables_exist(self):
        self.assertEqual(len(self.contrasts), 7)
        required = {
            "barrier_predicted_current_shift_dex",
            "electrostatic_barrier_delta_mV",
            "source_phin_signed_mean_mV",
            "source_electron_density_signed_mean_dex",
            "source_hole_density_signed_mean_dex",
            "source_barrier_electron_density_log_ratio_dex",
            "srh_integrated_log_magnitude_ratio_dex",
        }
        self.assertEqual({row["metric"] for row in self.contrasts}, required)

    def test_closed_topics_and_defaults_stayed_closed(self):
        self.assertTrue(self.report["acceptance"]["closed_topics_not_reopened"])
        self.assertFalse(self.evidence["default_physics_model_changed"])
        self.assertFalse(self.evidence["default_hfs_model_changed"])
        self.assertFalse(self.evidence["sg_or_contact_extraction_reinvestigated"])
        self.assertFalse(self.evidence["qf_packing_reinvestigated"])
        self.assertTrue(self.evidence["acceptance"]["all_checks_pass"])

    def test_frozen_artifact_hashes(self):
        for item in self.evidence["artifacts"]:
            self.assertEqual(sha256(REPO / item["path"]), item["sha256"])


if __name__ == "__main__":
    unittest.main()
