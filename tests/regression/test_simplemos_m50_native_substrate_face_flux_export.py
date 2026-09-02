import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
OUTPUT = ROOT / "native_substrate_face_flux_export"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM50NativeSubstrateFaceFluxExportTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract_path = ROOT / "simplemos_m50_native_substrate_face_flux_export_contract_v1.json"
        cls.freeze = read_json(
            ROOT / "simplemos_m50_native_substrate_face_flux_export_contract_freeze.json")
        cls.erratum_path = ROOT / "simplemos_m50_native_substrate_face_flux_export_contract_execution_erratum_v1.json"
        cls.erratum_freeze = read_json(
            ROOT / "simplemos_m50_native_substrate_face_flux_export_contract_execution_erratum_freeze.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m50_native_substrate_face_flux_export_evidence.json")
        cls.report = read_json(
            OUTPUT / "m50_native_substrate_face_flux_export_report.json")
        cls.flux = rows(OUTPUT / "m50_native_substrate_face_flux_ledger.csv")
        cls.replay = rows(OUTPUT / "m50_terminal_replay_ledger.csv")

    def test_contract_and_execution_erratum_were_frozen(self):
        self.assertEqual(self.freeze["status"], "frozen_before_execution")
        self.assertEqual(sha256(self.contract_path), self.freeze["contract_sha256"])
        self.assertEqual(self.erratum_freeze["status"],
                         "frozen_before_retry_execution")
        self.assertEqual(sha256(self.erratum_path),
                         self.erratum_freeze["erratum_sha256"])
        self.assertEqual(self.evidence["contract_sha256"],
                         self.freeze["contract_sha256"])

    def test_six_exact_states_and_terminal_replay_are_complete(self):
        self.assertEqual(len(self.flux), 6)
        self.assertEqual(len(self.replay), 72)
        self.assertEqual(
            {(row["device"], float(row["gate_voltage_V"])) for row in self.flux},
            {(device, gate) for device in ("n23", "n19")
             for gate in (0.0, 0.05, 0.1)})
        self.assertTrue(self.report["acceptance"]["terminal_replay_closes"])

    def test_native_orientation_and_bridge_are_explicit(self):
        allowed = {
            "contact_surface_normal_integral_matches_terminal",
            "negative_contact_surface_normal_integral_matches_terminal",
            "no_terminal_closing_orientation",
        }
        self.assertTrue({row["orientation"] for row in self.flux} <= allowed)
        self.assertIn(self.report["classification"], {
            "native_bridge_closed", "native_observable_mismatch",
            "replay_state_mismatch"})
        self.assertTrue(self.report["acceptance"][
            "native_matches_m49_exported_boundary"])
        self.assertLessEqual(
            self.report["native_observable"][
                "maximum_negative_native_to_m49_boundary_absolute_error_A_per_um"],
            1e-30)
        if self.report["classification"] == "native_bridge_closed":
            self.assertTrue(self.report["acceptance"][
                "native_to_same_run_terminal_closes"])
            self.assertTrue(self.report["acceptance"][
                "native_bridge_identity_closes"])

    def test_target_anchor_and_controls_are_reported(self):
        target = self.report["target"]
        self.assertEqual(target["device"], "n23")
        self.assertEqual(float(target["gate_voltage_V"]), 0.05)
        self.assertAlmostEqual(
            float(target["frozen_m48_vela_minus_sentaurus_A_per_um"]),
            -4.6947034945330634e-17, delta=1e-30)
        self.assertEqual(len(self.report["control_localization"]["controls"]), 3)

    def test_execution_did_not_reopen_closed_topics(self):
        self.assertTrue(self.evidence["new_sentaurus_execution"])
        self.assertFalse(self.evidence["new_vela_execution"])
        self.assertFalse(self.evidence["default_physics_model_changed"])
        self.assertFalse(self.evidence["closed_topics_reinvestigated"])

    def test_frozen_artifact_and_source_hashes_match(self):
        self.assertEqual(self.evidence["status"], "frozen")
        for item in self.evidence["artifacts"]:
            self.assertEqual(sha256(REPO / item["path"]), item["sha256"])
        for relative, expected in self.evidence["source_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
