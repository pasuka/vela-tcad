import csv
import hashlib
import json
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
OUTPUT = ROOT / "doping_well_terminal_attribution"


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM55DopingWellTerminalAttributionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contract_path = (
            ROOT / "simplemos_m55_doping_well_terminal_attribution_contract_v1.json")
        cls.freeze = read_json(
            ROOT / "simplemos_m55_doping_well_terminal_attribution_contract_freeze.json")
        cls.evidence = read_json(
            ROOT / "simplemos_m55_doping_well_terminal_attribution_evidence.json")
        cls.report = read_json(
            OUTPUT / "m55_doping_well_terminal_attribution_report.json")
        cls.replay = rows(OUTPUT / "m55_default_terminal_replay_ledger.csv")
        cls.wells = rows(OUTPUT / "m55_doping_well_ledger.csv")
        cls.terminals = rows(OUTPUT / "m55_terminal_component_attribution_ledger.csv")
        cls.pairs = rows(OUTPUT / "m55_nwell_pair_summary.csv")

    def test_contract_was_frozen_before_execution(self):
        self.assertEqual(self.freeze["status"], "frozen_before_execution")
        self.assertEqual(sha256(self.contract_path), self.freeze["contract_sha256"])
        self.assertEqual(self.evidence["contract_sha256"],
                         self.freeze["contract_sha256"])

    def test_required_matrices_and_ledgers_are_complete(self):
        self.assertEqual(len(self.replay), 816)
        self.assertEqual(len(self.wells), 144)
        self.assertEqual(len(self.terminals), 576)
        self.assertEqual(len(self.pairs), 8)
        self.assertEqual({row["contact"] for row in self.wells},
                         {"source", "drain", "substrate"})
        self.assertEqual({float(row["gate_voltage_V"]) for row in self.wells},
                         {0.0, 0.05, 0.1})

    def test_default_terminal_and_shared_state_replays_are_exact(self):
        replay = self.report["default_terminal_replay"]
        self.assertEqual(replay["maximum_absolute_difference_A_per_um"], 0.0)
        self.assertEqual(replay["maximum_symmetric_relative_difference"], 0.0)
        self.assertTrue(replay["all_within_tolerance"])
        self.assertTrue(self.report["shared_state_replay"][
            "all_pointwise_identical"])
        self.assertEqual(self.report["shared_state_replay"]["mismatch_count"], 0)

    def test_nwell_changes_contact_associated_well_support(self):
        finding = self.report["doping_well"]
        self.assertTrue(finding[
            "topology_changed_across_at_least_one_nwell_pair"])
        self.assertAlmostEqual(
            finding["maximum_matched_pair_relative_area_difference"],
            0.1255721945407226, places=15)
        detail = finding["maximum_area_difference_pair"]
        self.assertEqual((detail["low_device"], detail["high_device"]),
                         ("n20", "n24"))
        self.assertEqual(detail["maximum_area_difference_contact"], "drain")
        self.assertLess(detail["maximum_area_difference_high_well_area_cm2"],
                        detail["maximum_area_difference_low_well_area_cm2"])

    def test_srh_charge_cancels_and_total_shift_is_surface_redistribution(self):
        attribution = self.report["terminal_attribution"]
        contract = read_json(self.contract_path)
        self.assertLessEqual(
            attribution["maximum_srh_electron_hole_cancellation_A_per_um"],
            contract["acceptance"][
                "srh_electron_hole_cancellation_absolute_tolerance_A_per_um"])
        self.assertLessEqual(
            attribution["maximum_total_surface_identity_residual_A_per_um"],
            contract["acceptance"][
                "total_surface_identity_absolute_tolerance_A_per_um"])
        target = {(row["contact"], row["component"]): row
                  for row in self.report["target"]["terminal_rows"]}
        self.assertAlmostEqual(
            target[("drain", "total")]["default_minus_direct_A_per_um"],
            -1.3364110459112917e-15, places=27)
        self.assertAlmostEqual(
            target[("substrate", "electron")][
                "absolute_generation_fraction_of_shift"],
            0.32597351370040956, places=15)

    def test_classification_is_declared_and_accepted(self):
        declared = read_json(self.contract_path)["analysis"]["classifications"]
        self.assertIn(self.report["classification"], declared)
        self.assertEqual(
            self.report["classification"],
            "nwell_changes_well_topology_surface_redistribution")
        self.assertEqual(self.evidence["classification"],
                         self.report["classification"])
        self.assertEqual(self.report["status"], "accepted")
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])

    def test_execution_scope_and_frozen_hashes(self):
        self.assertTrue(self.evidence["new_sentaurus_execution"])
        self.assertFalse(self.evidence["new_vela_execution"])
        self.assertFalse(self.evidence["historical_artifacts_rewritten"])
        self.assertFalse(self.evidence["closed_topics_reinvestigated"])
        for relative, expected in self.evidence["artifacts"].items():
            self.assertEqual(sha256(REPO / relative), expected)
        for relative, expected in self.evidence["source_hashes"].items():
            self.assertEqual(sha256(REPO / relative), expected)


if __name__ == "__main__":
    unittest.main()
