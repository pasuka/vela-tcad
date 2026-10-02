from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "simplemos_m8b_hfs_controls",
    ROOT / "scripts" / "run_simplemos_m8b_hfs_controls.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
CONTRACT_PATH = (
    ROOT / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_m8b_hfs_controls_contract_v1.json"
)
EVIDENCE_PATH = (
    ROOT / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_m8b_hfs_controls_evidence.json"
)


class SimpleMosM8BHfsControlsTest(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = MODULE.read_json(CONTRACT_PATH)
        MODULE.validate_contract(self.contract)
        self.variants = {item["id"]: item
                         for item in self.contract["variants"]}

    def test_contract_freezes_single_variable_matrix(self) -> None:
        self.assertEqual(["n17", "n21"], [
            item["id"] for item in self.contract["devices"]])
        self.assertEqual(5, len(self.variants))
        self.assertEqual([0.05, 1.0], self.contract[
            "bias_matrix"]["drain_voltages_V"])
        self.assertEqual("forbidden", self.contract[
            "comparison"]["interpolation"])

    def test_explicit_gradqf_is_a_null_control(self) -> None:
        deck = MODULE.sentaurus_deck(
            "fixture", 0.05, self.variants["explicit_gradqf"], 50)
        self.assertIn(
            "Mobility(PhuMob HighFieldSaturation(GradQuasiFermi) Enormal)",
            deck)
        self.assertNotIn("ComputeGradQuasiFermiAtContacts", deck)
        self.assertNotIn("RefDens_", deck)

    def test_each_control_changes_only_declared_hfs_setting(self) -> None:
        eparallel = MODULE.sentaurus_deck(
            "fixture", 0.05, self.variants["eparallel"], 50)
        self.assertIn("HighFieldSaturation(Eparallel)", eparallel)
        self.assertIn("Recombination(SRH(DopingDependence))", eparallel)

        contacts = MODULE.sentaurus_deck(
            "fixture", 0.05, self.variants["qf_at_contacts"], 50)
        self.assertIn(
            "ComputeGradQuasiFermiAtContacts=UseQuasiFermi", contacts)

        boundary = MODULE.sentaurus_deck(
            "fixture", 0.05, self.variants["no_parallel_boundary"], 50)
        self.assertIn("-ParallelToInterfaceInBoundaryLayer", boundary)

        refdens = MODULE.sentaurus_deck(
            "fixture", 0.05, self.variants["refdens_efield_1e8"], 50)
        self.assertIn(
            "RefDens_eGradQuasiFermi_ElectricField_HFS=1e8", refdens)
        self.assertIn(
            "RefDens_hGradQuasiFermi_ElectricField_HFS=1e8", refdens)

    def test_parallel_executor_rejects_zero_workers(self) -> None:
        with self.assertRaisesRegex(ValueError, "at least one"):
            MODULE.run_sentaurus(
                {}, Path("unused"), "host", "ssh", "scp", "remote", 0,
                Path("unused_raw"))

    def test_prepared_cases_use_isolated_parameter_directories(self) -> None:
        deck = MODULE.sentaurus_deck(
            "fixture", 0.05, self.variants["explicit_gradqf"], 50).replace(
                'Grid="input_fps.tdr"', 'Grid="../../input_fps.tdr"')
        self.assertIn('Grid="../../input_fps.tdr"', deck)

    def test_declared_retry_changes_numerics_not_physics_or_bias_lattice(self) -> None:
        retries = self.contract["numerical_retry_policy"]["cases"]
        self.assertEqual(
            ["n17_eparallel_vd_1", "n21_eparallel_vd_1"],
            [item["case"] for item in retries])
        for retry in retries:
            deck = MODULE.sentaurus_deck(
                retry["case"], 1.0, self.variants["eparallel"], 50, retry)
            self.assertIn("HighFieldSaturation(Eparallel)", deck)
            self.assertIn("Iterations=100", deck)
            self.assertIn("MinStep=1e-08", deck)
            self.assertIn("LineSearchDamping=0.01", deck)
            self.assertEqual(1, deck.count("LineSearchDamping=0.01"))
            self.assertIn(
                "CurrentPlot(Time=(Range=(0 1) Intervals=50))", deck)

    def test_frozen_evidence_records_complete_control_diagnosis(self) -> None:
        evidence = MODULE.read_json(EVIDENCE_PATH)
        self.assertEqual("complete", evidence["status"])
        self.assertEqual(20, evidence["execution"]["curves"])
        self.assertEqual(1020, evidence["execution"][
            "total_direct_bias_points"])
        self.assertEqual("forbidden", evidence["execution"]["interpolation"])
        self.assertTrue(evidence["findings"]["null_control_pass"])

        report_artifact = evidence["artifacts"]["comparison_report"]
        report_path = ROOT / report_artifact["path"]
        self.assertEqual(report_artifact["sha256"],
                         MODULE.sha256(report_path))
        report = MODULE.read_json(report_path)
        self.assertEqual(20, len(report["cases"]))
        for case in report["cases"]:
            for key in ("comparison_csv", "response_csv"):
                path = ROOT / case[key]
                self.assertEqual(case[f"{key}_sha256"],
                                 MODULE.sha256(path))
        self.assertEqual(4, len(evidence["artifacts"]["figures"]))
        for figure in evidence["artifacts"]["figures"]:
            path = ROOT / figure["path"]
            self.assertEqual(figure["sha256"], MODULE.sha256(path))
        for relative, expected in evidence["implementation_sha256"].items():
            self.assertEqual(expected, MODULE.sha256(ROOT / relative))


if __name__ == "__main__":
    unittest.main()
