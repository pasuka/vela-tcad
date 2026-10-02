from __future__ import annotations

import importlib.util
import hashlib
import json
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "simplemos_m8", ROOT / "scripts" / "run_simplemos_m8_original_matrix.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
CONTRACT_PATH = (
    ROOT / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_m8_original_physics_contract_v1.json"
)
EVIDENCE_PATH = (
    ROOT / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_m8_original_physics_evidence.json"
)


class SimpleMosM8OriginalMatrixTest(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = MODULE.read_json(CONTRACT_PATH)

    def test_contract_freezes_original_physics_and_eight_device_matrix(self) -> None:
        MODULE.validate_contract(self.contract)
        physics = self.contract["physics"]
        self.assertEqual("OldSlotboom", physics["effective_intrinsic_density"])
        self.assertEqual(
            ["PhuMob", "HighFieldSaturation", "Enormal"],
            physics["mobility"])
        self.assertEqual("SRH(DopingDependence)", physics["recombination"])
        srh = physics["srh_scharfetter_defaults"]
        self.assertEqual(1.0e-5, srh["electron"]["tau_max_s"])
        self.assertEqual(3.0e-6, srh["hole"]["tau_max_s"])
        self.assertEqual(1.0e16, srh["electron"]["reference_doping_cm3"])
        self.assertEqual(list(range(17, 25)), [
            item["workbench_process_node"] for item in self.contract["devices"]
        ])
        combinations = {
            (item["NWell_cm3"], item["GOxTime_min"], item["LDD_Dose_cm2"])
            for item in self.contract["devices"]
        }
        self.assertEqual(8, len(combinations))
        self.assertEqual("n17", self.contract["execution_gates"]["nominal_device"])

    def test_contract_requires_exact_gate_lattice_without_interpolation(self) -> None:
        lattice = self.contract["bias_matrix"]["gate_lattice"]["values_V"]
        self.assertEqual([round(index * 0.05, 12) for index in range(51)], lattice)
        self.assertEqual([0.05, 1.0], self.contract["bias_matrix"]["drain_voltages_V"])
        self.assertEqual("forbidden", self.contract["comparison"]["interpolation"])
        self.assertTrue(self.contract["execution_gates"][
            "nominal_two_curve_gate_must_pass_before_matrix"])

    def test_sentaurus_deck_uses_original_physics_and_exact_51_points(self) -> None:
        deck = MODULE.sentaurus_deck("n17_vd_0p05", 0.05, 50)
        self.assertIn("EffectiveIntrinsicDensity(OldSlotboom)", deck)
        self.assertIn("Mobility(PhuMob HighFieldSaturation Enormal)", deck)
        self.assertIn("Recombination(SRH(DopingDependence))", deck)
        self.assertIn("CurrentPlot(Time=(Range=(0 1) Intervals=50))", deck)
        self.assertNotIn("DopingDependence HighFieldSaturation Enormal", deck)

    def test_sentaurus_bundle_keeps_each_tdr_distinct(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_simplemos_m8_") as temporary:
            root = Path(temporary)
            tdrs = {}
            for device in self.contract["devices"][:2]:
                path = root / f"{device['id']}.tdr"
                path.write_bytes(f"fixture-{device['id']}".encode("ascii"))
                tdrs[device["id"]] = path
            devices = self.contract["devices"][:2]
            manifest = MODULE.prepare_sentaurus_bundle(
                self.contract, devices, tdrs, root / "output")
            self.assertEqual(2, len(manifest["inputs"]))
            self.assertEqual(4, len(manifest["cases"]))
            self.assertEqual(
                {"n17_vd_0p05", "n17_vd_1", "n18_vd_0p05", "n18_vd_1"},
                {item["case"] for item in manifest["cases"]})
            self.assertNotEqual(
                manifest["inputs"][0]["tdr_sha256"],
                manifest["inputs"][1]["tdr_sha256"])

    def test_vela_base_uses_simplemos_default_srh_lifetimes(self) -> None:
        config = MODULE.m4.base_config(
            Path("mesh.json"), Path("doping.csv"), Path("materials.json"),
            {"model": "phumob_field_lombardi"})
        MODULE.apply_simplemos_srh_defaults(config, self.contract)
        srh_config = config["solver"]["srh_doping_dependence"]
        self.assertEqual(1.0e-5, srh_config["electron"]["tau_max_s"])
        self.assertEqual(3.0e-6, srh_config["hole"]["tau_max_s"])

    def test_parallel_execution_rejects_zero_workers(self) -> None:
        with self.assertRaisesRegex(ValueError, "at least one"):
            MODULE.execute_vela([], Path("unused"), Path("runner"), jobs=0)

    def test_tdr_spec_and_device_selection_reject_ambiguous_inputs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_simplemos_m8_") as temporary:
            path = Path(temporary) / "input.tdr"
            path.write_bytes(b"fixture")
            with self.assertRaisesRegex(ValueError, "duplicate TDR"):
                MODULE.parse_tdr_specs([f"n17={path}", f"n17={path}"])
            with self.assertRaisesRegex(ValueError, "unknown M8 device"):
                MODULE.selected_devices(self.contract, ["n99"])

    def test_non_nominal_devices_require_current_two_curve_gate(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_simplemos_m8_") as temporary:
            root = Path(temporary)
            devices = MODULE.selected_devices(self.contract, ["n18"])
            with self.assertRaisesRegex(RuntimeError, "run and compare n17 first"):
                MODULE.require_nominal_gate(CONTRACT_PATH, devices, root)
            gate_dir = root / "comparisons"
            gate_dir.mkdir()
            (gate_dir / MODULE.NOMINAL_GATE_FILE).write_text(json.dumps({
                "status": "pass",
                "contract_sha256": MODULE.sha256(CONTRACT_PATH),
                "cases": [{"status": "pass"}, {"status": "pass"}],
            }), encoding="utf-8")
            MODULE.require_nominal_gate(CONTRACT_PATH, devices, root)

    def test_checked_in_evidence_records_complete_accepted_matrix(self) -> None:
        evidence = MODULE.read_json(EVIDENCE_PATH)
        self.assertEqual("accepted", evidence["status"])
        self.assertEqual(
            hashlib.sha256(CONTRACT_PATH.read_bytes()).hexdigest(),
            evidence["contract_sha256"])
        self.assertEqual(8, len(evidence["input_tdrs"]))
        reference_manifest_path = EVIDENCE_PATH.parent / evidence[
            "reference_manifest"]
        reference_manifest = MODULE.read_json(reference_manifest_path)
        self.assertEqual(
            evidence["reference_manifest_sha256"],
            hashlib.sha256(reference_manifest_path.read_bytes()).hexdigest())
        self.assertEqual(16, len(reference_manifest["artifacts"]))
        self.assertTrue(all(item["point_count"] == 51
                            for item in reference_manifest["artifacts"]))
        for artifact in reference_manifest["artifacts"]:
            curve = reference_manifest_path.parent / artifact["path"]
            self.assertEqual(artifact["sha256"], MODULE.sha256(curve))

        comparison_path = EVIDENCE_PATH.parent / evidence["comparison"]["report"]
        comparison = MODULE.read_json(comparison_path)
        self.assertEqual("pass", comparison["status"])
        self.assertEqual(
            evidence["comparison"]["report_sha256"],
            hashlib.sha256(comparison_path.read_bytes()).hexdigest())
        self.assertEqual(16, len(comparison["cases"]))
        passed = {item["device"] for item in comparison["cases"]
                  if item["status"] == "pass"}
        failed = {item["device"] for item in comparison["cases"]
                  if item["status"] == "fail"}
        self.assertEqual({f"n{index}" for index in range(17, 25)}, passed)
        self.assertEqual(set(), failed)
        self.assertTrue(all(item["trend_match"] for item in comparison["cases"]))
        repair_path = EVIDENCE_PATH.parent / evidence["repair"]["evidence"]
        self.assertEqual(
            evidence["repair"]["evidence_sha256"], MODULE.sha256(repair_path))
        for item in comparison["cases"]:
            comparison_csv = Path(item["comparison_csv"])
            self.assertFalse(comparison_csv.is_absolute())
            if not comparison_csv.is_absolute():
                comparison_csv = ROOT / comparison_csv
            self.assertEqual(item["comparison_csv_sha256"],
                             MODULE.sha256(comparison_csv))
        self.assertEqual("generated", evidence["visualization"]["status"])
        self.assertEqual(4, len(evidence["visualization"]["figures"]))
        for figure in evidence["visualization"]["figures"]:
            figure_path = ROOT / figure["path"]
            self.assertTrue(figure_path.is_file())
            self.assertEqual(figure["sha256"], MODULE.sha256(figure_path))
        for relative, expected_hash in evidence["implementation_sha256"].items():
            self.assertEqual(expected_hash, MODULE.sha256(ROOT / relative))


if __name__ == "__main__":
    unittest.main()
