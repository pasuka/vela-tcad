from __future__ import annotations

import importlib.util
import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "simplemos_m4", ROOT / "scripts" / "run_simplemos_m4_controlled_matrix.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
COMPARE_SPEC = importlib.util.spec_from_file_location(
    "simplemos_m4_compare",
    ROOT / "scripts" / "compare_simplemos_m4_controlled_matrix.py")
assert COMPARE_SPEC and COMPARE_SPEC.loader
COMPARE = importlib.util.module_from_spec(COMPARE_SPEC)
COMPARE_SPEC.loader.exec_module(COMPARE)
CONTRACT_PATH = (
    ROOT / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_m4_controlled_mobility_contract_v1.json"
)
EVIDENCE_PATH = (
    ROOT / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_m4_controlled_mobility_evidence.json"
)


class SimpleMosM4ControlledMatrixTest(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = MODULE.read_json(CONTRACT_PATH)

    def test_contract_freezes_cumulative_a0_a3_ladder_and_exact_lattice(self) -> None:
        MODULE.validate_contract(self.contract)
        self.assertEqual(
            ["constant", "masetti", "masetti_field", "masetti_field_lombardi"],
            [item["vela_mobility_model"] for item in self.contract["variants"]])
        lattice = self.contract["bias_matrix"]["gate_lattice"]["values_V"]
        self.assertEqual(51, len(lattice))
        self.assertEqual([0.0, 0.05, 0.1], lattice[:3])
        self.assertEqual(2.5, lattice[-1])
        self.assertEqual("forbidden", self.contract["comparison"]["interpolation"])

    def test_checked_in_evidence_matches_contract_and_reference_manifest(self) -> None:
        evidence = MODULE.read_json(EVIDENCE_PATH)
        reference_manifest_path = (
            EVIDENCE_PATH.parent / evidence["reference_manifest"])
        reference_manifest = MODULE.read_json(reference_manifest_path)
        self.assertEqual("accepted", evidence["status"])
        self.assertEqual(
            hashlib.sha256(CONTRACT_PATH.read_bytes()).hexdigest(),
            evidence["contract_sha256"])
        self.assertEqual(
            hashlib.sha256(reference_manifest_path.read_bytes()).hexdigest(),
            evidence["reference_manifest_sha256"])
        self.assertEqual(8, len(reference_manifest["artifacts"]))
        self.assertTrue(all(item["point_count"] == 51
                            for item in reference_manifest["artifacts"]))
        lattice = self.contract["bias_matrix"]["gate_lattice"]["values_V"]
        tolerance = self.contract["comparison"]["exact_bias_tolerance_V"]
        for artifact in reference_manifest["artifacts"]:
            curve = reference_manifest_path.parent / artifact["path"]
            self.assertEqual(artifact["sha256"], MODULE.sha256(curve))
            points = COMPARE.exact_curve(
                curve, lattice, "gate_voltage_V",
                "drain_total_current_A_per_um", tolerance, False)
            self.assertEqual(51, len(points))
        cases = evidence["comparison"]["cases"]
        self.assertEqual(8, len(cases))
        self.assertTrue(all(item["status"] == "pass" for item in cases))
        self.assertEqual(
            max(item["max_log10_ratio"] for item in cases),
            evidence["comparison"]["maximum_absolute_log10_ratio_observed"])
        comparison_report = EVIDENCE_PATH.parent / evidence["comparison"]["report"]
        self.assertEqual(
            evidence["comparison"]["report_sha256"],
            hashlib.sha256(comparison_report.read_bytes()).hexdigest())
        figure_dir = ROOT / "docs" / "validation" / "figures" / "simplemos_m4"
        for name, expected_hash in evidence["visualization"]["figures"].items():
            self.assertEqual(
                expected_hash,
                hashlib.sha256((figure_dir / name).read_bytes()).hexdigest())
        for relative, expected_hash in evidence["implementation_sha256"].items():
            self.assertEqual(
                expected_hash,
                hashlib.sha256((ROOT / relative).read_bytes()).hexdigest())

    def test_sentaurus_decks_change_only_declared_mobility_ladder(self) -> None:
        decks = [MODULE.sentaurus_deck(
            item["id"].lower(), 0.05,
            item["sentaurus_mobility_models"], 50)
            for item in self.contract["variants"]]
        self.assertNotIn("Mobility(", decks[0])
        self.assertIn("Mobility( DopingDependence )", decks[1])
        self.assertIn(
            "Mobility( DopingDependence HighFieldSaturation )", decks[2])
        self.assertIn(
            "Mobility( DopingDependence HighFieldSaturation Enormal )", decks[3])
        for deck in decks:
            self.assertIn("OldSlotboom NoFermi", deck)
            self.assertIn("SRH(DopingDependence)", deck)
            self.assertIn("CurrentPlot(Time=(Range=(0 1) Intervals=50))", deck)

    def test_a3_vd1_recovery_changes_only_auxiliary_drain_step(self) -> None:
        default = MODULE.auxiliary_drain_schedule(self.contract, "A2", 1.0)
        recovery = MODULE.auxiliary_drain_schedule(self.contract, "A3", 1.0)
        self.assertEqual(0.001, default["initial_step_V"])
        self.assertEqual(0.005, recovery["initial_step_V"])
        self.assertEqual("none", recovery["acceptance_role"])
        lattice = self.contract["bias_matrix"]["gate_lattice"]
        self.assertEqual(51, lattice["point_count"])
        self.assertEqual(0.05, lattice["step_V"])

    def test_exact_plt_curve_rejects_missing_lattice_point(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_simplemos_m4_") as temporary:
            path = Path(temporary) / "curve.plt"
            path.write_text(
                'Info { datasets = ["gate OuterVoltage" "drain TotalCurrent"] }\n'
                "Data { 0 1e-12 0.1 2e-12 }\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "exactly one row"):
                MODULE.exact_plt_curve(path, [0.0, 0.05, 0.1], 1e-10)

    def test_neutral_mesh_and_vela_workflows_preserve_exact_biases(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_simplemos_m4_") as temporary:
            root = Path(temporary)
            neutral = root / "neutral"
            neutral.mkdir()
            (neutral / "nodes.csv").write_text(
                "id,x_um,y_um\n0,0,0\n1,1,0\n2,0,1\n", encoding="utf-8")
            (neutral / "elements.csv").write_text(
                "id,node0,node1,node2,region,material\n"
                "0,0,1,2,Silicon_1,Si\n", encoding="utf-8")
            (neutral / "contacts.csv").write_text(
                "name,node_ids,region\n"
                "source,0;1,Silicon_1\n"
                "drain,1;2,Silicon_1\n"
                "gate,0;2,Silicon_1\n"
                "substrate,0;1;2,Silicon_1\n", encoding="utf-8")
            (neutral / "doping.csv").write_text(
                "node_id,donors_cm3,acceptors_cm3\n"
                "0,0,1e17\n1,1e20,0\n2,0,1e17\n", encoding="utf-8")
            materials = root / "materials.json"
            materials.write_text(json.dumps({"materials": [{"name": "Si"}]}))

            manifest = MODULE.prepare_vela_workflows(
                self.contract, neutral, materials, root / "run")
            self.assertEqual(24, len(manifest["stages"]))
            mesh = MODULE.read_json(Path(manifest["mesh"]))
            self.assertEqual(3, len(mesh["nodes"]))
            self.assertEqual(1, len(mesh["triangles"]))
            gate_stages = [item for item in manifest["stages"]
                           if item["phase"] == "gate_sweep"]
            self.assertEqual(8, len(gate_stages))
            for stage in gate_stages:
                config = MODULE.read_json(Path(stage["config"]))
                self.assertEqual(
                    self.contract["bias_matrix"]["gate_lattice"]["values_V"],
                    config["sweep"]["bias_points"])
                self.assertEqual(
                    "sentaurus_default",
                    config["solver"]["srh_density_coupling"])
                self.assertNotIn("line_search_mode", config["solver"])
            drain_stages = [item for item in manifest["stages"]
                            if item["phase"] == "drain_ramp"]
            drain_steps = []
            for stage in drain_stages:
                config = MODULE.read_json(Path(stage["config"]))
                drain_steps.append(abs(config["sweep"]["step"]))
                self.assertLessEqual(config["sweep"]["max_step"], 0.01)
                self.assertEqual("block_filter", config["solver"]["line_search_mode"])
                self.assertEqual(
                    "sentaurus_default",
                    config["solver"]["srh_density_coupling"])
                self.assertTrue(config["simplemos_m4"][
                    "does_not_change_acceptance_lattice"])
            self.assertEqual([0.001] * 7 + [0.005], sorted(drain_steps))

    def write_curve(self, path: Path, scale: float, *, candidate: bool,
                    converged: str = "1") -> None:
        lattice = self.contract["bias_matrix"]["gate_lattice"]["values_V"]
        if candidate:
            header = "bias_V,current_total_A_per_um,converged\n"
            rows = "".join(
                f"{bias},{scale * (1e-16 + (bias + 0.1) ** 4 * 1e-6)},{converged}\n"
                for bias in lattice)
        else:
            header = "gate_voltage_V,drain_total_current_A_per_um\n"
            rows = "".join(
                f"{bias},{1e-16 + (bias + 0.1) ** 4 * 1e-6}\n"
                for bias in lattice)
        path.write_text(header + rows, encoding="utf-8")

    def test_exact_comparison_passes_without_interpolation(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_simplemos_m4_") as temporary:
            root = Path(temporary)
            reference = root / "reference.csv"
            candidate = root / "candidate.csv"
            self.write_curve(reference, 1.0, candidate=False)
            self.write_curve(candidate, 1.1, candidate=True)
            report = COMPARE.compare_case(reference, candidate, self.contract)
            self.assertEqual("pass", report["status"])
            self.assertEqual(51, report["point_count"])
            self.assertAlmostEqual(
                math.log10(1.1),
                report["maximum_absolute_log10_ratio_above_floor"])

    def test_comparison_rejects_nonconverged_candidate(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_simplemos_m4_") as temporary:
            root = Path(temporary)
            reference = root / "reference.csv"
            candidate = root / "candidate.csv"
            self.write_curve(reference, 1.0, candidate=False)
            self.write_curve(candidate, 1.0, candidate=True, converged="0")
            with self.assertRaisesRegex(ValueError, "non-converged"):
                COMPARE.compare_case(reference, candidate, self.contract)


if __name__ == "__main__":
    unittest.main()
