#!/usr/bin/env python3
"""Regression coverage for Templates/LDMOS phase-2/3 evidence tooling."""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from prepare_templates_ldmos_phase23 import (  # noqa: E402
    SENTAURUS_G3_DRAIN_PREBIAS_POINTS_V,
    classical_solver,
    derive_polysi_flatband,
    exact_bias_points,
    prepare as prepare_phase23,
)
from prepare_templates_ldmos_wp15_diagnostics import (  # noqa: E402
    prepare as prepare_wp15_diagnostics,
)
from prepare_templates_ldmos_sentaurus_ablation import (  # noqa: E402
    prepare as prepare_ablation,
)
from summarize_templates_ldmos_phase23 import (  # noqa: E402
    field_metrics,
    percentile,
)
from summarize_templates_ldmos_sentaurus_ablation import (  # noqa: E402
    compare as compare_ablation,
    crossing,
    wall_seconds,
)


def physics_contract() -> dict:
    return {
        "carrier_statistics": {"model": "fermi_dirac"},
        "bandgap_narrowing": {
            "model": "old_slotboom",
            "reference_doping_cm3": 1.0e17,
            "coefficient_eV": 0.009,
            "smoothing": 0.5,
            "offset_eV": 0.0,
            "fermi_statistics_correction": True,
        },
        "recombination": {
            "mechanisms": ["srh", "auger"],
            "taun_s": 1.0e-5,
            "taup_s": 3.0e-6,
            "auger_cn_m6_per_s": 2.9e-43,
            "auger_cp_m6_per_s": 1.028e-43,
            "srh_doping_dependence": {
                "enabled": True,
                "concentration_basis": "total_impurity",
                "density_coupling": "sentaurus_default",
                "temperature_dependence": True,
                "temperature_K": 300.0,
                "reference_temperature_K": 300.0,
                "electron_temperature_exponent": -1.5,
                "hole_temperature_exponent": -1.5,
                "electron": {"tau_min_s": 0.0, "tau_max_s": 1.0e-5,
                             "reference_doping_cm3": 1.0e16, "gamma": 1.0},
                "hole": {"tau_min_s": 0.0, "tau_max_s": 3.0e-6,
                         "reference_doping_cm3": 1.0e16, "gamma": 1.0},
            },
        },
        "mobility": {
            "doping_concentration_basis": "net_doping",
            "high_field_driving_force": "quasi_fermi_gradient",
            "high_field_gradient_discretization": "edge_projection",
            "electron": {
                "mu_const_cm2_per_V_s": 1417.0,
                "mu_min_cm2_per_V_s": 52.2,
                "mu_min2_cm2_per_V_s": 52.2,
                "mu1_cm2_per_V_s": 43.4,
                "pc_cm3": 0.0,
                "reference_doping_cm3": 9.68e16,
                "cs_cm3": 3.43e20,
                "doping_exponent": 0.68,
                "masetti_beta": 2.0,
                "saturation_velocity_cm_per_s": 1.07e7,
                "high_field_exponent": 1.109,
            },
            "hole": {
                "mu_const_cm2_per_V_s": 470.5,
                "mu_min_cm2_per_V_s": 44.9,
                "mu_min2_cm2_per_V_s": 0.0,
                "mu1_cm2_per_V_s": 29.0,
                "pc_cm3": 9.23e16,
                "reference_doping_cm3": 2.23e17,
                "cs_cm3": 6.10e20,
                "doping_exponent": 0.719,
                "masetti_beta": 2.0,
                "saturation_velocity_cm_per_s": 8.37e6,
                "high_field_exponent": 1.213,
            },
        },
    }


class Phase23DeckTest(unittest.TestCase):
    def test_polysi_mapping_uses_the_sealed_boundary_value(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mesh = root / "mesh.json"
            mesh.write_text(json.dumps({"contacts": [{
                "name": "gate", "region_id": 7, "node_ids": [2, 4]
            }]}), encoding="utf-8")
            fields = root / "export" / "fields"
            fields.mkdir(parents=True)
            with (fields / "ElectrostaticPotential_region7.csv").open(
                    "w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=["node_id", "component0"])
                writer.writeheader()
                writer.writerows([
                    {"node_id": 2, "component0": "0.5609636105202798"},
                    {"node_id": 4, "component0": "0.5609636105202798"},
                ])
            mapping = derive_polysi_flatband(mesh, root / "export")
            self.assertAlmostEqual(mapping["vela_flatband_voltage_V"], -0.5609636105202798)
            self.assertEqual(mapping["node_count"], 2)
            self.assertEqual(mapping["spread_V"], 0.0)

    def test_polysi_mapping_rejects_a_nonconstant_gate(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "mesh.json").write_text(json.dumps({"contacts": [{
                "name": "gate", "region_id": 1, "node_ids": [0, 1]
            }]}), encoding="utf-8")
            fields = root / "export" / "fields"
            fields.mkdir(parents=True)
            (fields / "ElectrostaticPotential_region1.csv").write_text(
                "node_id,component0\n0,0.5\n1,0.6\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "not constant"):
                derive_polysi_flatband(root / "mesh.json", root / "export")

    def test_classical_solver_materializes_units_and_physics_layer(self) -> None:
        solver = classical_solver(physics_contract(), high_field=True)
        self.assertEqual(solver["mobility"]["model"], "masetti_field")
        self.assertAlmostEqual(solver["mobility"]["electron_mumin1_m2_V_s"], 0.00522)
        self.assertAlmostEqual(
            solver["mobility"]["electron_saturation_velocity_m_s"], 1.07e5)
        self.assertEqual(solver["impact_ionization"]["model"], "none")
        self.assertEqual(solver["line_search_mode"], "merit")
        self.assertEqual(solver["damping_factor"], 1.0)
        self.assertNotIn("damping_psi", solver)
        self.assertEqual(solver["recombination"], ["srh", "auger"])
        self.assertTrue(solver["srh_doping_dependence"]["enabled"])
        self.assertEqual(
            solver["srh_doping_dependence"]["electron"]["reference_doping_m3"],
            1.0e22,
        )

    def test_currentplot_points_must_be_strictly_increasing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            curve = Path(directory) / "curve.csv"
            curve.write_text("bias_V,current_A_per_um\n0,0\n0.5,1\n1,2\n", encoding="utf-8")
            self.assertEqual(exact_bias_points(curve), [0.0, 0.5, 1.0])
            curve.write_text("bias_V,current_A_per_um\n0,0\n0,1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "strictly increasing"):
                exact_bias_points(curve)

    def test_phase23_zero_bias_handoff_uses_poisson_and_near_frozen_qf(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stage1 = root / "stage1"
            exact = stage1 / "vela_exact_topology"
            qualification = stage1 / "qualification"
            oracle = root / "oracle"
            contracts = root / "contracts"
            output = root / "output"
            for path in (exact, qualification / "sentaurus_eq_0v_export" / "fields",
                         oracle, contracts):
                path.mkdir(parents=True, exist_ok=True)
            (exact / "mesh.json").write_text(json.dumps({
                "contacts": [{"name": "gate", "region_id": 7, "node_ids": [0]}]
            }), encoding="utf-8")
            (exact / "doping.csv").write_text("node_id,net_doping_m3\n0,0\n",
                                                encoding="utf-8")
            (qualification / "sentaurus_eq_0v_export" / "fields" /
             "ElectrostaticPotential_region7.csv").write_text(
                "node_id,component0\n0,0.5\n", encoding="utf-8")
            (qualification / "sentaurus_eq_0v_state.csv").write_text(
                "node_id,psi,phin,phip,n,p\n", encoding="utf-8")
            (oracle / "IdVg_n2_des_drain_curve.csv").write_text(
                "bias_V,current_A_per_um\n0,0\n1,1\n", encoding="utf-8")
            (contracts / "physics_contract.json").write_text(json.dumps({
                "materials_file": "materials.json",
                **physics_contract(),
            }), encoding="utf-8")
            (contracts / "materials.json").write_text("{}\n", encoding="utf-8")

            manifest = prepare_phase23(stage1, oracle, contracts, output)
            poisson = json.loads(Path(
                manifest["decks"]["g_contact_polysi_poisson_eq"]
            ).read_text())
            equilibrium = json.loads(Path(
                manifest["decks"]["g_contact_polysi_eq"]
            ).read_text())
            repeat = json.loads(Path(
                manifest["decks"]["g_contact_polysi_eq_repeat"]
            ).read_text())
            prebias = json.loads(Path(
                manifest["decks"]["g3_drain_prebias"]
            ).read_text())
            sentaurus_path = json.loads(Path(
                manifest["decks"]["g3_drain_prebias_sentaurus_path"]
            ).read_text())
            prebias_repeat = json.loads(Path(
                manifest["decks"]["g3_drain_prebias_repeat"]
            ).read_text())
            idvg = json.loads(Path(
                manifest["decks"]["g3_idvg"]
            ).read_text())

            self.assertEqual(poisson["solver"]["method"], "poisson_only")
            self.assertTrue(equilibrium["sweep"]["initial_state_file"].endswith(
                "g_contact_polysi_poisson_eq_state.csv"))
            self.assertEqual(
                equilibrium["solver"]["quasi_fermi_update_limit_V"], 1.0e-12)
            self.assertEqual(
                repeat["solver"]["quasi_fermi_update_limit_V"], 1.0e-12)
            self.assertEqual(
                prebias["solver"]["quasi_fermi_update_limit_V"], 0.1)
            self.assertEqual(
                sentaurus_path["sweep"]["bias_points"],
                SENTAURUS_G3_DRAIN_PREBIAS_POINTS_V,
            )
            self.assertEqual(prebias["sweep"]["bias_points"], [0.0, 0.1])
            self.assertTrue(
                prebias_repeat["sweep"]["initial_state_file"].endswith(
                    "g3_drain_prebias_state.csv"
                )
            )
            self.assertTrue(idvg["sweep"]["initial_state_file"].endswith(
                "g3_drain_prebias_repeat_state.csv"
            ))

    def test_wp15_matrix_preserves_physics_and_exposes_solver_controls(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            base = {
                "simulation_type": "dc_sweep",
                "mesh_file": "mesh.json",
                "node_doping_file": "doping.csv",
                "materials_file": "materials.json",
                "output_csv": "base.csv",
                "contacts": [{"name": "drain", "type": "ohmic", "bias": 0.0}],
                "solver": {
                    "method": "newton",
                    "line_search_mode": "block_filter",
                    "stall_residual_floor": 1.0e-7,
                    "mobility": {"model": "masetti"},
                },
                "sweep": {
                    "bias_points": [0.0],
                    "initial_state_file": "state.csv",
                    "write_state_file": "base_state.csv",
                    "diagnostics": {
                        "terminal_balance": {"enabled": True, "csv_file": "terminal.csv"},
                        "srh_balance": {"enabled": True, "csv_file": "srh.csv"},
                    },
                },
            }
            base_path = root / "base.json"
            base_path.write_text(json.dumps(base), encoding="utf-8")
            paths = prepare_wp15_diagnostics(base_path, root / "matrix")
            baseline = json.loads(Path(paths["baseline_block_filter"]).read_text())
            merit = json.loads(Path(paths["merit"]).read_text())
            guard = json.loads(
                Path(paths["merit_contact_basin_reclose_guard"]).read_text()
            )
            frozen_qf = json.loads(
                Path(paths["poisson_handoff_newton_qf_1e_12"]).read_text()
            )
            frozen_qf_repeat = json.loads(
                Path(paths["poisson_handoff_newton_qf_1e_12_repeat"]).read_text()
            )
            frozen_qf_guard = json.loads(
                Path(paths["poisson_handoff_newton_qf_1e_12_guard"]).read_text()
            )
            last_bias_repeat = json.loads(
                Path(paths["merit_contact_basin_last_bias_repeat"]).read_text()
            )
            unscaled_repeat = json.loads(Path(
                paths["merit_contact_basin_no_row_scaling_repeat"]
            ).read_text())
            self.assertEqual(baseline["solver"]["mobility"], merit["solver"]["mobility"])
            self.assertEqual(baseline["contacts"], merit["contacts"])
            self.assertEqual(merit["solver"]["line_search_mode"], "merit")
            self.assertEqual(guard["solver"]["max_iter"], 0)
            self.assertEqual(guard["solver"]["quasi_fermi_reference"], "contact_basin")
            self.assertEqual(
                frozen_qf["solver"]["quasi_fermi_update_limit_V"], 1.0e-12
            )
            self.assertTrue(
                frozen_qf["sweep"]["initial_state_file"].endswith(
                    "poisson_only_state.csv"
                )
            )
            self.assertTrue(
                frozen_qf_repeat["sweep"]["initial_state_file"].endswith(
                    "poisson_handoff_newton_qf_1e_12_state.csv"
                )
            )
            self.assertEqual(frozen_qf_guard["solver"]["max_iter"], 0)
            self.assertEqual(last_bias_repeat["sweep"]["bias_points"], [0.0])
            self.assertTrue(
                last_bias_repeat["sweep"]["initial_state_file"].endswith(
                    "merit_contact_basin_state.csv"
                )
            )
            self.assertFalse(
                unscaled_repeat["solver"]["continuity_row_scaling"]["enabled"]
            )


class SentaurusAblationTest(unittest.TestCase):
    def test_chain_changes_exactly_one_declared_factor_per_parent(self) -> None:
        source_text = '''Electrode { { Name="gate" Material="PolySi"(N) Voltage=0 } }
Physics (Material="Silicon") {
  eQuantumPotential(density) hQuantumPotential(density)
  Mobility( HighFieldSaturation Enormal (IALMob(AutoOrientation)) )
}
Solve {
	Coupled { Poisson Electron Hole }
	Quasistationary(Goal { Name="drain" Voltage=0.1 }) {
		Coupled { Poisson Electron Hole }	}
}
'''
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "IdVg.cmd"
            source.write_text(source_text, encoding="utf-8")
            output = root / "out"
            manifest = prepare_ablation(source, output, -0.5609636105202798)
            records = {item["id"]: item for item in manifest["single_factor_chain"]}
            self.assertEqual(records["G4-no-highfield"]["parent"], "G3-no-IALMob")
            g3 = (output / "G3-no-IALMob.cmd").read_text(encoding="utf-8")
            self.assertNotIn("QuantumPotential", g3)
            self.assertNotIn("IALMob", g3)
            self.assertIn("HighFieldSaturation", g3)
            g4 = (output / "G4-no-highfield.cmd").read_text(encoding="utf-8")
            self.assertNotIn("HighFieldSaturation", g4)
            equilibrium = (output / "G4-equilibrium-capture.cmd").read_text(
                encoding="utf-8")
            self.assertNotIn("Quasistationary", equilibrium)
            self.assertIn('FilePrefix="G4_equilibrium_state"', equilibrium)
            self.assertEqual(
                records["G4-equilibrium-capture"]["parent"], "G4-no-highfield")
            control = (output / "G-contact-poly-barrier-control.cmd").read_text(
                encoding="utf-8")
            self.assertIn("Barrier= -0.560963610520279", control)
            self.assertNotIn('Material="PolySi"(N)', control)

    def test_missing_factor_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "IdVg.cmd"
            source.write_text("Physics { Fermi }\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "expected exactly one match"):
                prepare_ablation(source, root / "out", -0.5)


class Phase23SummaryMathTest(unittest.TestCase):
    def test_percentile_is_linear_between_order_statistics(self) -> None:
        self.assertEqual(percentile([0.0, 10.0], 0.95), 9.5)

    def test_field_metrics_use_common_semiconductor_nodes(self) -> None:
        reference = {
            0: {"psi": 0.0, "electrons_m3": 1.0e20, "holes_m3": 1.0e10},
            1: {"psi": 0.1, "electrons_m3": 1.0e19, "holes_m3": 1.0e11},
            2: {"psi": 9.0, "electrons_m3": 0.0, "holes_m3": 0.0},
        }
        candidate = {
            0: {"psi": 0.001, "electrons_m3": 1.0e20, "holes_m3": 1.0e10},
            1: {"psi": 0.102, "electrons_m3": 1.0e18, "holes_m3": 1.0e12},
            2: {"psi": 0.0, "electrons_m3": 0.0, "holes_m3": 0.0},
        }
        metrics = field_metrics(reference, candidate)
        self.assertEqual(metrics["semiconductor_node_count"], 2)
        self.assertAlmostEqual(metrics["psi_abs_error_V"]["median"], 0.0015)
        self.assertAlmostEqual(metrics["electron_abs_error_dex"]["max"], 1.0)


class SentaurusAblationSummaryTest(unittest.TestCase):
    def test_comparison_requires_the_exact_same_bias_grid(self) -> None:
        with self.assertRaisesRegex(ValueError, "bias grids differ"):
            compare_ablation([(0.0, 1.0), (1.0, 2.0)],
                             [(0.0, 1.0), (1.1, 2.0)])

    def test_comparison_reports_strong_inversion_without_curve_interpolation(self) -> None:
        metrics = compare_ablation(
            [(0.0, 1.0e-12), (4.0, 1.0), (5.0, 2.0)],
            [(0.0, 1.0e-12), (4.0, 1.1), (5.0, 1.8)],
        )
        strong = metrics["strong_inversion_relative_change_percent"]
        self.assertEqual(strong["count"], 2)
        self.assertAlmostEqual(strong["median"], 10.0)

    def test_crossing_and_timing_are_diagnostic_helpers(self) -> None:
        self.assertAlmostEqual(
            crossing([(0.0, 1.0e-12), (1.0, 1.0e-8)], 1.0e-10), 0.5)
        with tempfile.TemporaryDirectory() as directory:
            timing = Path(directory) / "timing.txt"
            timing.write_text(
                "Elapsed (wall clock) time (h:mm:ss or m:ss): 6:43.91\n",
                encoding="utf-8",
            )
            self.assertAlmostEqual(wall_seconds(timing), 403.91)


if __name__ == "__main__":
    unittest.main()
