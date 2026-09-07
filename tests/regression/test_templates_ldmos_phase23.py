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
from audit_templates_ldmos_fermi_bgn_mapping import (  # noqa: E402
    summarize_variant,
    variants as fermi_bgn_variants,
    write_reference_aligned_state,
)
from audit_templates_ldmos_g3_continuity_sg_contact import (  # noqa: E402
    Q_C,
    contact_summary,
    reconstruct_divergence,
    sentaurus_plt_current,
    sg_contact_cut,
    wp15_block_reclose_config,
)
from prepare_templates_ldmos_wp15_diagnostics import (  # noqa: E402
    prepare as prepare_wp15_diagnostics,
)
from prepare_templates_ldmos_sentaurus_ablation import (  # noqa: E402
    prepare as prepare_ablation,
)
from prepare_templates_ldmos_idvd_ablation import (  # noqa: E402
    prepare as prepare_idvd_ablation,
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
from summarize_templates_ldmos_idvd_ablation import (  # noqa: E402
    compare_curve as compare_idvd_curve,
    hrec_decision,
)
from run_templates_ldmos_stage4_d5 import (  # noqa: E402
    make_drain_zero_prebias,
    make_gate_prebias_vd0,
    make_gate8_prebias,
    make_idvd,
)
from analyze_templates_ldmos_stage4_d5 import (  # noqa: E402
    curve_error as stage4_curve_error,
    kcl_audit as stage4_kcl_audit,
    ratio_error as stage4_ratio_error,
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
            "sentaurus_dEg0_eV": -0.01595,
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
            "g3_model": "constant_field",
            "doping_dependence_enabled": False,
            "high_field_saturation_enabled": True,
            "doping_concentration_basis": "net_doping",
            "high_field_driving_force": "quasi_fermi_gradient",
            "high_field_gradient_discretization": "edge_projection",
            "contact_electric_field_fallback": True,
            "contact_electric_field_fallback_scope": "contact_node_cell",
            "contact_electric_field_fallback_mode": "cell_gradient_magnitude",
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
    def test_fermi_bgn_audit_matrix_is_single_factor_explicit(self) -> None:
        base = {
            "carrier_statistics": {"model": "fermi_dirac"},
            "bandgap_narrowing": {
                "model": "old_slotboom",
                "coefficient_eV": 0.009,
                "offset_eV": 0.0,
                "fermi_statistics_correction": True,
            },
        }
        matrix = dict(fermi_bgn_variants(base))
        self.assertEqual(len(matrix), 6)
        self.assertTrue(matrix["fermi_oldslotboom_correction"]
                        ["bandgap_narrowing"]["fermi_statistics_correction"])
        self.assertFalse(matrix["fermi_oldslotboom_no_correction"]
                         ["bandgap_narrowing"]["fermi_statistics_correction"])
        self.assertEqual(matrix["fermi_no_bgn"]["bandgap_narrowing"],
                         {"model": "none"})
        self.assertEqual(matrix["boltzmann_oldslotboom"]
                         ["carrier_statistics"]["model"], "boltzmann")

    def test_fermi_bgn_audit_summarizes_heavy_majority_mapping(self) -> None:
        def row(node: int, doping: float, electron_error: float,
                hole_error: float) -> dict[str, str]:
            return {
                "node_id": str(node), "net_doping_m3": str(doping),
                "constrained": "0", "production_residual": "2.0",
                "input_electron_density_m3": "100",
                "reconstructed_electron_density_m3": "10",
                "input_hole_density_m3": "100",
                "reconstructed_hole_density_m3": "10",
                "electron_qf_mapping_error_V": str(electron_error),
                "hole_qf_mapping_error_V": str(hole_error),
            }

        summary = summarize_variant(
            [row(1, 2.0e25, 0.01, 0.5), row(2, -3.0e25, 0.6, 0.02)],
            {1, 2},
        )
        self.assertEqual(summary["heavy_majority_density_abs_error_dex"]
                         ["median"], 1.0)
        self.assertAlmostEqual(
            summary["heavy_majority_qf_mapping_abs_error_V"]["median"],
            0.015,
        )

    def test_fermi_bgn_reference_alignment_uses_half_deg0_band_centre_shift(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "state.csv"
            destination = root / "aligned.csv"
            source.write_text(
                "node_id,psi,phin,phip,electrons_m3,holes_m3\n"
                "0,0.1,0.02,-0.03,1e20,2e10\n",
                encoding="utf-8")
            transform = write_reference_aligned_state(
                source, destination, -0.01595)
            with destination.open(newline="", encoding="utf-8") as stream:
                row = next(csv.DictReader(stream))
            self.assertAlmostEqual(float(row["phin"]), 0.027975)
            self.assertAlmostEqual(float(row["phip"]), -0.037975)
            self.assertAlmostEqual(transform["electron_qf_offset_V"], 0.007975)
            self.assertAlmostEqual(transform["hole_qf_offset_V"], -0.007975)

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
        self.assertEqual(solver["mobility"]["model"], "constant_field")
        self.assertAlmostEqual(
            solver["mobility"]["electron_saturation_velocity_m_s"], 1.07e7)
        self.assertAlmostEqual(
            solver["mobility"]["hole_saturation_velocity_m_s"], 8.37e6)
        self.assertTrue(
            solver["mobility"]["contact_electric_field_fallback"])
        self.assertEqual(
            solver["mobility"]["contact_electric_field_fallback_scope"],
            "contact_node_cell")
        self.assertNotIn("electron_mumin1_m2_V_s", solver["mobility"])
        self.assertNotIn("electron_cr_m3", solver["mobility"])
        low_field = classical_solver(physics_contract(), high_field=False)
        self.assertEqual(low_field["mobility"], {"model": "constant"})
        self.assertEqual(
            low_field["contact_boundary_reconstruction"], "legacy_node_local")
        self.assertEqual(solver["impact_ionization"]["model"], "none")
        self.assertEqual(solver["line_search_mode"], "block_filter")
        self.assertEqual(solver["quasi_fermi_reference"], "contact_basin")
        self.assertEqual(
            solver["contact_boundary_reconstruction"], "legacy_node_local")
        self.assertEqual(solver["damping_factor"], 1.0)
        self.assertEqual(solver["reltol"], 1.0e-10)
        self.assertEqual(solver["abstol"], 1.0e-14)
        self.assertEqual(solver["stall_residual_floor"], 1.0e-12)
        self.assertEqual(
            solver["linear_equilibration"], {"mode": "l2_row_column"})
        self.assertEqual(solver["block_absolute_convergence"], {
            "mode": "enforce",
            "psi_residual_ceiling": 5.0e-8,
            "electron_residual_ceiling": 1.0e-11,
            "hole_residual_ceiling": 3.0e-10,
        })
        self.assertEqual(
            solver["contact_majority_qf_branch_guard_contacts"], ["drain"])
        self.assertEqual(
            solver["contact_majority_qf_branch_drop_limit_V"], 5.0e-11)
        self.assertNotIn("damping_psi", solver)
        self.assertEqual(solver["recombination"], ["srh", "auger"])
        self.assertTrue(solver["srh_doping_dependence"]["enabled"])
        self.assertEqual(
            solver["bandgap_narrowing"]["reference_doping_m3"],
            1.0e17,
        )
        self.assertEqual(
            solver["srh_doping_dependence"]["electron"]["reference_doping_m3"],
            1.0e16,
        )
        self.assertEqual(
            solver["srh_doping_dependence"]["hole"]["reference_doping_m3"],
            1.0e16,
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
                "node_id,psi,phin,phip,n,p\n"
                "0,0.1,0.02,-0.03,1e20,2e10\n", encoding="utf-8")
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
            checkpoint = json.loads(Path(
                manifest["decks"]["g3_drain_prebias_vd0p023_checkpoint"]
            ).read_text())
            prebias_repeat = json.loads(Path(
                manifest["decks"]["g3_drain_prebias_repeat"]
            ).read_text())
            idvg = json.loads(Path(
                manifest["decks"]["g3_idvg"]
            ).read_text())

            self.assertEqual(poisson["solver"]["method"], "poisson_only")
            self.assertTrue(poisson["sweep"]["initial_state_file"].endswith(
                "sentaurus_eq_0v_state.csv"))
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
            self.assertEqual(
                checkpoint["sweep"]["bias_points"],
                SENTAURUS_G3_DRAIN_PREBIAS_POINTS_V[:9],
            )
            self.assertEqual(
                checkpoint["sweep"]["bias_points"][-1],
                0.0231559774221138,
            )
            self.assertNotIn("predictor", checkpoint["sweep"])
            self.assertTrue(
                checkpoint["sweep"]["write_state_file"].endswith(
                    "g3_drain_prebias_vd0p023_checkpoint_state.csv"
                )
            )
            self.assertEqual(prebias["sweep"]["bias_points"], [0.0, 0.1])
            self.assertTrue(
                prebias_repeat["sweep"]["initial_state_file"].endswith(
                    "g3_drain_prebias_sentaurus_path_state.csv"
                )
            )
            self.assertTrue(idvg["sweep"]["initial_state_file"].endswith(
                "g3_idvg_seed_repeat_state.csv"
            ))
            for deck in (prebias, sentaurus_path, checkpoint, prebias_repeat):
                self.assertNotIn("predictor", deck["sweep"])
                self.assertEqual(
                    deck["solver"]["quasi_fermi_reference"],
                    "contact_basin",
                )
                self.assertEqual(
                    deck["solver"]["line_search_mode"], "block_filter")
                self.assertEqual(
                    deck["solver"]["block_absolute_convergence"]["mode"],
                    "enforce",
                )
                self.assertEqual(
                    deck["solver"]["contact_majority_qf_branch_guard_contacts"],
                    ["drain"],
                )
                self.assertTrue(
                    deck["sweep"]["diagnostics"]["newton_history"]["enabled"])
            self.assertNotIn("predictor", idvg["sweep"])
            self.assertEqual(
                idvg["solver"]["block_absolute_convergence"]["mode"],
                "enforce",
            )
            self.assertNotIn(
                "contact_majority_qf_branch_drop_limit_V", idvg["solver"])
            self.assertNotIn(
                "contact_majority_qf_branch_guard_contacts", idvg["solver"])
            self.assertEqual(
                manifest["wp15_contract"]["contact_majority_qf_branch_guard"],
                "deep_off_seed_and_drain_prebias_only",
            )
            self.assertEqual(
                manifest["wp15_contract"]["linear_equilibration"],
                "l2_row_column_on_all_G3_decks",
            )
            self.assertEqual(
                manifest["wp15_contract"]["quasi_fermi_reference"],
                "contact_basin_on_all_G3_decks",
            )
            self.assertEqual(
                manifest["wp15_contract"]["line_search_mode"],
                "block_filter_on_all_G3_decks",
            )
            self.assertEqual(
                manifest["wp15_contract"]["contact_boundary_reconstruction"],
                "legacy_node_local_for_multipolarity_source_short",
            )

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

    def test_idvd_chain_is_strictly_cumulative_and_single_factor(self) -> None:
        source_text = '''Electrode {
 { Name= "drain" Voltage= 0 hRecVelocity= 1.93E6 }
 { Name= "source" Voltage= 0 hRecVelocity= 1.93E6 }
}
Thermode {
 { Name= "th_lat" Temperature= 300 }
}
Physics(Material="Silicon") {
 eQuantumPotential(density) hQuantumPotential(density)
 Mobility(HighFieldSaturation Enormal (IALMob(AutoOrientation)))
}
Solve {
 Coupled { Poisson Electron Hole Temperature }
 Coupled { Poisson Electron Hole Temperature }
 Coupled { Poisson Electron Hole Temperature }
 Coupled { Poisson Electron Hole Temperature }
}
'''
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "IdVd.cmd"
            parameter = root / "sdevice.par"
            source.write_text(source_text, encoding="utf-8")
            parameter.write_text("parameter\n", encoding="utf-8")
            manifest = prepare_idvd_ablation(source, parameter, root / "out")
            records = {item["id"]: item for item in manifest["single_factor_chain"]}
            self.assertEqual(records["D2-no-hRecVelocity"]["parent"], "D1-isothermal")
            d1 = (root / "out" / "D1-isothermal" / "IdVd.cmd").read_text()
            self.assertNotIn("Thermode", d1)
            self.assertNotIn("Hole Temperature", d1)
            self.assertIn("hRecVelocity", d1)
            d2 = (root / "out" / "D2-no-hRecVelocity" / "IdVd.cmd").read_text()
            self.assertNotIn("hRecVelocity", d2)
            self.assertIn("hQuantumPotential", d2)
            d5 = (root / "out" / "D5-no-IALMob" / "IdVd.cmd").read_text()
            self.assertNotIn("QuantumPotential", d5)
            self.assertNotIn("IALMob", d5)
            self.assertIn("HighFieldSaturation", d5)
            self.assertEqual(
                (root / "out" / "D5-no-IALMob" / "sdevice.par").read_text(),
                "parameter\n",
            )


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


class G3ContinuitySgContactAuditTest(unittest.TestCase):
    def test_wp15_reclose_uses_block_contract_and_drain_scoped_guard(self) -> None:
        baseline = {
            "solver": {},
            "contacts": [
                {"name": "gate"}, {"name": "drain"},
                {"name": "source"}, {"name": "substrate"},
            ],
        }
        config = wp15_block_reclose_config(
            baseline, Path("state.csv"), Path("out"), 0.0231559774221138
        )
        solver = config["solver"]
        self.assertEqual(solver["block_absolute_convergence"]["mode"], "enforce")
        self.assertEqual(
            solver["contact_majority_qf_branch_guard_contacts"], ["drain"]
        )
        self.assertEqual(
            solver["contact_majority_qf_branch_drop_limit_V"], 5.0e-11
        )
        self.assertNotIn("predictor", config["sweep"])

    def test_sg_cut_uses_production_contact_orientation(self) -> None:
        rows = [
            {
                "node0": "1", "node1": "2",
                "electron_particle_line_flux_per_m_s": "3.0",
                "hole_particle_line_flux_per_m_s": "1.0",
            },
            {
                "node0": "3", "node1": "1",
                "electron_particle_line_flux_per_m_s": "5.0",
                "hole_particle_line_flux_per_m_s": "2.0",
            },
            {
                "node0": "2", "node1": "3",
                "electron_particle_line_flux_per_m_s": "100.0",
                "hole_particle_line_flux_per_m_s": "100.0",
            },
        ]
        result = sg_contact_cut(rows, {1})
        self.assertEqual(result["crossing_edge_count"], 2)
        self.assertAlmostEqual(result["electron_A_per_um"], 2.0 * Q_C / 1.0e6)
        self.assertAlmostEqual(result["hole_A_per_um"], 1.0 * Q_C / 1.0e6)
        self.assertAlmostEqual(result["total_A_per_um"], 1.0 * Q_C / 1.0e6)

    def test_edge_divergence_matches_node0_plus_node1_minus_contract(self) -> None:
        rows = [
            {"node0": "0", "node1": "1", "electron_flux": "2.5"},
            {"node0": "1", "node1": "2", "electron_flux": "-1.0"},
        ]
        self.assertEqual(
            reconstruct_divergence(rows, "electron"),
            {0: 2.5, 1: -3.5, 2: 1.0},
        )

    def test_contact_summary_separates_stable_and_drift_diffusion_currents(self) -> None:
        row = {
            "current_contact": "drain", "edge_id": "7", "node0": "1", "node1": "2",
            "current_electron": "2e-6", "current_electron_long_double_reference": "2e-6",
            "current_electron_drift": "10", "current_electron_diffusion": "-9.999998",
            "current_hole": "0", "current_total": "2e-6", "phin0": "0.1",
            "phin1": "0.100000001", "psi0": "0.2", "psi1": "0.3",
            "n0": "1e20", "n1": "2e20",
        }
        result = contact_summary([row], "drain", {1})
        self.assertEqual(result["edge_count"], 1)
        self.assertAlmostEqual(result["electron_A_per_um"], 2.0e-12)
        self.assertGreater(result["drift_diffusion_cancellation_condition"], 1.0e6)

    def test_dfise_reference_selects_highest_time_at_exact_bias(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reference.plt"
            path.write_text(
                '''DF-ISE text
Info { datasets = [ "time" "drain OuterVoltage" "drain eCurrent"
  "drain hCurrent" "drain TotalCurrent" ] }
Data {
  0 2.31559774221138E-02 1.0E-15 2.0E-20 1.00002E-15
  1 2.31559774221138E-02 2.0E-15 3.0E-20 2.00003E-15
}
''',
                encoding="utf-8",
            )
            result = sentaurus_plt_current(path, 0.0231559774221138)
            self.assertEqual(result["matching_rows"], 2)
            self.assertEqual(result["time"], 1.0)
            self.assertAlmostEqual(result["drain_total_A_per_um"], 2.00003e-15)


class SentaurusAblationSummaryTest(unittest.TestCase):
    def test_stage4_kcl_uses_maximum_terminal_current_and_keeps_zero_bias(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "balance.csv"
            with path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream)
                writer.writerow(["point_index", "bias_V", "contact", "current_total_A_per_um"])
                for index, bias, currents in [(0, 0.0, [3e-26, 2e-26, 0.0, 0.0]),
                                              (1, 1.0, [-0.99, 1.0, 0.0, 0.0])]:
                    for contact, current in zip(("source", "drain", "gate", "substrate"), currents):
                        writer.writerow([index, bias, contact, current])
            result = stage4_kcl_audit(path)
            # Sum-absolute normalization would understate the 1 V error by ~2x.
            self.assertAlmostEqual(result["max_nonzero_bias_normalized_kcl_percent"], 1.0)
            zero = result["zero_bias_points"][0]
            self.assertAlmostEqual(zero["absolute_kcl_A_per_um"] / 1e-26, 5.0)
            self.assertAlmostEqual(zero["normalized_kcl_percent"], 500.0 / 3.0)
            self.assertEqual(result["max_normalized_kcl_percent"], zero["normalized_kcl_percent"])

    def test_stage4_kcl_rejects_missing_or_nonfinite_terminal_data(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "balance.csv"
            for currents in ([1.0, -1.0, 0.0], [1.0, float("nan"), 0.0, 0.0]):
                with self.subTest(currents=currents):
                    with path.open("w", newline="", encoding="utf-8") as stream:
                        writer = csv.writer(stream)
                        writer.writerow(["point_index", "bias_V", "contact", "current_total_A_per_um"])
                        for contact, current in zip(("drain", "source", "gate", "substrate"), currents):
                            writer.writerow([0, 1.0, contact, current])
                    with self.assertRaises(ValueError):
                        stage4_kcl_audit(path)

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

    def test_idvd_curve_metrics_and_hrec_gate_use_exact_points(self) -> None:
        parent = [(0.0, 0.0), (1.0, 1.0), (40.0, 2.0)]
        child = [(0.0, 0.0), (1.0, 1.005), (40.0, 2.01)]
        metric = compare_idvd_curve(parent, child)
        self.assertAlmostEqual(
            metric["low_vd_differential_resistance_change_percent"],
            abs(1.0 / 1.005 - 1.0) * 100.0,
        )
        decision = hrec_decision({
            "Vg4": metric,
            "Vg8": metric,
            "gate_ratio_change_percent": {"count": 2, "median": 0.0,
                                            "p95": 0.0, "max": 0.0},
        })
        self.assertFalse(decision["implementation_required_by_curve"])
        with self.assertRaisesRegex(ValueError, "exact Id-Vd grids differ"):
            compare_idvd_curve(parent, [(0.0, 0.0), (1.1, 1.0), (40.0, 2.0)])

    def test_stage4_d5_configs_preserve_qualified_contract(self) -> None:
        base = {
            "solver": {"impact_ionization": {"model": "none"}},
            "contacts": [
                {"name": "gate", "bias": 0.0}, {"name": "drain", "bias": 0.1},
                {"name": "source", "bias": 0.0}, {"name": "substrate", "bias": 0.0},
            ],
            "sweep": {"bias_points": [0.0], "diagnostics": {}},
            "mesh_geometry": {
                "carrier_transport_couple_profile": "templates_ldmos_external_averagebox"
            },
            "discretization": {"poisson_charge_volume_policy": "material_local"},
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gate8 = make_gate8_prebias(base, root / "vg5.csv", root)
            self.assertEqual(gate8["sweep"]["bias_points"][-1], 8.0)
            self.assertNotIn("predictor", gate8["sweep"])
            gate_vd0 = make_gate_prebias_vd0(base, root / "eq.csv", root)
            self.assertEqual(gate_vd0["sweep"]["bias_points"][8], 4.0)
            self.assertEqual(gate_vd0["sweep"]["bias_points"][-1], 8.0)
            self.assertEqual(next(item for item in gate_vd0["contacts"]
                                  if item["name"] == "drain")["bias"], 0.0)
            idvd = make_idvd(base, 4.0, root / "vg4.csv", [0.0, 1.0, 2.0], root)
            self.assertEqual(idvd["sweep"]["contact"], "drain")
            self.assertTrue({0.0, 0.01, 0.8, 1.0, 2.0}.issubset(
                set(idvd["sweep"]["bias_points"])))
            self.assertEqual(idvd["sweep"]["min_step"], 1.0e-3)
            self.assertEqual(idvd["sweep"]["max_retries"], 12)
            self.assertEqual(idvd["sweep"]["initial_step"], 2.5e-3)
            self.assertEqual(idvd["sweep"]["growth_factor"], 1.0)
            self.assertTrue(
                idvd["solver"]["quasi_fermi_recenter_on_initial_state"]
            )
            self.assertFalse(
                idvd["solver"]["mobility"]["jacobian_field_derivatives"]
            )
            self.assertEqual(next(item for item in idvd["contacts"]
                                  if item["name"] == "gate")["bias"], 4.0)
            self.assertEqual(
                idvd["solver"]["contact_majority_qf_branch_guard_contacts"],
                ["drain"],
            )
            self.assertEqual(
                idvd["solver"]["contact_majority_qf_branch_drop_limit_V"],
                5.0e-11,
            )
            zero = make_drain_zero_prebias(base, 4.0, root / "vg4.csv", root)
            self.assertEqual(zero["sweep"]["bias_points"],
                             [0.1, 0.075, 0.05, 0.025, 0.0])
            self.assertLess(zero["sweep"]["step"], 0.0)

    def test_stage4_d5_metrics_use_exact_grid_and_endpoint_gate_ratio(self) -> None:
        reference = [(0.0, 0.0), (1.0, 1.0), (40.0, 2.0)]
        candidate = [(0.0, 0.0), (1.0, 1.1), (40.0, 2.2)]
        metric = stage4_curve_error(reference, candidate)
        self.assertAlmostEqual(metric["relative_error_percent"]["median"], 10.0)
        ratios = stage4_ratio_error(
            {"Vg4": reference, "Vg8": [(0.0, 0.0), (1.0, 2.0), (40.0, 4.0)]},
            {"Vg4": candidate, "Vg8": [(0.0, 0.0), (1.0, 2.42), (40.0, 4.84)]},
        )
        self.assertAlmostEqual(ratios["endpoint"], 10.0)
        with_bridge = [(0.0, 0.0), (0.1, 0.02), (1.0, 1.1), (40.0, 2.2)]
        self.assertAlmostEqual(
            stage4_curve_error(reference, with_bridge)["relative_error_percent"]["median"],
            10.0,
        )


if __name__ == "__main__":
    unittest.main()
