#!/usr/bin/env python3
"""Regression coverage for Templates/LDMOS WP0 and phase-0 contracts."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from run_templates_ldmos_sentaurus_vm import (  # noqa: E402
    BUNDLE_NAMES,
    REQUIRED_SOURCE_FILES,
    main as runner_main,
    parse_stages,
    prepare_bundle,
)
from analyze_templates_ldmos_structure import (  # noqa: E402
    discretization_draft,
    material_class,
    percentile,
    separate_thermal_contacts,
    triangle_angles,
    triangle_area,
)
from prepare_templates_ldmos_state_decks import (  # noqa: E402
    build_bv,
    build_idvd,
    build_idvg,
)
from audit_templates_ldmos_oracle import curve_metrics, strong_log_errors  # noqa: E402
from run_templates_ldmos_cost_probe import mesh_pattern, prepare_deck  # noqa: E402
from convert_tcad_export import parse_edge_node_ids  # noqa: E402
from classify_templates_ldmos_states import select_bv, terminal_values  # noqa: E402
from templates_ldmos_contracts import (  # noqa: E402
    SCHEMA_FILES,
    canonical_round_trip,
    draft_governance_contracts,
    migrate_discretization_draft,
    migrate_materials_v0,
    read_json,
    render_summary,
    validate_document,
)
from summarize_templates_ldmos_restart_qualification import (  # noqa: E402
    CURRENT_RESOLUTION_FLOOR_A_PER_UM,
    state_max_abs_deltas,
    summarize as summarize_restart_qualification,
    terminal_kcl,
)
from prepare_templates_ldmos_wp175_contracts import prepare as prepare_wp175  # noqa: E402
from prepare_templates_ldmos_restart_qualification import (  # noqa: E402
    prepare as prepare_restart_qualification,
)


def write_source(root: Path) -> None:
    root.mkdir(parents=True)
    (root / "sprocess_fps.cmd").write_text("struct tdr=n@node@\n", encoding="utf-8")
    common = (
        'File { Grid="@tdr@" Parameters="@parameter@" Output="@log@" '
        'Current="@plot@" Plot="@tdrdat@" }\nSave(FilePrefix="n@node@_state")\n'
    )
    for name in ("IdVg_des.cmd", "IdVd_des.cmd", "BVdss_des.cmd"):
        (root / name).write_text(common, encoding="utf-8")
    (root / "sdevice.par").write_text('Material="Silicon" {}\n', encoding="utf-8")
    (root / "gtree.dat").write_text("sprocess sprocess\n", encoding="utf-8")


class TemplatesLdmosContractsTest(unittest.TestCase):
    def test_all_contract_schemas_are_draft_2020_12_documents(self) -> None:
        for filename in SCHEMA_FILES.values():
            with self.subTest(filename=filename):
                schema = read_json(REPO / "schemas" / filename)
                self.assertEqual(schema["$schema"], "https://json-schema.org/draft/2020-12/schema")
                self.assertFalse(schema["additionalProperties"])

    def test_draft_governance_contracts_validate_and_reject_unknown_keys(self) -> None:
        contracts = draft_governance_contracts("270f525")
        self.assertEqual(
            set(contracts),
            {"known_difference_ledger.json", "threshold_freeze.json", "budget_freeze.json"},
        )
        for document in contracts.values():
            validate_document(document)
        invalid = json.loads(json.dumps(contracts["budget_freeze.json"]))
        invalid["unexpected"] = True
        with self.assertRaisesRegex(ValueError, "unexpected key"):
            validate_document(invalid)

    def test_tracked_budget_snapshot_is_double_approved_and_valid(self) -> None:
        budget = read_json(
            REPO / "reference_tcad" / "templates_ldmos_sentaurus2022" /
            "budget_freeze.json"
        )
        validate_document(budget)
        self.assertEqual(budget["approval"]["status"], "approved")
        self.assertTrue(budget["approval"]["benchmark_owner"])
        self.assertTrue(budget["approval"]["independent_reviewer"])
        self.assertTrue(budget["approval"]["approved_at"])

    def test_validation_summary_renderer_uses_machine_readable_gates(self) -> None:
        summary = {
            "schema": "vela.templates_ldmos.validation_summary.v1",
            "benchmark": "sentaurus_t2022_03_sp2_templates_ldmos",
            "generated_at": "2026-08-26T00:00:00Z",
            "status": "not_run",
            "highest_level": "none",
            "gates": [{"id": "oracle", "status": "not_run", "summary": "dry run"}],
            "evidence": [],
            "limitations": ["No live run."],
        }
        rendered = render_summary(summary)
        self.assertIn("| `oracle` | `not_run` | dry run |", rendered)
        self.assertIn("No live run.", rendered)

    def test_wp175_golden_contracts_validate_roundtrip_and_cross_reference(self) -> None:
        root = (
            REPO / "reference_tcad" / "templates_ldmos_sentaurus2022" / "contracts"
        )
        documents = {
            path.name: read_json(path)
            for path in root.glob("*.json")
        }
        self.assertEqual(
            set(documents),
            {"materials.json", "physics_contract.json", "discretization_contract.json"},
        )
        for document in documents.values():
            self.assertEqual(canonical_round_trip(document), document)
        physics = documents["physics_contract.json"]
        discretization = documents["discretization_contract.json"]
        self.assertEqual(physics["materials_file"], "materials.json")
        self.assertEqual(physics["revision"], 4)
        self.assertEqual(
            physics["bandgap_narrowing"]["sentaurus_dEg0_eV"], -0.01595)
        self.assertIn(
            "preserve physical contact quasi-Fermi potentials",
            physics["bandgap_narrowing"]["sentaurus_qf_reference_mapping"],
        )
        self.assertEqual(physics["mobility"]["g3_model"], "constant_field")
        self.assertFalse(physics["mobility"]["doping_dependence_enabled"])
        self.assertTrue(physics["mobility"]["high_field_saturation_enabled"])
        self.assertEqual(
            physics["discretization_profile"], discretization["profile_name"])
        self.assertFalse(physics["impact_ionization"]["enabled"])
        self.assertEqual(
            discretization["avalanche_profile"], "not_authorized_in_phase_a")

    def test_wp175_schema_rejects_unknown_keys_ranges_and_wrong_units(self) -> None:
        materials = read_json(
            REPO / "reference_tcad" / "templates_ldmos_sentaurus2022" /
            "contracts" / "materials.json"
        )
        invalid = json.loads(json.dumps(materials))
        invalid["materials"][0]["electron_mobility_cm2_per_V_s"] = 1.0e12
        with self.assertRaisesRegex(ValueError, "above maximum"):
            validate_document(invalid)
        invalid = json.loads(json.dumps(materials))
        invalid["materials"][0]["mobility"] = 1350.0
        with self.assertRaisesRegex(ValueError, "unexpected key"):
            validate_document(invalid)
        invalid = json.loads(json.dumps(materials))
        invalid["unit_system"]["mobility"] = "m^2/(V*s)"
        with self.assertRaisesRegex(ValueError, "expected constant"):
            validate_document(invalid)

    def test_wp175_version_migrations_require_explicit_units_and_are_lossless(self) -> None:
        legacy = {"materials": [{
            "name": "Si", "eps_r": 11.7, "ni": 1.0e16,
            "mun": 0.135, "mup": 0.048, "Nc_m3": 2.8e25,
            "Nv_m3": 1.04e25,
        }]}
        with self.assertRaisesRegex(ValueError, "source_unit_system"):
            migrate_materials_v0(legacy, "ambiguous")
        migrated = migrate_materials_v0(legacy, "legacy_si")
        silicon = migrated["materials"][0]
        self.assertEqual(silicon["intrinsic_carrier_density_cm3"], 1.0e10)
        self.assertEqual(silicon["electron_mobility_cm2_per_V_s"], 1350.0)
        self.assertEqual(silicon["conduction_band_density_of_states_cm3"], 2.8e19)
        self.assertEqual(canonical_round_trip(migrated), migrated)

        draft = {
            "schema": "vela.templates_ldmos.discretization_contract.v1-draft-unvalidated",
            "profile_name": "templates_ldmos_exact_topology_legacy_candidate",
            "applicable_mesh": "exact mesh", "current_support": "scharfetter_gummel_edge_flux",
            "control_volume": "barycentric", "field_recovery": "cell_reconstructed",
            "volume_source_mapping": "cell_reconstructed",
            "contact_edge_integration": "exact_imported_contact_boundary_edges",
            "obtuse_policy": "accepted_and_listed_with_barycentric_positive_control_volumes",
            "require_non_obtuse": False,
            "non_delaunay_policy": "accepted_only_for_structural_import; physics parity unresolved",
            "avalanche_profile": "legacy_cell_reconstructed",
            "forbidden_inference": "Do not apply PN2D bundle globally.",
            "physics_use_authorized": False,
            "status": "draft_pending_wp1_75_schema_and_independent_approval",
        }
        migrated_discretization = migrate_discretization_draft(draft)
        self.assertTrue(migrated_discretization["physics_use_authorized"])
        self.assertEqual(
            migrated_discretization["avalanche_profile"],
            "not_authorized_in_phase_a",
        )

    def test_wp175_preparer_validates_and_hashes_run_local_contracts(self) -> None:
        with tempfile.TemporaryDirectory(prefix="templates_ldmos_wp175_") as directory:
            stage1 = Path(directory)
            contracts = stage1 / "contracts"
            contracts.mkdir()
            draft = discretization_draft({"mesh_quality": {"obtuse_triangle_count": 4}})
            (contracts / "discretization_contract.json").write_text(
                json.dumps(draft), encoding="utf-8")
            report = prepare_wp175(stage1)
            self.assertEqual(report["status"], "pass")
            self.assertEqual(len(report["contracts"]), 3)
            output = contracts / "wp175"
            for item in report["contracts"]:
                self.assertTrue((output / item["path"]).is_file())
                self.assertEqual(len(item["sha256"]), 64)


class TemplatesLdmosRunnerTest(unittest.TestCase):
    def test_bundle_materialization_resolves_only_declared_workbench_tokens(self) -> None:
        with tempfile.TemporaryDirectory(prefix="templates_ldmos_bundle_") as directory:
            root = Path(directory)
            source, bundle = root / "source", root / "bundle"
            write_source(source)
            report = prepare_bundle(source, bundle)
            self.assertEqual(set(BUNDLE_NAMES.values()) | {"materialization.json"},
                             {path.name for path in bundle.iterdir()})
            self.assertEqual(len(report["files"]), len(BUNDLE_NAMES))
            for path in bundle.glob("*.cmd"):
                self.assertNotIn("@", path.read_text(encoding="utf-8"))
            self.assertIn("n1_fps.tdr", (bundle / "IdVg.cmd").read_text(encoding="utf-8"))
            self.assertIn("n4_des.plt", (bundle / "IdVd.cmd").read_text(encoding="utf-8"))
            self.assertIn("n6_des.tdr", (bundle / "BVdss.cmd").read_text(encoding="utf-8"))

    def test_stage_parser_preserves_dependency_order(self) -> None:
        self.assertEqual(parse_stages("bv,sprocess,idvg"), ["sprocess", "idvg", "bv"])
        with self.assertRaises(ValueError):
            parse_stages("unknown")

    def test_dry_run_writes_manifests_governance_contracts_and_reports(self) -> None:
        with tempfile.TemporaryDirectory(prefix="templates_ldmos_dry_") as directory:
            root = Path(directory)
            source = root / "source"
            write_source(source)
            exit_code = runner_main([
                "--source-dir", str(source),
                "--staging-root", str(root / "staging"),
                "--run-id", "dry_run",
                "--stages", "sprocess,idvg,idvd,bv",
            ])
            self.assertEqual(exit_code, 0)
            run_dir = root / "staging" / "dry_run"
            for name in (
                "source_manifest.json", "run_manifest.json", "artifact_manifest.json",
                "known_difference_ledger.json", "threshold_freeze.json", "budget_freeze.json",
            ):
                validate_document(read_json(run_dir / "manifest" / name))
            summary = read_json(run_dir / "reports" / "validation_summary.json")
            validate_document(summary)
            self.assertEqual(summary["status"], "not_run")
            self.assertTrue((run_dir / "reports" / "validation_summary.md").is_file())
            self.assertEqual(
                set(REQUIRED_SOURCE_FILES),
                {item["path"] for item in read_json(run_dir / "manifest" / "source_manifest.json")["records"]},
            )


class TemplatesLdmosStructureAuditTest(unittest.TestCase):
    def test_geometry_helpers_cover_right_and_obtuse_triangles(self) -> None:
        right = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
        obtuse = [(0.0, 0.0), (1.0, 0.0), (0.1, 0.01)]
        self.assertAlmostEqual(triangle_area(right), 0.5)
        self.assertAlmostEqual(sum(triangle_angles(right)), 180.0)
        self.assertGreater(max(triangle_angles(obtuse)), 90.0)
        self.assertEqual(percentile([0.0, 1.0, 2.0], 0.5), 1.0)

    def test_material_classification_and_discretization_are_fail_closed(self) -> None:
        self.assertEqual(material_class("Silicon"), "transport_semiconductor")
        self.assertEqual(material_class("PolySilicon"), "electrostatic_only_semiconductor_poly")
        self.assertEqual(material_class("Oxide"), "dielectric")
        draft = discretization_draft({"mesh_quality": {"obtuse_triangle_count": 4}})
        self.assertFalse(draft["physics_use_authorized"])
        self.assertFalse(draft["require_non_obtuse"])
        self.assertEqual(draft["avalanche_profile"], "legacy_cell_reconstructed")
        self.assertIn("PN2D", draft["forbidden_inference"])

    def test_state_decks_add_plot_only_controls_at_declared_points(self) -> None:
        idvg, vg_states = build_idvg(
            "\tCoupled { Poisson Electron Hole }\n"
            "\t\tCurrentPlot( Time= (Range=(0 1) Intervals= 30)  )\n",
            1.5,
        )
        self.assertIn("state_idvg_eq_0V", idvg)
        self.assertEqual(vg_states["biases_V"], [0.0, 1.5, 2.5, 5.0])
        idvd, vd_states = build_idvd(
            "\t\tCurrentPlot( Time= (Range= (0 1) Intervals= 30) )\n"
            "\t\tCurrentPlot( Time= (Range= (0 1) Intervals= 30)  )\n"
        )
        self.assertIn("state_idvd_vg4", idvd)
        self.assertIn("state_idvd_vg8", idvd)
        self.assertEqual(vd_states["drain_biases_V"][-1], 40.0)
        targets = [
            {"role": "pre_iadapt", "time": 1.0, "voltage_V": 49.0,
             "current_A_per_um": 5e-13},
            {"role": "criterion_post", "time": 2.0, "voltage_V": 50.0,
             "current_A_per_um": 1.1e-8},
        ]
        bv, states = build_bv(
            "\t){ Coupled { Poisson Electron Hole Temperature } }\n", targets,
        )
        self.assertIn("state_bv_path", bv)
        self.assertEqual(states["requested_times"], [1.0, 2.0])
        self.assertEqual(states["expected_state_count"], 2)

    def test_oracle_curve_and_log_audit_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory(prefix="templates_ldmos_audit_") as directory:
            root = Path(directory)
            curve = root / "curve.csv"
            curve.write_text(
                "bias_V,current_total_A_per_um\n0,1e-12\n1,1e-9\n2,1e-6\n",
                encoding="utf-8",
            )
            metrics = curve_metrics(curve)
            self.assertTrue(metrics["finite"])
            self.assertTrue(metrics["bias_nondecreasing"])
            log = root / "run.out"
            log.write_text(
                "No errors reported\nError smaller than 1 (1e-4).\n"
                "Fatal: controlled failure\n",
                encoding="utf-8",
            )
            errors = strong_log_errors(log)
            self.assertEqual(len(errors), 1)
            self.assertEqual(errors[0]["line"], 3)

    def test_cost_probe_matrix_pattern_is_exact_for_tri3_graph(self) -> None:
        mesh = {
            "nodes": [{"id": value} for value in range(4)],
            "triangles": [
                {"node_ids": [0, 1, 2]},
                {"node_ids": [1, 3, 2]},
            ],
            "contacts": [{"node_ids": [0, 1]}],
        }
        pattern = mesh_pattern(mesh)
        self.assertEqual(pattern["unique_edges"], 5)
        self.assertEqual(pattern["poisson_structural_nnz_upper_bound"], 14)
        self.assertEqual(pattern["coupled_dd_unknowns_estimate"], 12)

    def test_cost_probe_uses_true_linear_poisson_not_dc_initial_coupled(self) -> None:
        with tempfile.TemporaryDirectory(prefix="templates_ldmos_cost_") as directory:
            root = Path(directory)
            deck_path = prepare_deck({
                "simulation_type": "dc_sweep",
                "mesh_file": "old.json",
                "node_doping_file": "doping.csv",
                "doping": [{"region": "Si", "donors": 1e20, "acceptors": 1e15}],
                "contacts": [{"name": "drain", "bias": 0.0}],
                "solver": {"method": "poisson_only"},
                "sweep": {"mode": "iv"},
            }, root, "probe", 0.4)
            deck = json.loads(deck_path.read_text(encoding="utf-8"))
            self.assertEqual(deck["simulation_type"], "poisson")
            self.assertNotIn("sweep", deck)
            self.assertNotIn("solver", deck)
            self.assertNotIn("node_doping_file", deck)
            self.assertEqual(deck["doping"][0]["donors"], 0.0)
            self.assertEqual(deck["contacts"][0]["bias"], 0.4)

    def test_thermal_contact_is_not_written_as_an_electrical_bias(self) -> None:
        with tempfile.TemporaryDirectory(prefix="templates_ldmos_thermal_") as directory:
            root = Path(directory)
            (root / "simulation_iv.json").write_text(json.dumps({
                "_comment": "generated",
                "contacts": [
                    {"name": "gate", "bias": 0.0},
                    {"name": "drain", "bias": 0.0},
                    {"name": "th_lat", "bias": 0.0},
                ],
            }), encoding="utf-8")
            separate_thermal_contacts(root, {"contacts": [
                {
                    "name": "gate", "owner_region": "Oxide_1.1",
                    "thermal_only": False,
                },
                {"name": "drain", "thermal_only": False},
                {"name": "th_lat", "thermal_only": True},
            ]})
            deck = json.loads((root / "simulation_iv.json").read_text(encoding="utf-8"))
            self.assertEqual([item["name"] for item in deck["contacts"]], ["gate", "drain"])
            self.assertEqual(deck["contacts"][0]["type"], "metal_gate")
            self.assertEqual(deck["contacts"][0]["flatband_voltage"], 0.0)
            thermal = json.loads((root / "thermal_contacts.json").read_text(encoding="utf-8"))
            self.assertFalse(thermal["electrical_use_authorized"])

    def test_converter_preserves_explicit_contact_edges(self) -> None:
        self.assertEqual(parse_edge_node_ids("1-2;2-3"), [[1, 2], [2, 3]])
        self.assertEqual(parse_edge_node_ids("1-2|2-3"), [[1, 2], [2, 3]])
        with self.assertRaisesRegex(ValueError, "invalid contact edge pair"):
            parse_edge_node_ids("1-2-3")

    def test_representative_state_classifier_uses_exact_saved_bv_states(self) -> None:
        document = {
            "geometry": {"regions": [{"index": 6, "name": "drain"}]},
            "fields": [
                {"name": "ContactExternalVoltage", "region": 6, "raw_values": [12.0]},
                {"name": "ContactCurrentFlux", "region": 6, "raw_values": [1e-10]},
            ],
        }
        self.assertEqual(terminal_values(document)["drain"]["voltage_V"], 12.0)
        records = [
            {"name": f"state_{index}.tdr", "terminals": {"drain": {
                "voltage_V": float(index), "current_A_per_um": current,
            }}}
            for index, current in enumerate((1e-14, 5e-13, 7e-13, 1e-10, 9e-9, 1.1e-8))
        ]
        selected = select_bv(records, 6.5e-13, 1e-8)
        self.assertEqual(selected["pre_iadapt"]["state"], "state_1.tdr")
        self.assertEqual(selected["near_iadapt"]["state"], "state_2.tdr")
        self.assertEqual(selected["criterion_pre"]["state"], "state_4.tdr")
        self.assertEqual(selected["criterion_post"]["state"], "state_5.tdr")


class TemplatesLdmosRestartQualificationTest(unittest.TestCase):
    @staticmethod
    def _write_state(path: Path, psi: float) -> None:
        path.write_text(
            "node_id,psi,phin,phip,electrons_m3,holes_m3\n"
            f"0,{psi:.17g},0,0,1e20,2e20\n"
            f"1,{(psi + 0.1):.17g},0.01,-0.01,2e20,1e20\n",
            encoding="utf-8",
        )

    @staticmethod
    def _write_sweep(path: Path, current: float) -> None:
        path.write_text(
            "converged,iterations,newton_convergence_reason,"
            "final_psi_residual_norm,current_total_A_per_um\n"
            f"1,4,reltol,1e-10,{current:.17g}\n",
            encoding="utf-8",
        )

    @staticmethod
    def _write_terminal(path: Path, current: float, imbalance: float) -> None:
        path.write_text(
            "contact,current_total_A_per_um\n"
            f"source,{-current + imbalance:.17g}\n"
            f"drain,{current:.17g}\n"
            "gate,0\nsubstrate,0\n",
            encoding="utf-8",
        )

    def _write_fixture(self, root: Path) -> None:
        prepared = prepare_restart_qualification(root.parent)
        self.assertEqual(Path(next(iter(prepared.values()))).parent, root)
        states = {
            "sentaurus_eq_0v_state.csv": 0.0,
            "frozen_eq_0v_roundtrip.csv": 0.0,
            "reclose_eq_0v_state.csv": 0.2,
            "reclose_eq_0v_repeat_state.csv": 0.2 + 1.0e-6,
            "reclose_idvg_vg0_vd0p1_state.csv": 0.3,
            "reclose_idvg_vg0_vd0p1_repeat_state.csv": 0.3 + 100.0e-6,
            "reclose_idvg_vg0_vd0p1_repeat2_state.csv": 0.3 + 100.0e-6 + 1.0e-9,
            "reclose_idvd_vg4_vd0p1_state.csv": 0.4,
            "reclose_idvd_vg4_vd0p1_repeat_state.csv": 0.4 + 1.0e-9,
        }
        for name, psi in states.items():
            self._write_state(root / name, psi)
        sweeps = {
            "reclose_eq_0v.csv": 0.0,
            "reclose_eq_0v_repeat.csv": 0.0,
            "reclose_idvg_vg0_vd0p1_repeat.csv": 1.0e-19,
            "reclose_idvg_vg0_vd0p1_repeat2.csv": 1.00000001e-19,
            "reclose_idvd_vg4_vd0p1.csv": 1.0e-5,
            "reclose_idvd_vg4_vd0p1_repeat.csv": 1.00000001e-5,
        }
        for name, current in sweeps.items():
            self._write_sweep(root / name, current)
        self._write_terminal(
            root / "reclose_idvd_vg4_vd0p1_terminal_balance.csv", 1.0e-5, 1.0e-12)
        self._write_terminal(
            root / "reclose_idvd_vg4_vd0p1_repeat_terminal_balance.csv",
            1.00000001e-5, 1.0e-13)

    def test_state_and_terminal_metrics_are_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory(prefix="templates_ldmos_restart_metric_") as directory:
            root = Path(directory)
            self._write_state(root / "first.csv", 0.0)
            self._write_state(root / "second.csv", 2.0e-6)
            self.assertAlmostEqual(
                state_max_abs_deltas(root / "first.csv", root / "second.csv")["psi"],
                2.0e-6,
            )
            self._write_terminal(root / "terminal.csv", 1.0e-5, 1.0e-12)
            self.assertAlmostEqual(
                terminal_kcl(root / "terminal.csv")["normalized_imbalance"],
                1.0e-7,
                places=10,
            )

    def test_restart_qualification_distinguishes_settled_and_unresolved_states(self) -> None:
        with tempfile.TemporaryDirectory(prefix="templates_ldmos_restart_summary_") as directory:
            root = Path(directory) / "qualification"
            self._write_fixture(root)
            summary = summarize_restart_qualification(root)
            self.assertEqual(summary["status"], "pass")
            idvg = summary["cases"]["idvg_low_drain_reclose"]
            self.assertGreater(
                idvg["max_abs_imported_to_settled_delta"]["psi"], 10.0e-6)
            self.assertLess(
                idvg["max_abs_settled_reclose_delta"]["psi"], 10.0e-6)
            self.assertFalse(idvg["current_resolved"])
            self.assertLess(
                abs(idvg["sweep_points"][0]["current_total_A_per_um"]),
                CURRENT_RESOLUTION_FLOOR_A_PER_UM,
            )
            self._write_state(
                root / "reclose_idvd_vg4_vd0p1_repeat_state.csv", 0.401)
            failed = summarize_restart_qualification(root)
            self.assertEqual(failed["status"], "fail")
            self.assertEqual(
                next(gate for gate in failed["gates"]
                     if gate["id"] == "idvd_prebias_same_bias_psi")["status"],
                "fail",
            )

    def test_restart_deck_generator_freezes_dependencies_and_solver_controls(self) -> None:
        with tempfile.TemporaryDirectory(prefix="templates_ldmos_restart_decks_") as directory:
            stage1 = Path(directory)
            written = prepare_restart_qualification(stage1)
            self.assertEqual(len(written), 8)
            qualification = stage1 / "qualification"
            frozen = json.loads((qualification / "frozen_eq_0v.json").read_text())
            self.assertEqual(frozen["solver"]["method"], "frozen_state")
            self.assertFalse(frozen["sweep"]["frozen_state_compute_current"])
            repeat2 = json.loads(
                (qualification / "reclose_idvg_vg0_vd0p1_repeat2.json").read_text())
            self.assertEqual(
                repeat2["sweep"]["initial_state_file"],
                "reclose_idvg_vg0_vd0p1_repeat_state.csv",
            )
            self.assertEqual(repeat2["contacts"][0]["type"], "metal_gate")
            self.assertEqual(repeat2["contacts"][0]["flatband_voltage"], 0.0)
            self.assertEqual(
                repeat2["solver"]["bandgap_narrowing"]["model"], "old_slotboom")
            self.assertTrue(repeat2["solver"]["continuity_row_scaling"]["enabled"])


if __name__ == "__main__":
    unittest.main()
