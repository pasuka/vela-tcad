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
from templates_ldmos_contracts import (  # noqa: E402
    SCHEMA_FILES,
    draft_governance_contracts,
    read_json,
    render_summary,
    validate_document,
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
        bv, states = build_bv("\t){ Coupled { Poisson Electron Hole Temperature } }\n")
        self.assertIn("state_bv_path", bv)
        self.assertEqual(states["path_samples"], 31)

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
                    {"name": "drain", "bias": 0.0},
                    {"name": "th_lat", "bias": 0.0},
                ],
            }), encoding="utf-8")
            separate_thermal_contacts(root, {"contacts": [
                {"name": "drain", "thermal_only": False},
                {"name": "th_lat", "thermal_only": True},
            ]})
            deck = json.loads((root / "simulation_iv.json").read_text(encoding="utf-8"))
            self.assertEqual([item["name"] for item in deck["contacts"]], ["drain"])
            thermal = json.loads((root / "thermal_contacts.json").read_text(encoding="utf-8"))
            self.assertFalse(thermal["electrical_use_authorized"])

    def test_converter_preserves_explicit_contact_edges(self) -> None:
        self.assertEqual(parse_edge_node_ids("1-2;2-3"), [[1, 2], [2, 3]])
        self.assertEqual(parse_edge_node_ids("1-2|2-3"), [[1, 2], [2, 3]])
        with self.assertRaisesRegex(ValueError, "invalid contact edge pair"):
            parse_edge_node_ids("1-2-3")


if __name__ == "__main__":
    unittest.main()
