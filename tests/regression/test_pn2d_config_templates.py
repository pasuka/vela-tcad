#!/usr/bin/env python3
"""Regression tests for versioned PN2D simulation templates."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path, PureWindowsPath
from typing import Any


REPO = Path(__file__).resolve().parents[2]
if str(REPO / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO / "scripts"))

from generate_pn2d_config import (  # noqa: E402
    BV_AVALANCHE_CURRENT_SUPPORT_PROFILES,
    BV_CONFIGURATION_PROFILES,
    TEMPLATES,
    TemplateError,
    render_named_template,
    validate_pn2d_config,
    write_rendered_config,
)


def all_strings(value: Any):
    if isinstance(value, dict):
        for item in value.values():
            yield from all_strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from all_strings(item)
    elif isinstance(value, str):
        yield value


class Pn2dConfigTemplatesTest(unittest.TestCase):
    def test_iv_defaults_capture_qualified_forward_models(self) -> None:
        config, manifest = render_named_template("pn2d_iv")
        self.assertEqual(config["sweep"]["mode"], "iv")
        self.assertEqual(config["sweep"]["initial_step"], 1.0e-3)
        self.assertEqual(config["sweep"]["min_step"], 1.0e-8)
        self.assertEqual(config["sweep"]["max_step"], 0.025)
        self.assertEqual(config["sweep"]["growth_factor"], 1.3)
        self.assertEqual(config["solver"]["mobility"]["model"], "masetti")
        self.assertEqual(
            config["solver"]["mobility"]["doping_concentration_basis"],
            "cell_reconstructed_total_impurity",
        )
        self.assertEqual(config["solver"]["impact_ionization"]["model"], "none")
        self.assertEqual(manifest["template_version"], 1)

    def test_bv_defaults_capture_sentaurus_aligned_models(self) -> None:
        config, manifest = render_named_template("pn2d_bv")
        self.assertEqual(config["sweep"]["mode"], "bv_reverse")
        self.assertEqual(config["sweep"]["initial_step"], 1.0e-4)
        self.assertEqual(config["sweep"]["min_step"], 1.0e-10)
        self.assertEqual(config["sweep"]["max_step"], 0.05)
        self.assertEqual(config["sweep"]["growth_factor"], 1.2)
        self.assertEqual(config["solver"]["max_iter"], 80)
        self.assertEqual(config["solver"]["mobility"]["model"], "masetti_field")
        self.assertEqual(
            config["solver"]["mobility"]["doping_concentration_basis"],
            "net_doping",
        )
        self.assertEqual(
            config["solver"]["impact_ionization"]["model"], "van_overstraeten"
        )
        impact = config["solver"]["impact_ionization"]
        for name, value in BV_AVALANCHE_CURRENT_SUPPORT_PROFILES[
            "element_edge_sg_gss_laux"
        ].items():
            self.assertEqual(impact[name], value)
        self.assertEqual(
            config["mesh_geometry"],
            {"node_volume_policy": "mixed_voronoi", "require_non_obtuse": True},
        )
        self.assertEqual(
            manifest["parameters"]["avalanche_current_support_profile"],
            "element_edge_sg_gss_laux",
        )
        self.assertEqual(manifest["template_version"], 3)
        self.assertEqual(manifest["overrides"], {})
        self.assertEqual(
            manifest["resolved_profile"],
            BV_CONFIGURATION_PROFILES["element_edge_sg_gss_laux"],
        )
        self.assertEqual(
            config["sweep"]["write_state_every_point_prefix"],
            "pn2d_bv_states/state",
        )
        newton_history = config["sweep"]["diagnostics"]["newton_history"]
        self.assertEqual(
            newton_history["attempts_csv_file"],
            "pn2d_bv_newton_attempts.csv",
        )
        self.assertEqual(
            newton_history["iterations_csv_file"],
            "pn2d_bv_newton_iterations.csv",
        )

    def test_bv_legacy_profile_is_an_atomic_rollback(self) -> None:
        config, manifest = render_named_template(
            "pn2d_bv",
            {"avalanche_current_support_profile": "legacy_cell_reconstructed"},
        )
        impact = config["solver"]["impact_ionization"]
        for name, value in BV_AVALANCHE_CURRENT_SUPPORT_PROFILES[
            "legacy_cell_reconstructed"
        ].items():
            self.assertEqual(impact[name], value)
        self.assertEqual(
            config["mesh_geometry"],
            {"node_volume_policy": "barycentric", "require_non_obtuse": False},
        )
        self.assertEqual(
            manifest["overrides"],
            {"avalanche_current_support_profile": "legacy_cell_reconstructed"},
        )

    def test_bv_sg_laux_default_profile_is_atomic(self) -> None:
        config, manifest = render_named_template(
            "pn2d_bv",
            {"avalanche_current_support_profile": "element_edge_sg_gss_laux"},
        )
        impact = config["solver"]["impact_ionization"]
        for name, value in BV_AVALANCHE_CURRENT_SUPPORT_PROFILES[
            "element_edge_sg_gss_laux"
        ].items():
            self.assertEqual(impact[name], value)
        self.assertEqual(
            config["mesh_geometry"],
            {"node_volume_policy": "mixed_voronoi", "require_non_obtuse": True},
        )
        self.assertEqual(
            manifest["overrides"],
            {"avalanche_current_support_profile": "element_edge_sg_gss_laux"},
        )

    def test_bv_legacy_config_remains_valid_but_mixed_profiles_fail_closed(self) -> None:
        legacy, _ = render_named_template(
            "pn2d_bv",
            {"avalanche_current_support_profile": "legacy_cell_reconstructed"},
        )
        validate_pn2d_config(legacy, "pn2d_bv")
        mixed = json.loads(json.dumps(legacy))
        mixed["solver"]["impact_ionization"]["current_approximation"] = (
            "element_edge_sg_gss_laux"
        )
        with self.assertRaisesRegex(TemplateError, "complete atomic profile"):
            validate_pn2d_config(mixed, "pn2d_bv")
        omitted = json.loads(json.dumps(legacy))
        del omitted["solver"]["impact_ionization"]["source_mapping_mode"]
        with self.assertRaisesRegex(TemplateError, "complete atomic profile"):
            validate_pn2d_config(omitted, "pn2d_bv")

        half_migrated = json.loads(json.dumps(legacy))
        half_migrated["mesh_geometry"] = {
            "node_volume_policy": "mixed_voronoi",
            "require_non_obtuse": True,
        }
        with self.assertRaisesRegex(TemplateError, "complete atomic profile"):
            validate_pn2d_config(half_migrated, "pn2d_bv")

        default, _ = render_named_template("pn2d_bv")
        default["mesh_geometry"]["node_volume_policy"] = "barycentric"
        with self.assertRaisesRegex(TemplateError, "complete atomic profile"):
            validate_pn2d_config(default, "pn2d_bv")

        wrong_type, _ = render_named_template("pn2d_bv")
        wrong_type["mesh_geometry"]["require_non_obtuse"] = 1
        with self.assertRaisesRegex(TemplateError, "must be boolean"):
            validate_pn2d_config(wrong_type, "pn2d_bv")

    def test_iv_template_does_not_inherit_bv_mesh_policy(self) -> None:
        config, _ = render_named_template("pn2d_iv")
        self.assertNotIn("mesh_geometry", config)

    def test_bvmethods_external_resistor_defaults_match_validated_method(self) -> None:
        config, manifest = render_named_template(
            "bvmethods_nmos_external_resistor"
        )
        sweep = config["sweep"]
        circuit = sweep["external_circuit"]
        self.assertEqual(sweep["mode"], "bv_reverse")
        self.assertEqual((sweep["start"], sweep["stop"], sweep["step"]),
                         (406.0, 1206.0, 200.0))
        self.assertEqual(circuit["mode"], "series_resistor")
        self.assertEqual(circuit["resistance_ohm_um"], 1.0e7)
        self.assertEqual(circuit["current_direction"], 1.0)
        self.assertEqual(circuit["initial_inner_voltage_V"], 5.9)
        self.assertEqual(circuit["max_inner_voltage_step_V"], 0.025)
        self.assertNotIn("voltage_to_current", sweep)
        self.assertEqual(
            config["solver"]["impact_ionization"]["coupling_mode"],
            "self_consistent",
        )
        self.assertEqual(manifest["template_version"], 1)

    def test_bvmethods_voltage_to_current_defaults_match_validated_method(self) -> None:
        config, manifest = render_named_template(
            "bvmethods_nmos_voltage_to_current"
        )
        sweep = config["sweep"]
        control = sweep["voltage_to_current"]
        self.assertEqual((sweep["start"], sweep["stop"], sweep["step"]),
                         (5.9, 6.0, 0.025))
        self.assertEqual(control["switch_voltage_V"], 6.0)
        self.assertEqual(control["current_direction"], 1.0)
        self.assertEqual(
            control["current_points_A_per_um"], [4.0e-5, 6.0e-5, 1.0e-4]
        )
        self.assertEqual(control["max_inner_voltage_step_V"], 0.0125)
        self.assertNotIn("external_circuit", sweep)
        self.assertEqual(manifest["template_version"], 1)

    def test_bvmethods_templates_reject_invalid_boundary_controls(self) -> None:
        external, _ = render_named_template("bvmethods_nmos_external_resistor")
        external["sweep"]["external_circuit"]["resistance_ohm_um"] = 0.0
        with self.assertRaisesRegex(TemplateError, "resistance must be positive"):
            validate_pn2d_config(external, "bvmethods_nmos_external_resistor")

        switched, _ = render_named_template("bvmethods_nmos_voltage_to_current")
        switched["sweep"]["voltage_to_current"]["switch_voltage_V"] = 5.95
        with self.assertRaisesRegex(TemplateError, "must equal"):
            validate_pn2d_config(switched, "bvmethods_nmos_voltage_to_current")

    def test_defaults_contain_no_absolute_paths_or_placeholders(self) -> None:
        for template in TEMPLATES:
            config, _ = render_named_template(template)
            for value in all_strings(config):
                self.assertNotIn("${", value)
                self.assertFalse(Path(value).is_absolute(), value)
                self.assertFalse(PureWindowsPath(value).is_absolute(), value)

    def test_unknown_and_invalid_overrides_are_rejected(self) -> None:
        with self.assertRaisesRegex(TemplateError, "unknown"):
            render_named_template("pn2d_iv", {"typo": 1})
        with self.assertRaisesRegex(TemplateError, "type number"):
            render_named_template("pn2d_iv", {"stop_voltage": "20"})
        with self.assertRaisesRegex(TemplateError, "must be relative"):
            render_named_template(
                "pn2d_iv", {"mesh_file": r"D:\external\mesh.json"}
            )
        with self.assertRaisesRegex(TemplateError, "positive sweep.step"):
            render_named_template("pn2d_bv", {"stop_voltage": 20.0})
        with self.assertRaisesRegex(TemplateError, "must be one of"):
            render_named_template(
                "pn2d_bv",
                {"avalanche_current_support_profile": "half_migrated"},
            )

    def test_rendering_is_byte_deterministic_and_manifest_is_separate(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_template_") as td:
            root = Path(td)
            first = root / "first.json"
            second = root / "second.json"
            write_rendered_config("pn2d_iv", first)
            write_rendered_config("pn2d_iv", second)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            first_manifest = json.loads(
                first.with_suffix(".manifest.json").read_text(encoding="utf-8")
            )
            rendered = json.loads(first.read_text(encoding="utf-8"))
            self.assertEqual(first_manifest["template"], "pn2d_iv")
            self.assertEqual(first_manifest["overrides"], {})
            self.assertNotIn("template_schema", rendered)

    def test_cli_writes_valid_bv_config_and_manifest(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_template_cli_") as td:
            output = Path(td) / "simulation.json"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "generate_pn2d_config.py"),
                    "--template",
                    "pn2d_bv",
                    "--output",
                    str(output),
                    "--set",
                    "write_vtk=true",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            config = json.loads(output.read_text(encoding="utf-8"))
            self.assertTrue(config["sweep"]["write_vtk"])
            self.assertTrue(output.with_suffix(".manifest.json").is_file())

    def test_json_schema_artifact_is_draft_2020_12(self) -> None:
        schema = json.loads(
            (REPO / "configs" / "schema" / "vela-simulation.schema.json")
            .read_text(encoding="utf-8")
        )
        self.assertEqual(
            schema["$schema"], "https://json-schema.org/draft/2020-12/schema"
        )
        self.assertIn("sweep", schema["required"])
        self.assertIn("initial_step", schema["properties"]["sweep"]["required"])
        sweep_properties = schema["properties"]["sweep"]["properties"]
        self.assertIn("write_state_every_point_prefix", sweep_properties)
        solver_properties = schema["properties"]["solver"]["properties"]
        self.assertIn("block_absolute_convergence", solver_properties)
        self.assertIn(
            "contact_majority_qf_branch_guard_contacts", solver_properties
        )
        newton_history = sweep_properties["diagnostics"]["properties"]["newton_history"]
        self.assertIn("attempts_csv_file", newton_history["properties"])
        self.assertIn("iterations_csv_file", newton_history["properties"])


if __name__ == "__main__":
    unittest.main()
