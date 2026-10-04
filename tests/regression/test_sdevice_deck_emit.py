"""WP4 contracts, physical mappings and real-runner smoke checks."""
from __future__ import annotations
from copy import deepcopy
import csv
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))
from vela_sdevice.cmd_parser import parse_cmd_text
from vela_sdevice.deck_emit import (DeckEmitError, EmitOptions, TDR_FIELDS, emit_stage_plan,
                                    expand_times, rebase_stage, v2_key, write_stage_plan)
from vela_sdevice.materials import (active_models_from_ir, default_models_par_path,
                                   emit_materials, write_materials)
from vela_sdevice import schema_validator as sv

FIXTURES = REPO / "tests/fixtures/sdevice"
SUPPORTED = [f for f in json.loads((FIXTURES / "fixtures.json").read_text())["fixtures"]
             if f["compare"] == "full"]


def accepted(relative: str) -> dict:
    path = FIXTURES / relative
    text = path.read_text(encoding="utf-8")
    if not re.search(r"Parameter\s*=", text):
        text = re.sub(r"\bFile\s*\{", 'File { Parameter="models.par" ', text, count=1)
    ir = parse_cmd_text(text, path_name=path.name)
    assert ir["verdict"]["accepted"], ir["errors"]
    return ir


def qualification(ir: dict, materials: tuple = ("Si",)) -> dict:
    return {"schema": "vela.mesh_qualification.v1", "source": {"grid_file": "mesh.grd",
        "format": "dfise", "sha256": "0" * 64}, "dimension": 2,
        "coordinate_unit": {"declared": "um", "applied_scale_to_um": 1.0},
        "bounding_box_um": {"min": [0.0, 0.0], "max": [1.0, 1.0]},
        "nodes": {"count": 4}, "elements": {"count": 2, "by_kind": {"triangle": 2},
            "converted_quadrilaterals": 0},
        "regions": [{"name": "bulk" if m == "Si" else "oxide", "material": m,
                     "source_material": "Silicon" if m == "Si" else "Oxide", "element_count": 2} for m in materials],
        "materials": list(materials),
        "contacts": [{"name": e["name"], "node_count": 2, "edge_count": 1,
                      "adjacent_materials": ["Si"]} for e in ir["electrodes"]],
        "angles": {"max_deg": 90.0, "min_deg": 45.0, "obtuse_element_count": 0, "obtuse_fraction": 0.0},
        "doping": {"file": "doping.csv", "columns": ["node_id", "donors_cm3", "acceptors_cm3"],
            "source_fields": ["DonorConcentration", "AcceptorConcentration"]},
        "outputs": {"mesh_json": "mesh.json", "doping_csv": "doping.csv"},
        "qualification": {"two_dimensional": True, "coordinates_plausible": True,
            "triangles_only": True, "non_obtuse": True,
            "recommended_node_volume_policy": "barycentric", "allow_obtuse": True}}


def material_result(ir: dict, materials: tuple = ("Si",)):
    return emit_materials(default_models_par_path(), active_models=active_models_from_ir(ir),
        temperature_K=ir["physics"]["temperature_K"] or 300.0, mesh_materials=materials,
        allow_lossy=True, effective_models=ir["physics"]["effective"])


def emit(ir: dict, **kwargs):
    return emit_stage_plan(ir, mesh_qualification=kwargs.pop("mesh_qualification", qualification(ir)),
        materials_report=kwargs.pop("materials_report", material_result(ir).report),
        work_dir="device_vela", cmd_basename=ir["source"]["path"],
        cmd_sha256=ir["source"]["sha256"], options=kwargs.pop("options", EmitOptions(allow_lossy=True)), **kwargs)


class DeckEmitterTests(unittest.TestCase):
    def setUp(self):
        self.ir = accepted("pn2d_sentaurus2018/pn2d_iv_sdevice.cmd")

    def assert_error(self, identifier, function, *args, **kwargs):
        with self.assertRaises(DeckEmitError) as caught:
            function(*args, **kwargs)
        self.assertEqual(caught.exception.message_id, identifier)

    def test_accepted_fixture_goldens_and_full_schema_validation(self):
        try:
            import jsonschema
        except ImportError:
            jsonschema = None
        schema = json.loads((REPO / "configs/schema/vela-simulation.schema.json").read_text())
        for fixture in SUPPORTED:
            with self.subTest(deck=fixture["deck"]):
                ir = accepted(fixture["deck"])
                materials = ("Si", "SiO2") if fixture["case"] == "bvmethods_sentaurus2018" else ("Si",)
                result = emit(ir, mesh_qualification=qualification(ir, materials),
                              materials_report=material_result(ir, materials).report)
                folder = FIXTURES / Path(fixture["deck"]).with_suffix("")
                self.assertEqual(result.stage_plan, json.loads((folder / "expected_stage_plan.json").read_text()))
                if jsonschema is not None:
                    jsonschema.Draft202012Validator(sv.load_schema(REPO / "schemas/vela.stage_plan.v1.schema.json")).validate(result.stage_plan)
                for i, deck in enumerate(result.decks):
                    self.assertEqual(deck, json.loads((folder / f"expected_stage{i:02d}.json").read_text()))
                    sv.validate(schema, deck)
                    if jsonschema is not None:
                        jsonschema.Draft202012Validator(schema).validate(deck)

    def test_time_voltage_conversion_and_reverse_sweep(self):
        ir = deepcopy(self.ir)
        ir["stages"][-1]["outputs"]["plot_files"][0]["time"]["intervals"] = 400
        result = emit(ir)
        sweep = result.decks[-1]["sweep"]
        self.assertEqual(sweep["max_step"], 0.25)
        self.assertEqual(sweep["initial_step"], 0.01)
        self.assertEqual(sweep["shrink_factor"], 0.5)
        self.assertEqual(len(sweep["vtk_bias_points"]), 401)
        self.assertEqual(sweep["vtk_bias_points"][0], 0)
        self.assertEqual(sweep["vtk_bias_points"][-1], 10)
        ir["stages"][-1]["goal"]["voltage_V"] = -10
        reverse = emit(ir).decks[-1]["sweep"]
        self.assertEqual(reverse["step"], -0.25)
        self.assertEqual(reverse["vtk_bias_points"][-1], -10)
        self.assertEqual(reverse["max_step"], 0.25)

    def test_initial_stage_merge_and_state_chain(self):
        result = emit(self.ir)
        self.assertEqual(len(result.decks), 2)
        stages = result.stage_plan["stages"]
        self.assertEqual(stages[0]["kind"], "equilibrium")
        self.assertIsNone(stages[0]["input_state"])
        self.assertEqual(result.decks[0]["sweep"]["initialization"], {"mode": "poisson_block"})
        self.assertEqual(stages[1]["input_state"], stages[0]["output_state"])
        self.assertNotIn("initialization", result.decks[1]["sweep"])

    def test_multistage_rebase_updates_later_contacts_and_step_sizes(self):
        ir = accepted("genius_bjt_sentaurus2022/bjt_m0_des.cmd")
        result = emit(ir)
        stages = result.stage_plan["stages"]
        self.assertEqual(stages[2]["electrode_voltages_start_V"], stages[1]["electrode_voltages_goal_V"])
        rebased = rebase_stage(result, 2, {"collector": 3.2, "base": .65})
        self.assertEqual(rebased.stage_plan["stages"][2]["sweep"]["start_V"], 3.2)
        self.assertAlmostEqual(rebased.decks[2]["sweep"]["max_step"], .01)
        self.assertLess(rebased.decks[2]["sweep"]["step"], 0)
        self.assertEqual(rebased.decks[2]["contacts"][1]["bias"], .65)
        self.assertEqual(result.stage_plan["stages"][2]["sweep"]["start_V"], 0)

    def test_global_plot_indices_and_currentplot_required_union(self):
        ir = deepcopy(self.ir)
        second = deepcopy(ir["stages"][-1])
        second.update(index=3, depends_on=[2])
        second["goal"]["voltage_V"] = 20
        second["outputs"]["current_plot_times"] = {"range": [0, 1], "intervals": 400}
        ir["stages"].append(second)
        result = emit(ir)
        last = result.stage_plan["stages"][-1]
        self.assertEqual(last["vtu_write"]["series"][0]["global_index_start"], 201)
        self.assertEqual(len(last["plt_write"]), 401)
        self.assertEqual(len(last["required_t"]), 401)
        self.assertEqual(result.decks[-1]["sweep"]["step"], .025)

    def test_zero_span_and_time_errors(self):
        ir = deepcopy(self.ir)
        ir["stages"][-1]["goal"]["voltage_V"] = 0
        result = emit(ir)
        self.assertEqual(result.decks[-1]["sweep"]["bias_points"], [0])
        self.assertTrue(any(d.get("message_id") == "W211" for d in result.decisions))
        self.assert_error("E407", expand_times, {"list": [1.01]})
        self.assert_error("E407", expand_times, {"range": [1, 0], "intervals": 2})
        self.assert_error("E407", expand_times, {"list": [float("nan")]})

    def test_sidecar_whitelist_deep_merge_and_stage_selection(self):
        sidecar = {"schema": "vela.sdevice_sidecar.v1", "apply_to": [1],
                   "solver": {"max_iter": 120, "normalization": {"reference_concentration_cm3": 1e16}},
                   "sweep": {"diagnostics": {"newton_history": {"enabled": True}}}}
        result = emit(self.ir, sidecar=sidecar)
        self.assertEqual(result.decks[0]["solver"]["max_iter"], 100)
        self.assertEqual(result.decks[1]["solver"]["max_iter"], 120)
        self.assertIn("solver.max_iter", result.stage_plan["stages"][1]["sidecar_overrides"])
        for bad in ({"contacts": []}, {"solver": {"mobility": {"model": "constant"}}},
                    {"sweep": {"diagnostics": {}, "stop": 7}}, {"solver": {"tau_n_s": 1e-7}},
                    {"solver": {"max_iter": 0}}, {"apply_to": [9]}):
            self.assert_error("E409", emit, self.ir, sidecar={"schema": "vela.sdevice_sidecar.v1", **bad})

    def test_mixed_voronoi_accepts_right_angles_and_rejects_obtuse_mesh(self):
        sidecar = {"schema": "vela.sdevice_sidecar.v1", "mesh_geometry": {"node_volume_policy": "mixed_voronoi"}}
        mesh = qualification(self.ir)
        result = emit(self.ir, sidecar=sidecar, mesh_qualification=mesh)
        self.assertTrue(result.decks[0]["mesh_geometry"]["require_non_obtuse"])
        mesh["angles"]["max_deg"] = 90.01
        mesh["angles"]["obtuse_element_count"] = 1
        mesh["qualification"]["non_obtuse"] = False
        self.assert_error("E408", emit, self.ir, sidecar=sidecar, mesh_qualification=mesh)

    def test_explicit_avalanche_profile_is_atomic_and_preserves_physics(self):
        from generate_pn2d_config import BV_CONFIGURATION_PROFILES
        ir = accepted("pn2d_sentaurus2018/pn2d_bv_sdevice.cmd")
        baseline = emit(ir)
        for name, bundle in BV_CONFIGURATION_PROFILES.items():
            sidecar = {"schema": "vela.sdevice_sidecar.v1", "avalanche_current_support_profile": name}
            result = emit(ir, sidecar=sidecar)
            for before, after in zip(baseline.decks, result.decks):
                self.assertEqual(after["mesh_geometry"], bundle["mesh_geometry"])
                for key, value in bundle["impact_ionization"].items():
                    self.assertEqual(after["solver"]["impact_ionization"][key], value)
                self.assertEqual(before["solver"]["impact_ionization"]["electron_a_low_cm_inv"],
                                 after["solver"]["impact_ionization"]["electron_a_low_cm_inv"])
                self.assertEqual(before["solver"]["mobility"], after["solver"]["mobility"])
                self.assertEqual(before["sweep"]["stop"], after["sweep"]["stop"])
            self.assert_error("E409", emit, ir, sidecar={**sidecar, "mesh_geometry": {
                "node_volume_policy": "barycentric" if name == "element_edge_sg_gss_laux" else "mixed_voronoi"}})
        self.assert_error("E409", emit, self.ir, sidecar={"schema": "vela.sdevice_sidecar.v1",
                          "avalanche_current_support_profile": "element_edge_sg_gss_laux"})

    def test_gate_barrier_workfunction_and_material(self):
        ir = deepcopy(self.ir)
        ir["electrodes"][0]["attributes"] = {"barrier_V": .55}
        self.assertEqual(emit(ir).decks[0]["contacts"][0]["work_function_eV"], .55)
        # Independent 300 K reference from the current WP3 default Si library.
        reference_eV = 4.611003610676079
        ir["electrodes"][0]["attributes"] = {"workfunction_eV": reference_eV}
        self.assertAlmostEqual(emit(ir).decks[0]["contacts"][0]["work_function_eV"], 0, delta=1e-9)
        ir["electrodes"][0]["attributes"] = {"workfunction_eV": 4.6228018936}
        self.assertAlmostEqual(emit(ir).decks[0]["contacts"][0]["work_function_eV"],
                               0.011798282923921377, delta=1e-9)
        report = material_result(ir).report
        report["electrodes"].append({"material": "Gold", "work_function_eV": 5.1})
        ir["electrodes"][0]["attributes"] = {"material": "Gold"}
        self.assertAlmostEqual(emit(ir, materials_report=report).decks[0]["contacts"][0]["work_function_eV"], 5.1 - reference_eV, delta=1e-9)
        self.assert_error("E407", emit, ir)

    def test_schottky_uses_absolute_workfunction_and_newton(self):
        ir = accepted("schottky_charon_sentaurus2018/schottky_n_des.cmd")
        result = emit(ir)
        attr = next(e["attributes"] for e in ir["electrodes"] if e["attributes"].get("schottky"))
        contact = next(c for c in result.decks[0]["contacts"] if c["type"] == "schottky")
        self.assertEqual(contact["work_function_eV"], attr["workfunction_eV"])
        self.assertEqual(result.decks[0]["solver"]["method"], "newton")
        self.assert_error("E407", emit, ir, sidecar={"schema": "vela.sdevice_sidecar.v1", "solver": {"method": "gummel"}})

    def test_material_units_and_math_heuristics(self):
        ir = accepted("pn2d_sentaurus2018/pn2d_bv_sdevice.cmd")
        report = material_result(ir).report
        result = emit(ir, materials_report=report)
        for row in report["solver_parameters"]:
            if row["section"] == "Auger" and row["parameter"] != "A+B(T/300)+C(T/300)^2":
                continue
            node = result.decks[-1]
            for key in v2_key(row["target"]).split("."):
                node = node[key]
            expected = row["value_si"] if row["target"].startswith("solver.band_to_band.") else row["value_par"]
            self.assertAlmostEqual(node / expected, 1) if expected else self.assertEqual(node, 0)
        self.assertNotIn("reltol", result.decks[0]["solver"])
        self.assertTrue(any(d["keyword"] == "digits" for d in result.decisions))
        self.assertEqual(v2_key("solver.mobility.electron_mumin1_m2_V_s"), "solver.mobility.electron_mumin1_cm2_V_s")
        self.assertEqual(v2_key("solver.band_to_band.a_m_inv"), "solver.band_to_band.a_m_inv")

    def test_rejected_and_unsupported_ir_are_fail_closed(self):
        ir = deepcopy(self.ir)
        ir["verdict"] = {"accepted": False, "exit_code": 1, "fatal_ids": ["E106"]}
        self.assert_error("E106", emit, ir)
        ir = deepcopy(self.ir)
        ir["physics"]["effective"][0].update(status="unsupported", finding_id="E106", vela=None)
        self.assert_error("E106", emit, ir)

    def test_input_purity_and_unknown_electrode(self):
        ir = deepcopy(self.ir)
        mesh = qualification(ir)
        report = material_result(ir).report
        before = deepcopy((ir, mesh, report))
        emit(ir, mesh_qualification=mesh, materials_report=report)
        self.assertEqual((ir, mesh, report), before)
        mesh["contacts"] = mesh["contacts"][1:]
        self.assert_error("E401", emit, ir, mesh_qualification=mesh, materials_report=report)

    def test_options_default_controls_and_explicit_plot_fields(self):
        ir = deepcopy(self.ir)
        ir["plot"]["explicit"] = True
        ir["plot"]["fields"] = [{"name": "ElectricField", "data_name": "ElectricField", "component": "vector", "line": 1}]
        ir["stages"][-1]["t_control"] = dict.fromkeys(
            ("initial_step", "min_step", "max_step", "increment", "decrement")) | {"do_zero": True}
        result = emit(ir, options=EmitOptions(threads=3, exit_on_failure=True,
                                              plot_final_only=True, vtu_binary=True))
        self.assertEqual(result.stage_plan["threads"], 3)
        self.assertTrue(result.stage_plan["exit_on_failure"])
        self.assertFalse(result.stage_plan["stages"][1]["vtu_write"]["per_accepted_step"])
        sweep = result.decks[1]["sweep"]
        self.assertEqual(sweep["vtk_fields"], ["ElectricField/Vector"])
        self.assertTrue(sweep["vtk_binary"])
        self.assertEqual(sweep["initial_step"], 1)
        self.assertEqual(sweep["min_step"], .01)
        self.assertEqual(sweep["max_step"], 10)

    def test_material_lifetimes_auger_and_si_band_to_band(self):
        ir = deepcopy(self.ir)
        report = material_result(ir).report
        row = {"material": "Si", "section": "Auger", "parameter": "A+B(T/300)+C(T/300)^2",
               "target": "solver.auger_cn_m6_per_s", "value_si": 2.9e-43,
               "value_par": 2.9e-31, "status": "frozen_at_temperature"}
        report["solver_parameters"].append(row)
        report["solver_parameters"].append({**row, "section": "Band2Band", "parameter": "A",
            "target": "solver.band_to_band.a_e2", "value_si": 1.9e15, "value_par": 1.9e9})
        solver = emit(ir, materials_report=report).decks[0]["solver"]
        self.assertAlmostEqual(solver["auger_cn_cm6_per_s"] / 2.9e-31, 1)
        self.assertEqual(solver["band_to_band"]["a_e2"], 1.9e15)
        self.assertEqual(solver["tau_n_s"], solver["srh_doping_dependence"]["electron"]["tau_max_s"])
        self.assertEqual(solver["tau_p_s"], solver["srh_doping_dependence"]["hole"]["tau_max_s"])

    def test_srh_parameter_declarations_do_not_enable_doping_dependence(self):
        for option, enabled in (("", False), ("(DopingDependence)", True)):
            with self.subTest(option=option):
                ir = parse_cmd_text('''File { Grid="mesh.grd" Parameter="models.par" }
                    Electrode { {Name="left" Voltage=0} {Name="right" Voltage=0} }
                    Physics { Recombination(SRH''' + option + ''') }
                    Solve { Coupled {Poisson Electron Hole} }''', path_name="srh.cmd")
                solver = emit(ir).decks[0]["solver"]
                self.assertFalse(solver["contact_boundary_minority_electron_relaxation"])
                self.assertEqual(solver["srh_doping_dependence"]["enabled"], enabled)
                self.assertEqual(solver["tau_n_s"], 1e-5)
                self.assertEqual(solver["tau_p_s"], 3e-6)

    def test_drift_diffusion_driving_force_defaults_and_explicit_override(self):
        for force, expected in (("", "quasi_fermi_gradient"), ("(ElectricField)", "electric_field")):
            ir = parse_cmd_text('''File { Grid="mesh.grd" Parameter="models.par" }
                Electrode { {Name="left" Voltage=0} {Name="right" Voltage=0} }
                Physics { Mobility(DopingDependence HighFieldSaturation''' + force + ''')
                  Recombination(Avalanche(VanOverstraeten ''' + ("ElectricField" if force else "") + ''')) }
                Solve { Coupled {Poisson Electron Hole} }''', path_name="force.cmd")
            solver = emit(ir).decks[0]["solver"]
            self.assertEqual(solver["mobility"]["high_field_driving_force"], expected)
            self.assertEqual(solver["impact_ionization"]["driving_force"], expected)
            self.assertEqual(solver["impact_ionization"]["generation"], "current_density")
            self.assertEqual(solver["impact_ionization"]["current_approximation"], "density_gradient")

    def test_coupled_continuation_preserves_state_and_resolves_tiny_currents(self):
        result = emit(self.ir)
        self.assertEqual(result.decks[0]["solver"]["method"], "gummel_newton")
        for deck in result.decks[1:]:
            self.assertEqual(deck["solver"]["method"], "newton")
            self.assertTrue(deck["solver"]["warm_start"])
            self.assertEqual(deck["solver"]["quasi_fermi_reference"], "contact_majority")
            self.assertEqual(deck["solver"]["carrier_row_convergence"]["min_newton_max_iter"],
                             deck["solver"]["max_iter"])

    def test_save_load_restores_voltages_and_rebase_preserves_snapshot(self):
        ir = parse_cmd_text('''File { Grid="mesh.grd" Parameter="models.par" }
            Electrode { {Name="left" Voltage=0} {Name="right" Voltage=0} }
            Solve { Coupled {Poisson Electron Hole} Save(FilePrefix="eq")
              Quasistationary(Goal {Name="left" Voltage=0.1}) {Coupled {Poisson Electron Hole}}
              Load(FilePrefix="eq")
              Quasistationary(Goal {Name="left" Voltage=0.2}) {Coupled {Poisson Electron Hole}} }''', path_name="save.cmd")
        result = emit(ir)
        self.assertEqual(result.stage_plan["stages"][0]["output_state"], "eq_des.h5")
        self.assertEqual(result.stage_plan["stages"][2]["input_state"], "eq_des.h5")
        self.assertEqual(result.decks[2]["solver"]["method"], "frozen_state")
        rebased = rebase_stage(result, 1, {"left": .02})
        self.assertEqual(rebased.stage_plan["stages"][3]["sweep"]["start_V"], 0)
        ir["solve"]["statements"][3]["file_prefix"] = "absent"
        self.assert_error("E407", emit, ir)

    def test_fixed_charge_units_and_scopes(self):
        ir = deepcopy(self.ir)
        ir["physics"]["fixed_charges"] = [
            {"scope": "Region", "target": "bulk", "line": 1, "conc_cm3": 2e15},
            {"scope": "MaterialInterface", "target": "Silicon/Oxide", "line": 1, "conc_cm2": -2e10}]
        result = emit(ir, mesh_qualification=qualification(ir, ("Si", "SiO2")),
                      materials_report=material_result(ir, ("Si", "SiO2")).report)
        self.assertEqual(result.decks[0]["regions"], [{"name": "bulk", "fixed_charge_cm3": 2e15}])
        self.assertEqual(result.decks[0]["interfaces"], [{"regions": ["bulk", "oxide"], "fixed_charge_cm2": -2e10}])

    def test_cpp_catalogue_and_schema_unique_items(self):
        cpp = (REPO / "src/io/VtuWriter.cpp").read_text()
        block = cpp.split("static const std::vector<TdrFieldSpec> catalogue = {")[1].split("};")[0]
        self.assertEqual(tuple(re.findall(r'\{"([A-Za-z0-9]+)"', block)), TDR_FIELDS)
        schema = {"$defs": {"names": {"type": "array", "uniqueItems": True}}, "$ref": "#/$defs/names"}
        self.assertTrue(sv.is_valid(schema, [1, True]))
        self.assertFalse(sv.is_valid(schema, [1, 1.0]))
        self.assertFalse(sv.is_valid(schema, [{"a": 1}, {"a": 1}]))

    def test_signed_current_goal_and_missing_initial_current(self):
        ir = deepcopy(self.ir)
        stage = ir["stages"][-1]
        stage["goal"] = {"contact": "Anode", "current_A_per_um": -1e-7}
        stage["current_mode_contacts"] = ["Anode"]
        stage["outputs"]["plot_files"] = []
        result = emit(ir)
        control = result.decks[-1]["sweep"]["voltage_to_current"]
        self.assertEqual(control["current_points_A_per_um"], [1e-7])
        self.assertEqual(control["current_direction"], -1)
        self.assertEqual(result.stage_plan["stages"][-1]["sweep"]["stop_current_A_per_um"], -1e-7)
        stage["outputs"]["current_plot_times"] = {"list": [0, 1]}
        self.assert_error("E407", emit, ir)

    def test_malformed_input_has_a_numbered_error(self):
        self.assert_error("E407", emit_stage_plan, {}, mesh_qualification={}, materials_report={},
                          work_dir="device_vela", cmd_basename="bad.cmd", cmd_sha256="0" * 64)
        report = material_result(self.ir).report
        del report["solver_parameters"]
        self.assert_error("E407", emit, self.ir, materials_report=report)

    def test_writer_and_cli_roundtrip(self):
        result = emit(self.ir)
        with tempfile.TemporaryDirectory(dir=REPO / "build") as temp:
            root = Path(temp)
            write_stage_plan(root, result)
            self.assertEqual(json.loads((root / "stage_plan.json").read_text()), result.stage_plan)
            for name, data in (("ir", self.ir), ("mesh", qualification(self.ir)), ("report", material_result(self.ir).report)):
                (root / f"{name}.json").write_text(json.dumps(data))
            env = dict(os.environ, PYTHONPATH=str(REPO / "scripts"))
            done = subprocess.run([sys.executable, "-m", "vela_sdevice.deck_emit", "--ir", str(root / "ir.json"),
                "--mesh-qualification", str(root / "mesh.json"), "--materials-report", str(root / "report.json"),
                "--out", str(root / "cli"), "--allow-lossy"], cwd=REPO, env=env, capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertEqual(json.loads((root / "cli/stage_plan.json").read_text()), result.stage_plan)


class RunnerSmokeTests(unittest.TestCase):
    def test_biased_p_contact_pins_minority_quasi_fermi_and_density(self):
        runner = os.environ.get("VELA_EXAMPLE_RUNNER")
        if not runner or not Path(runner).is_file():
            self.skipTest("VELA_EXAMPLE_RUNNER is not available")
        ir = parse_cmd_text('''File { Grid="mesh.grd" Parameter="models.par" }
            Electrode { {Name="left" Voltage=0} {Name="right" Voltage=0} }
            Physics { Recombination(SRH) }
            Plot { eDensity hDensity eQuasiFermi hQuasiFermi }
            Solve { Coupled {Poisson Electron Hole}
              Quasistationary(Goal {Name="left" Voltage=0.2}) {Coupled {Poisson Electron Hole}} }''',
            path_name="ohmic.cmd")
        materials = material_result(ir)
        result = emit(ir, materials_report=materials.report)
        mesh = {"nodes": [{"id": i, "x": x, "y": y} for i, (x, y) in
                          enumerate(((0, 0), (1, 0), (1, 1), (0, 1)))],
            "triangles": [{"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                          {"id": 1, "region_id": 0, "node_ids": [0, 2, 3]}],
            "regions": [{"id": 0, "name": "bulk", "material": "Si", "cell_ids": [0, 1]}],
            "contacts": [{"id": 0, "name": "left", "region_id": 0, "node_ids": [0, 3], "edge_node_ids": [[0, 3]]},
                         {"id": 1, "name": "right", "region_id": 0, "node_ids": [1, 2], "edge_node_ids": [[1, 2]]}]}
        with tempfile.TemporaryDirectory(dir=REPO / "build") as temp:
            root = Path(temp)
            (root / "mesh.json").write_text(json.dumps(mesh))
            (root / "doping.csv").write_text("node_id,donors_cm3,acceptors_cm3\n" +
                                           "".join(f"{i},0,1e17\n" for i in range(4)))
            write_materials(root, materials)
            write_stage_plan(root, result)
            for stage in result.stage_plan["stages"]:
                done = subprocess.run([str(Path(runner).resolve()), "--config", stage["deck"], "--log", "off"],
                                      cwd=root, capture_output=True, text=True, timeout=60)
                self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
            import xml.etree.ElementTree as ET
            from sdevice_acceptance import array_values, vtu_metadata
            final = next(p for p in root.glob("stage01*.vtu") if abs(vtu_metadata(p)["Bias_V"] - .2) < 1e-9)
            fields = {v.attrib["Name"]: array_values(v) for v in ET.parse(final).findall(".//PointData/DataArray")}
            for name in ("eQuasiFermiPotential", "hQuasiFermiPotential"):
                for node in (0, 3):
                    self.assertEqual(len(fields[name][node]), 1)
                    self.assertAlmostEqual(fields[name][node][0], .2, delta=1e-12)
            for name in ("eDensity", "hDensity"):
                self.assertEqual(len(fields[name][0]), 1)
                self.assertEqual(len(fields[name][1]), 1)
                self.assertAlmostEqual(fields[name][0][0] / fields[name][1][0], 1., delta=1e-10)

    def test_equilibrium_voltage_and_current_goals_with_state_chain(self):
        runner = os.environ.get("VELA_EXAMPLE_RUNNER")
        if not runner or not Path(runner).is_file():
            self.skipTest("VELA_EXAMPLE_RUNNER is not available")
        ir = parse_cmd_text('''File { Grid="mesh.grd" Parameter="models.par" }
            Electrode { {Name="left" Voltage=0} {Name="right" Voltage=0} }
            Solve { Coupled {Poisson} Coupled {Poisson Electron Hole} Save(FilePrefix="eq")
              Quasistationary(InitialStep=0.1 MaxStep=0.25 Goal {Name="left" Voltage=0.001}) {
                Coupled {Poisson Electron Hole}}
              Set("left" mode current)
              Quasistationary(Goal {Name="left" Current=1e-7}) {Coupled {Poisson Electron Hole}}
              Load(FilePrefix="eq") }''', path_name="smoke.cmd")
        materials = material_result(ir)
        result = emit(ir, materials_report=materials.report,
                      sidecar={"schema": "vela.sdevice_sidecar.v1", "apply_to": [0, 1, 2], "solver": {"method": "newton"}})
        with tempfile.TemporaryDirectory(dir=REPO / "build") as temp:
            root = Path(temp)
            mesh = {"nodes": [{"id": 0, "x": 0, "y": 0}, {"id": 1, "x": 1, "y": 0},
                              {"id": 2, "x": 1, "y": 1}, {"id": 3, "x": 0, "y": 1}],
                    "triangles": [{"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                                  {"id": 1, "region_id": 0, "node_ids": [0, 2, 3]}],
                    "regions": [{"id": 0, "name": "bulk", "material": "Si", "cell_ids": [0, 1]}],
                    "contacts": [{"id": 0, "name": "left", "region_id": 0, "node_ids": [0, 3], "edge_node_ids": [[0, 3]]},
                                 {"id": 1, "name": "right", "region_id": 0, "node_ids": [1, 2], "edge_node_ids": [[1, 2]]}]}
            (root / "mesh.json").write_text(json.dumps(mesh))
            (root / "doping.csv").write_text("node_id,donors_cm3,acceptors_cm3\n" + "".join(f"{i},1e15,0\n" for i in range(4)))
            write_materials(root, materials)
            write_stage_plan(root, result)
            for stage in result.stage_plan["stages"]:
                done = subprocess.run([str(Path(runner).resolve()), "--config", stage["deck"], "--log", "off"],
                                      cwd=root, capture_output=True, text=True, timeout=120)
                self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
                self.assertTrue((root / stage["output_state"]).is_file())
                with (root / stage["csv"]).open() as stream:
                    rows = list(csv.DictReader(stream))
                self.assertGreater(len(rows), 0)
                self.assertTrue(any(root.glob("*.vtu")))
                if stage.get("sweep", {}).get("mode") == "current":
                    self.assertAlmostEqual(float(rows[-1]["contact_left_total_A_per_um"]), 1e-7, delta=5e-10)


if __name__ == "__main__":
    unittest.main()
