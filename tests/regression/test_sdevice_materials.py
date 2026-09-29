"""WP3 materials emitter tests (``vela_sdevice.materials``).

Acceptance (plan section 10, WP3):

1. PN2D ``models.par`` + IV active set reproduces the closed-form silicon
   scalars to 1e-12 and the checked-in IV reference ``ni`` to the documented
   tolerance (see ``ReferenceNiTest`` for why 1e-12 is not attainable there);
2. the generated default ``models.par`` round-trips the C++ built-ins;
3. a mesh material the parameter file does not define is E403;
4. ``unsupported_*`` active parameters are fatal, ``approximated`` ones are
   fatal without ``allow_lossy`` and recorded with it;
5. the output validates against the C3 schema.
"""

from __future__ import annotations

import json
import math
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

import generate_default_models_par as gen  # noqa: E402
import sentaurus_parameter_map as pmap  # noqa: E402
from vela_sdevice import schema_validator as sv  # noqa: E402
from vela_sdevice import materials as mat  # noqa: E402

FIXTURES = REPO / "tests" / "fixtures" / "sdevice"
PN2D_PAR = FIXTURES / "pn2d_sentaurus2018" / "models.par"
BV_PAR = FIXTURES / "bvmethods_sentaurus2018" / "models.par"
IV_IR = FIXTURES / "pn2d_sentaurus2018" / "pn2d_iv_sdevice.ir.json"
BV_IR = FIXTURES / "pn2d_sentaurus2018" / "pn2d_bv_sdevice.ir.json"
IV_REFERENCE = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source" / "pn2d_sentaurus2018_iv_materials.json"
SCHEMA = sv.load_schema(REPO / "schemas" / "vela.templates_ldmos.materials.v1.schema.json")
WORKFLOW_DOC = REPO / "docs" / "sentaurus_sdevice_workflow.md"
CATALOGUE_ROW = re.compile(r"^\| (?P<id>[WE][1-8]\d{2}) \| (?P<severity>warning|fatal) \| (?P<exit>[0-8]|\u2014) \| `(?P<template>[^`]*)` \|")

KB = 1.380649e-23
Q = 1.602176634e-19
H = 6.62607015e-34
M0 = 9.1093837015e-31
KT300 = KB / Q * 300.0


def _tmpdir() -> tempfile.TemporaryDirectory:
    base = REPO / "build"
    return tempfile.TemporaryDirectory(dir=base if base.is_dir() else None, prefix="sdevice_materials_")


def _ir(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _emit(par: Path, ir_path: Path, materials=("Si",), **kwargs) -> mat.EmitResult:
    ir = _ir(ir_path)
    return mat.emit_materials(par, active_models=mat.active_models_from_ir(ir),
                              temperature_K=ir["physics"]["temperature_K"], mesh_materials=list(materials),
                              effective_models=ir["physics"]["effective"], **kwargs)


def _material(result: mat.EmitResult, name: str) -> dict:
    return next(m for m in result.materials["materials"] if m["name"] == name)


def _write_par(directory: Path, text: str, name: str = "user.par") -> Path:
    path = directory / name
    path.write_text(text, encoding="utf-8")
    return path


def pn2d_silicon_closed_form(fold_old_slotboom: bool, T: float = 300.0) -> dict:
    """Independent evaluation of the PN2D Silicon parameters (values copied from models.par)."""
    eg0, alpha, beta, tpar = 1.16964, 4.73e-4, 636.0, 0.0
    eg_zero = eg0 + alpha * tpar ** 2 / (beta + tpar)
    eg = eg_zero - alpha * T * T / (beta + T)
    me = ((6.0 * 0.1905 * eg_zero / eg) ** 2 * 0.9163) ** (1.0 / 3.0)
    num = 0.443587 + 0.003609528 * T + 0.0001173515 * T ** 2 + 1.263218e-6 * T ** 3 + 3.025581e-9 * T ** 4
    den = 1.0 + 0.004683382 * T + 0.0002286895 * T ** 2 + 7.469271e-7 * T ** 3 + 1.727481e-9 * T ** 4
    mh = (num / den) ** (2.0 / 3.0)
    pref = 2.0 * (2.0 * math.pi * M0 * KB * T / H ** 2) ** 1.5 / 1.0e6
    nc, nv = pref * me ** 1.5, pref * mh ** 1.5
    gap = eg + (-0.01595 if fold_old_slotboom else 0.0)
    ni = math.sqrt(nc * nv) * math.exp(-gap / (2.0 * KB / Q * T))
    return {"bandgap_eV": eg, "Nc": nc, "Nv": nv, "ni": ni}


class Pn2dSiliconTest(unittest.TestCase):
    """Acceptance 1 and 5: PN2D parameter dump with the IV active set."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.result = _emit(PN2D_PAR, IV_IR)
        cls.si = _material(cls.result, "Si")

    def test_output_validates_against_c3(self) -> None:
        self.assertEqual([], sv.collect_errors(SCHEMA, self.result.materials))
        self.assertEqual("vela.templates_ldmos.materials.v1", self.result.materials["schema"])
        self.assertEqual(mat.UNIT_SYSTEM, self.result.materials["unit_system"])

    def test_active_models(self) -> None:
        self.assertEqual(["DopingDependence", "OldSlotboom", "SRH"], self.result.report["active_models"])

    def test_scalars_match_closed_form(self) -> None:
        expected = pn2d_silicon_closed_form(fold_old_slotboom=True)
        self.assertEqual("transport_semiconductor", self.si["role"])
        self.assertEqual(11.7, self.si["eps_r"])
        self.assertEqual(4.05, self.si["electron_affinity_eV"])
        self.assertEqual(1417.0, self.si["electron_mobility_cm2_per_V_s"])
        self.assertEqual(470.5, self.si["hole_mobility_cm2_per_V_s"])
        self.assertEqual(300.0, self.si["temperature_K"])
        for key, field in (("bandgap_eV", "bandgap_eV"), ("Nc", "conduction_band_density_of_states_cm3"),
                           ("Nv", "valence_band_density_of_states_cm3"), ("ni", "intrinsic_carrier_density_cm3")):
            self.assertAlmostEqual(expected[key] / self.si[field], 1.0, delta=1e-12, msg=field)
        # dEg0(OldSlotboom) lives in ni only; the gap stays the plain Varshni value.
        self.assertAlmostEqual(1.1241592307692307, self.si["bandgap_eV"], delta=1e-15)
        self.assertIn("dEg0(OldSlotboom)=-0.01595 eV folded into ni", self.si["provenance"])

    def test_old_slotboom_fold_is_the_only_ni_difference(self) -> None:
        ir = _ir(IV_IR)
        without = mat.emit_materials(PN2D_PAR, active_models=["DopingDependence", "SRH"], temperature_K=300.0,
                                     mesh_materials=["Si"])
        ratio = self.si["intrinsic_carrier_density_cm3"] / _material(without, "Si")["intrinsic_carrier_density_cm3"]
        self.assertAlmostEqual(math.exp(0.01595 / (2.0 * KT300)), ratio, delta=1e-12)
        self.assertEqual(ir["physics"]["effective"], self.result.report["effective_models"])

    def test_library_dump_needs_no_lossy_flag(self) -> None:
        report = self.result.report
        self.assertEqual([], report["lossy"])
        inert = {(x["section"], x["parameter"]) for x in report["inert"]}
        self.assertEqual({("Bandgap", "Bgn2Chi"), ("Bandgap", "alpha2"), ("Bandgap", "beta2"), ("Bandgap", "EgMin"),
                          ("Bandgap", "dEgMin"), ("ConstantMobility", "Exponent"), ("ConstantMobility", "mutunnel")},
                         inert)
        self.assertNotIn("unmapped", report["counts"])
        self.assertGreater(report["counts"]["inactive"], 1000)

    def test_solver_parameters_are_routed_in_si(self) -> None:
        rows = {(r["target"]): r for r in self.result.report["solver_parameters"]}
        self.assertAlmostEqual(52.2e-4, rows["solver.mobility.electron_mumin1_m2_V_s"]["value_si"], delta=1e-18)
        self.assertEqual(1.0e23, rows["solver.bandgap_narrowing.reference_doping_m3"]["value_si"])
        self.assertEqual(1.0e22, rows["solver.srh_doping_dependence.electron.reference_doping_m3"]["value_si"])
        self.assertEqual(1.0e-5, rows["solver.srh_doping_dependence.electron.tau_max_s"]["value_si"])
        self.assertEqual(3.0e-6, rows["solver.srh_doping_dependence.hole.tau_max_s"]["value_si"])
        self.assertEqual("hole", rows["solver.mobility.hole_pc_m3"]["carrier"])
        self.assertEqual({"exact"}, {r["status"] for r in rows.values()})
        self.assertFalse(any(r["target"].startswith("solver.impact_ionization") for r in rows.values()))

    def test_bv_active_set_routes_high_field_and_van_overstraeten(self) -> None:
        result = _emit(PN2D_PAR, BV_IR)
        self.assertEqual(["DopingDependence", "HighFieldSaturation", "OldSlotboom", "SRH", "VanOverstraeten"],
                         result.report["active_models"])
        rows = {r["target"]: r for r in result.report["solver_parameters"]}
        a_low = rows["solver.impact_ionization.electron_a_low_m_inv"]
        self.assertEqual("m^-1", a_low["unit_si"])
        self.assertAlmostEqual(a_low["value_par"] * 1.0e2, a_low["value_si"], delta=1e-6)
        self.assertEqual("V/m", rows["solver.impact_ionization.switch_field_V_m"]["unit_si"])
        self.assertEqual("eV", rows["solver.impact_ionization.phonon_energy_eV"]["unit_si"])
        vsat = rows["solver.mobility.electron_saturation_velocity_m_s"]
        self.assertEqual("m/s", vsat["unit_si"])
        self.assertAlmostEqual(vsat["value_par"] * 1.0e-2, vsat["value_si"], delta=1e-9)
        self.assertEqual("frozen_at_temperature", vsat["status"])
        inert = {(x["section"], x["parameter"]) for x in result.report["inert"]}
        self.assertIn(("HighFieldDependence", "betaexp"), inert)
        self.assertIn(("HighFieldDependence", "vsatexp"), inert)
        self.assertEqual([], result.report["lossy"])


class ReferenceNiTest(unittest.TestCase):
    """The checked-in IV reference ``ni`` (1.4638914958767616e10 cm^-3).

    The June 2026 PN2D work derived that literal from ``models.par`` (plan
    2026-06-13, commit e2d9a1c8) but did not record the computation; it also
    equals the Sentaurus contact-edge ``ni_eff`` (~1.65562e10, doping 1e17)
    divided by the OldSlotboom factor at Nref, i.e. Sentaurus' internal value.
    No evaluation of the documented laws with CODATA (1969-2018), rounded or
    Sentaurus (q = 1.602192e-19) constants reproduces it closer than ~2.7e-4.
    The emitter therefore pins the closed form at 1e-12 (``Pn2dSiliconTest``)
    and checks the reference at 5e-4, freezing the observed deviation.
    """

    def test_reference_ni_within_documented_tolerance(self) -> None:
        reference = json.loads(IV_REFERENCE.read_text(encoding="utf-8"))
        ref_ni = next(m for m in reference["materials"] if m["name"] == "Si")["ni"]
        si = _material(_emit(PN2D_PAR, IV_IR), "Si")
        deviation = si["intrinsic_carrier_density_cm3"] / ref_ni - 1.0
        self.assertLess(abs(deviation), 5e-4, f"ni deviates {deviation:.3e} from the IV reference")
        self.assertAlmostEqual(-2.74e-4, deviation, delta=2e-6, msg="documented deviation changed")


class DefaultLibraryTest(unittest.TestCase):
    """Acceptance 2: ``share/vela/materials/models.par`` (``vela device -P``).

    The library carries the published silicon defaults Sentaurus also ships,
    so it must reproduce the PN2D ``sdevice -P`` dump physics exactly and
    provide a fallback section for every model the emitter supports.
    """

    ALL_MODELS = ["DopingDependence", "HighFieldSaturation", "SRH", "TempDependence", "Auger", "OldSlotboom",
                  "VanOverstraeten", "eQuantumPotential"]

    @classmethod
    def setUpClass(cls) -> None:
        cls.result = mat.emit_materials(mat.default_models_par_path(), active_models=[], temperature_K=300.0,
                                        mesh_materials=["Si", "SiO2", "Nitride"])

    def test_generated_file_is_current(self) -> None:
        self.assertEqual(0, gen.main(["--check"]))
        text = mat.default_models_par_path().read_text(encoding="utf-8")
        self.assertNotIn("Synopsys", text)
        self.assertNotIn("Copyright", text)

    def test_silicon_matches_the_sentaurus_default_dump(self) -> None:
        for models in ([], ["OldSlotboom"]):
            library = _material(mat.emit_materials(mat.default_models_par_path(), active_models=models,
                                                   temperature_K=300.0, mesh_materials=["Si"]), "Si")
            dump = _material(mat.emit_materials(PN2D_PAR, active_models=models, temperature_K=300.0,
                                                mesh_materials=["Si"]), "Si")
            for key, value in dump.items():
                if key != "provenance":
                    self.assertEqual(value, library[key], f"{key} with {models}")
        si = _material(self.result, "Si")
        expected = pn2d_silicon_closed_form(fold_old_slotboom=False)
        self.assertAlmostEqual(expected["ni"] / si["intrinsic_carrier_density_cm3"], 1.0, delta=1e-12)
        self.assertAlmostEqual(1.0750e10, si["intrinsic_carrier_density_cm3"], delta=1e6)
        self.assertEqual(1417.0, si["electron_mobility_cm2_per_V_s"])
        self.assertEqual([], sv.collect_errors(SCHEMA, self.result.materials))

    def test_dielectrics_follow_the_sentaurus_defaults(self) -> None:
        sio2 = _material(self.result, "SiO2")
        nitride = _material(self.result, "Nitride")
        self.assertEqual({"name": "SiO2", "role": "dielectric", "eps_r": 3.9, "bandgap_eV": 9.0,
                          "electron_affinity_eV": 0.9, "temperature_K": 300.0},
                         {k: v for k, v in sio2.items() if k != "provenance"})
        self.assertEqual({"name": "Nitride", "role": "dielectric", "eps_r": 7.5, "bandgap_eV": 5.0,
                          "electron_affinity_eV": 1.9, "temperature_K": 300.0},
                         {k: v for k, v in nitride.items() if k != "provenance"})

    def test_every_supported_model_has_a_lossless_fallback(self) -> None:
        for temperature in (300.0, 350.0):
            result = mat.emit_materials(mat.default_models_par_path(), active_models=self.ALL_MODELS,
                                        temperature_K=temperature, mesh_materials=["Si"])
            self.assertEqual([], result.report["lossy"], temperature)
            self.assertNotIn("unmapped", result.report["counts"])
            targets = {r["target"].split(".")[1] for r in result.report["solver_parameters"]}
            self.assertEqual({"auger_cn_m6_per_s", "auger_cp_m6_per_s", "bandgap_narrowing",
                              "electron_quantum_potential", "impact_ionization", "mobility",
                              "srh_doping_dependence"}, targets)
            self.assertEqual(3.6, _material(result, "Si")["electron_quantum_gamma"])
        # Vela-representable substitutions for Sentaurus defaults are explicit in the library.
        text = mat.default_models_par_path().read_text(encoding="utf-8")
        self.assertRegex(text, r"Exponent\s*=\s*2\.2 ,\s*2\.2")
        self.assertRegex(text, r"\bH\s*=\s*0 ,\s*0")
        self.assertNotRegex(text, r"^\s*mutunnel\s*=", msg="mutunnel must not be a library parameter")

    def test_builtin_cpp_silicon_is_documented_as_legacy(self) -> None:
        """The C++ MaterialDatabase Si (ni 1e10, Nc 2.8e19, Nv 1.04e19, Eg 1.12) is not the library source."""
        si = _material(self.result, "Si")
        self.assertNotAlmostEqual(1.0e10, si["intrinsic_carrier_density_cm3"], delta=1e8)
        self.assertNotEqual(1.12, si["bandgap_eV"])
        self.assertIn("MaterialDatabase", gen.__doc__)

    def test_partial_user_file_falls_back_per_parameter(self) -> None:
        result = mat.emit_materials(BV_PAR, active_models=["SRH"], temperature_K=300.0, mesh_materials=["Si", "SiO2"])
        si_report = next(m for m in result.report["materials"] if m["name"] == "Si")
        self.assertFalse(si_report["defined_in_parameter_file"])
        self.assertEqual(["user_global", "vela_default"], si_report["sections"]["Scharfetter"]["origins"])
        self.assertEqual(["vela_default"], si_report["sections"]["Bandgap"]["origins"])
        taumax = next(r for r in result.report["solver_parameters"]
                      if r["target"] == "solver.srh_doping_dependence.electron.tau_max_s")
        self.assertEqual(1e-7, taumax["value_si"])
        self.assertEqual("user_global", taumax["origin"])
        self.assertEqual(_material(self.result, "Si")["conduction_band_density_of_states_cm3"],
                         _material(result, "Si")["conduction_band_density_of_states_cm3"])
        # The user's global Scharfetter block wins over the library's Silicon block parameter by parameter.
        self.assertNotEqual(1e-5, taumax["value_si"])


class FailClosedTest(unittest.TestCase):
    """Acceptance 3 and 4."""

    def test_missing_mesh_material_is_e403(self) -> None:
        with self.assertRaises(mat.MaterialsError) as ctx:
            mat.emit_materials(BV_PAR, active_models=["SRH"], temperature_K=300.0, mesh_materials=["Si", "HfO2"])
        self.assertEqual("E403", ctx.exception.message_id)
        self.assertEqual("Material 'HfO2' is neither built in nor defined in the parameter file",
                         ctx.exception.message)

    def test_user_defined_dielectric_and_library_fallback(self) -> None:
        text = ('Material = "HfO2" {\nEpsilon\n{\n  epsilon = 25  # [1]\n}\n}\n'
                'Material = "Nitride" {\nEpsilon\n{\n  epsilon = 7.0  # [1]\n}\n}\n')
        with _tmpdir() as tmp:
            par = _write_par(Path(tmp), text)
            result = mat.emit_materials(par, active_models=[], temperature_K=300.0,
                                        mesh_materials=["Si", "HfO2", "Nitride"])
        hfo2 = _material(result, "HfO2")
        self.assertEqual({"name": "HfO2", "role": "dielectric", "eps_r": 25.0, "temperature_K": 300.0},
                         {k: v for k, v in hfo2.items() if k != "provenance"})
        # The user's Epsilon overrides the library; the library's Bandgap fills in.
        nitride = _material(result, "Nitride")
        self.assertEqual({"name": "Nitride", "role": "dielectric", "eps_r": 7.0, "bandgap_eV": 5.0,
                          "electron_affinity_eV": 1.9, "temperature_K": 300.0},
                         {k: v for k, v in nitride.items() if k != "provenance"})
        self.assertEqual([], sv.collect_errors(SCHEMA, result.materials))

    def test_unsupported_active_parameter_is_e404(self) -> None:
        text = 'Material = "Silicon" {\nScharfetter\n{\n  taumin = 0 , 0  # [s]\n  taumax = 1e-5 , 3e-6  # [s]\n' \
               '  Nref = 1e16 , 1e16  # [cm^(-3)]\n  gamma = 1 , 1  # [1]\n  Etrap = 0.1  # [eV]\n}\n}\n'
        with _tmpdir() as tmp:
            par = _write_par(Path(tmp), text)
            with self.assertRaises(mat.MaterialsError) as ctx:
                mat.emit_materials(par, active_models=["SRH"], temperature_K=300.0, mesh_materials=["Si"])
            self.assertEqual("E404", ctx.exception.message_id)
            self.assertEqual("Parameter 'Si.Scharfetter.Etrap' of active model 'SRH' cannot be imported "
                             "(unsupported_model)", ctx.exception.message)
            # The same file is fine when SRH is not active: the parameter is inert.
            result = mat.emit_materials(par, active_models=[], temperature_K=300.0, mesh_materials=["Si"])
            self.assertIn("inactive", result.report["counts"])

    def test_approximated_needs_allow_lossy(self) -> None:
        text = 'Material = "Silicon" {\nConstantMobility:\n{\n  mumax = 1417 , 470.5  # [cm^2/(Vs)]\n' \
               '  Exponent = 2.5 , 2.2  # [1]\n}\n}\n'
        with _tmpdir() as tmp:
            par = _write_par(Path(tmp), text)
            # At 300 K the exponent is inert whatever its value.
            ok = mat.emit_materials(par, active_models=[], temperature_K=300.0, mesh_materials=["Si"])
            self.assertIn(("ConstantMobility", "Exponent"), {(x["section"], x["parameter"]) for x in ok.report["inert"]})
            with self.assertRaises(mat.MaterialsError) as ctx:
                mat.emit_materials(par, active_models=[], temperature_K=350.0, mesh_materials=["Si"])
            self.assertEqual("E404", ctx.exception.message_id)
            self.assertIn("(approximated)", ctx.exception.message)
            lossy = mat.emit_materials(par, active_models=[], temperature_K=350.0, mesh_materials=["Si"],
                                       allow_lossy=True)
            self.assertEqual([("ConstantMobility", "Exponent")],
                             [(x["section"], x["parameter"]) for x in lossy.report["lossy"]])
            self.assertTrue(lossy.report["allow_lossy"])
            self.assertEqual(1417.0, _material(lossy, "Si")["electron_mobility_cm2_per_V_s"])

    def test_bgn2chi_is_lossy_when_band_edges_are_observable(self) -> None:
        ir = _ir(IV_IR)
        models = mat.active_models_from_ir(ir)
        with self.assertRaises(mat.MaterialsError) as ctx:
            mat.emit_materials(PN2D_PAR, active_models=models, temperature_K=300.0, mesh_materials=["Si"],
                               band_edge_sensitive=True)
        self.assertEqual("E404", ctx.exception.message_id)
        self.assertIn("Si.Bandgap.Bgn2Chi", ctx.exception.message)
        result = mat.emit_materials(PN2D_PAR, active_models=models, temperature_K=300.0, mesh_materials=["Si"],
                                    band_edge_sensitive=True, allow_lossy=True)
        self.assertEqual([("Bandgap", "Bgn2Chi")], [(x["section"], x["parameter"]) for x in result.report["lossy"]])

    def test_missing_section_for_active_model_is_e404(self) -> None:
        with _tmpdir() as tmp:
            par = _write_par(Path(tmp), 'Material = "Silicon" {\nEpsilon\n{\n  epsilon = 11.7  # [1]\n}\n}\n')
            # Without the default library nothing can supply the SRH lifetimes.
            with self.assertRaises(mat.MaterialsError) as ctx:
                mat.emit_materials(par, active_models=["SRH"], temperature_K=300.0, mesh_materials=["Si"],
                                   default_par=Path(tmp) / "absent.par")
            self.assertEqual("E404", ctx.exception.message_id)
            self.assertEqual("Parameter 'Si.Scharfetter' of active model 'SRH' cannot be imported "
                             "(missing_section)", ctx.exception.message)
            # With the library the same file imports and records the fallback origin.
            result = mat.emit_materials(par, active_models=["SRH"], temperature_K=300.0, mesh_materials=["Si"])
            si_report = next(m for m in result.report["materials"] if m["name"] == "Si")
            self.assertEqual(["vela_default"], si_report["sections"]["Scharfetter"]["origins"])
            self.assertEqual(["user_material"], si_report["sections"]["Epsilon"]["origins"])

    def test_region_scope_with_active_parameters_is_e404(self) -> None:
        text = 'Region = "R.Si" {\nScharfetter\n{\n  taumax = 1e-6 , 1e-6  # [s]\n}\n}\n'
        with _tmpdir() as tmp:
            par = _write_par(Path(tmp), text)
            with self.assertRaises(mat.MaterialsError) as ctx:
                mat.emit_materials(par, active_models=["SRH"], temperature_K=300.0, mesh_materials=["Si"])
            self.assertEqual("E404", ctx.exception.message_id)
            self.assertIn("region:R.Si.Scharfetter.taumax", ctx.exception.message)
            self.assertIn("(unsupported_scope)", ctx.exception.message)
            result = mat.emit_materials(par, active_models=[], temperature_K=300.0, mesh_materials=["Si"])
            self.assertEqual([{"kind": "region", "name": "R.Si", "section": "Scharfetter", "line": 2}],
                             result.report["ignored_scopes"])

    def test_missing_or_unparsable_parameter_file(self) -> None:
        with _tmpdir() as tmp:
            with self.assertRaises(mat.MaterialsError) as ctx:
                mat.emit_materials(Path(tmp) / "absent.par", active_models=[], temperature_K=300.0,
                                   mesh_materials=["Si"])
            self.assertEqual("E303", ctx.exception.message_id)
            par = _write_par(Path(tmp), 'Material = "Silicon" {\nEpsilon\n{\n  epsilon 11.7\n}\n}\n')
            with self.assertRaises(mat.MaterialsError) as ctx:
                mat.emit_materials(par, active_models=[], temperature_K=300.0, mesh_materials=["Si"])
            self.assertEqual("E410", ctx.exception.message_id)

    def test_auger_polynomial_collapses_to_frozen_constants(self) -> None:
        text = ('Material = "Silicon" {\nAuger\n{\n  A = 0.67e-31 , 0.72e-31  # [cm^6/s]\n'
                '  B = 2.45e-31 , 4.50e-33  # [cm^6/s]\n  C = -2.2e-32 , 2.63e-32  # [cm^6/s]\n'
                '  H = 0 , 0  # [1]\n  N0 = 1e18 , 1e18  # [cm^-3]\n}\n}\n')
        with _tmpdir() as tmp:
            par = _write_par(Path(tmp), text)
            result = mat.emit_materials(par, active_models=["Auger"], temperature_K=300.0, mesh_materials=["Si"])
            rows = {(r["parameter"], r["carrier"]): r for r in result.report["solver_parameters"] if r["section"] == "Auger"}
            collapsed = rows[("A+B(T/300)+C(T/300)^2", "electron")]
            self.assertEqual("solver.auger_cn_m6_per_s", collapsed["target"])
            self.assertAlmostEqual(2.9e-43, collapsed["value_si"], delta=1e-55)  # Vela's Cn = 2.9e-31 cm^6/s
            self.assertAlmostEqual(1.028e-43, rows[("A+B(T/300)+C(T/300)^2", "hole")]["value_si"], delta=1e-55)
            self.assertEqual("frozen_at_temperature", collapsed["status"])
            self.assertEqual([], result.report["lossy"])
            inert = {(x["section"], x["parameter"]) for x in result.report["inert"]}
            self.assertIn(("Auger", "H"), inert)
            self.assertIn(("Auger", "N0"), inert)
            # H != 0 is a real carrier-density enhancement Vela lacks: fatal even with allow_lossy.
            par2 = _write_par(Path(tmp), text.replace("H = 0 , 0", "H = 3.46667 , 8.25688"), "enh.par")
            with self.assertRaises(mat.MaterialsError) as ctx:
                mat.emit_materials(par2, active_models=["Auger"], temperature_K=300.0, mesh_materials=["Si"],
                                   allow_lossy=True)
            self.assertIn("Si.Auger.H", ctx.exception.message)

    def test_electrode_work_function(self) -> None:
        text = 'Material = "Gold" {\nBandgap\n{\n  WorkFunction = 5.1  # [eV]\n}\n}\n'
        with _tmpdir() as tmp:
            par = _write_par(Path(tmp), text)
            result = mat.emit_materials(par, active_models=[], temperature_K=300.0, mesh_materials=["Si"],
                                        electrode_materials=["Gold"])
            self.assertEqual([("Gold", 5.1, "user_material")],
                             [(e["material"], e["work_function_eV"], e["origin"]) for e in result.report["electrodes"]])
            with self.assertRaises(mat.MaterialsError) as ctx:
                mat.emit_materials(par, active_models=[], temperature_K=300.0, mesh_materials=["Si"],
                                   electrode_materials=["Aluminum"])
            self.assertEqual("E403", ctx.exception.message_id)
            bare = _write_par(Path(tmp), 'Material = "Gold" {\nEpsilon\n{\n  epsilon = 1  # [1]\n}\n}\n', "bare.par")
            with self.assertRaises(mat.MaterialsError) as ctx:
                mat.emit_materials(bare, active_models=[], temperature_K=300.0, mesh_materials=["Si"],
                                   electrode_materials=["Gold"])
            self.assertEqual("E404", ctx.exception.message_id)
            self.assertIn("Gold.Bandgap.WorkFunction", ctx.exception.message)


class ContractTest(unittest.TestCase):
    def test_message_templates_match_catalogue(self) -> None:
        rows = {}
        for line in WORKFLOW_DOC.read_text(encoding="utf-8").splitlines():
            match = CATALOGUE_ROW.match(line)
            if match:
                rows[match["id"]] = match
        for message_id, (severity, template) in mat.TEMPLATES.items():
            self.assertIn(message_id, rows, f"{message_id} missing from docs/sentaurus_sdevice_workflow.md")
            self.assertEqual(rows[message_id]["severity"], severity, message_id)
            self.assertEqual(rows[message_id]["template"], template, message_id)

    def test_constants_match_physical_constants_header(self) -> None:
        header = (REPO / "include" / "vela" / "core" / "PhysicalConstants.h").read_text(encoding="utf-8")
        for name, key in (("kb", "kb_J_per_K"), ("q", "q_C"), ("h", "h_J_s"), ("m0", "m0_kg")):
            match = re.search(r"constexpr double " + name + r"\s*=\s*([\d.eE+-]+)", header)
            self.assertIsNotNone(match, name)
            self.assertEqual(float(match[1]), mat.CONSTANTS[key], name)

    def test_active_models_from_all_goldens(self) -> None:
        expected = {
            "pn2d_iv_sdevice.ir.json": ["DopingDependence", "OldSlotboom", "SRH"],
            "pn2d_bv_sdevice.ir.json": ["DopingDependence", "HighFieldSaturation", "OldSlotboom", "SRH",
                                        "VanOverstraeten"],
            "aba_poisson_sdevice.ir.json": ["Fermi", "VanOverstraeten"],
            "external_resistor_sdevice.ir.json": ["Band2BandTunneling", "DopingDependence", "Enormal",
                                                  "EnormalDependence", "Fermi", "HighFieldSaturation",
                                                  "OldSlotboom", "SRH", "VanOverstraeten"],
            "bjt_m0_des.ir.json": ["SRH"],
        }
        for name, models in expected.items():
            path = next(FIXTURES.rglob(name))
            self.assertEqual(models, mat.active_models_from_ir(_ir(path)), name)
        self.assertIn("TempDependence", mat.active_models_from_ir(
            {"physics": {"effective": [{"key": "Recombination.SRH", "parameters": {"temperature_dependence": True}}]}}))

    def test_matrix_statuses_are_the_only_inputs(self) -> None:
        """Every emitted status is a matrix status; refinements only demote to exact/frozen."""
        result = _emit(PN2D_PAR, BV_IR)
        statuses = set(result.report["counts"])
        self.assertTrue(statuses <= set(pmap.STATUSES) | {pmap.STATUS_INACTIVE})
        for item in result.report["inert"]:
            self.assertEqual(pmap.STATUS_EXACT, item["status"])
            self.assertNotEqual(item["matrix_status"], item["status"])
            self.assertIn("inert here", item["reason"])


class CliTest(unittest.TestCase):
    def run_cli(self, *args: str) -> subprocess.CompletedProcess:
        env = dict(os.environ, PYTHONPATH=str(REPO / "scripts"), PYTHONIOENCODING="utf-8")
        return subprocess.run([sys.executable, "-m", "vela_sdevice.materials", *args], capture_output=True,
                              text=True, env=env, cwd=str(REPO))

    def test_cli_writes_c3_and_report(self) -> None:
        with _tmpdir() as tmp:
            out = Path(tmp) / "out"
            proc = self.run_cli("--par", str(PN2D_PAR), "--ir", str(IV_IR), "--materials", "Si", "--out", str(out))
            self.assertEqual(0, proc.returncode, proc.stderr)
            summary = json.loads(proc.stdout)
            materials = json.loads((out / "materials.json").read_text(encoding="utf-8"))
            report = json.loads((out / "materials_report.json").read_text(encoding="utf-8"))
            self.assertEqual([], sv.collect_errors(SCHEMA, materials))
            self.assertEqual("vela.sdevice_materials_report.v1", report["schema"])
            self.assertEqual(summary["counts"], report["counts"])
            self.assertEqual("sdevice:models", materials["benchmark"])

    def test_cli_exit_codes(self) -> None:
        with _tmpdir() as tmp:
            out = Path(tmp) / "out"
            proc = self.run_cli("--par", str(BV_PAR), "--model", "SRH", "--materials", "Si", "HfO2", "--out", str(out))
            self.assertEqual(4, proc.returncode)
            self.assertIn("E403", proc.stderr)
            proc = self.run_cli("--par", str(Path(tmp) / "absent.par"), "--materials", "Si", "--out", str(out))
            self.assertEqual(3, proc.returncode)
            self.assertIn("E303", proc.stderr)
            proc = self.run_cli("--par", str(BV_PAR), "--out", str(out))
            self.assertEqual(2, proc.returncode)


if __name__ == "__main__":
    unittest.main()
