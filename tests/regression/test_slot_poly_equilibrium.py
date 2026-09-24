"""Synthetic restart invariants; these tests do not qualify a device solution."""
import copy
from tests.regression import slot_state_fixture as h5
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest


HELPER = Path(__file__).resolve().parents[2] / "scripts/slot_poly_equilibrium.py"


class PolyEquilibriumTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(HELPER.exists(), "guarded PolySilicon restart helper is missing")
        spec = importlib.util.spec_from_file_location("slot_poly_equilibrium", HELPER)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.write_restart = module.write_poly_equilibrium_restart
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "state.h5"
        self.output = self.root / "candidate.h5"
        self.mesh = {
            "nodes": [{"id": i, "x": i, "y": i % 2} for i in range(9)],
            "regions": [{"id": 0, "material": "PolySilicon"},
                        {"id": 1, "material": "Si"}, {"id": 2, "material": "SiO2"}],
            "triangles": [{"id": 0, "region_id": 0, "node_ids": [0, 1, 2]},
                          {"id": 1, "region_id": 0, "node_ids": [3, 4, 5]},
                          {"id": 2, "region_id": 1, "node_ids": [6, 7, 8]},
                          {"id": 3, "region_id": 2, "node_ids": [2, 5, 6]}],
            "contacts": [{"name": "gate", "node_ids": [0]},
                         {"name": "SLOT", "node_ids": [3]},
                         {"name": "drain", "node_ids": [8]}],
        }
        self.doc = {"mesh_file": "mesh.json", "materials_file": "materials.json",
                    "simulation_type": "dc_sweep", "sweep": {"contact": "drain"},
                    "contacts": [{"name": "gate", "type": "ohmic", "bias": 0.1},
                                 {"name": "SLOT", "type": "ohmic", "bias": -0.1},
                                 {"name": "drain", "type": "ohmic", "bias": 0.0}],
                    "solver": {"method": "newton", "recombination_materials": ["Si"],
                               "mobility": {"model": "ialmob", "ialmob": {
                                   "constant_material_mobility_cm2_per_Vs": {
                                       "PolySilicon": {"electron": 1417, "hole": 470.5}}}}}}
        (self.root / "materials.json").write_text(json.dumps({"materials": [
            {"name": "PolySilicon", "ni": 1e16, "mun": 0.1, "mup": 0.04}]}))
        self.fields = ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3",
                       "electron_quantum_potential_V", "electron_qf_increment_V",
                       "hole_qf_increment_V", "electron_qf_reference_V", "hole_qf_reference_V"]
        self.rows = [dict(zip(self.fields, map(str, [i, 0.3, 0.12, -0.08, 2e20, 3e10,
                                                   0.002, 0.02, 0.02, 0.1, -0.1])))
                     for i in range(9)]

    def persist(self):
        (self.root / "mesh.json").write_text(json.dumps(self.mesh))
        h5.write(self.source, self.rows, self.mesh)

    def run_helper(self, bias=2.0):
        self.persist()
        return self.write_restart(self.doc, self.source, self.output, self.root, bias)

    def read_output(self):
        return h5.read(self.output)

    def test_components_preserve_input_nonpoly_psi_and_quantum(self):
        self.persist()
        original = self.source.read_bytes()
        document = copy.deepcopy(self.doc)
        result = self.write_restart(self.doc, self.source, self.output, self.root, 2.0)
        out = self.read_output()
        self.assertEqual(self.source.read_bytes(), original)
        self.assertEqual(self.doc, document)
        self.assertEqual(out[6:], h5.read(self.source)[6:])
        for i in range(6):
            for field in ("psi", "electron_quantum_potential_V"):
                self.assertEqual(out[i][field], float(self.rows[i][field]))
            qf = 0.1 if i < 3 else -0.1
            self.assertEqual(float(out[i]["phin"]), qf)
            self.assertEqual(float(out[i]["phip"]), qf)
            for carrier in ("electron", "hole"):
                self.assertEqual(float(out[i][carrier + "_qf_reference_V"]), qf)
                self.assertEqual(float(out[i][carrier + "_qf_increment_V"]), 0.0)
        self.assertEqual([c["nodes"] for c in result["components"]], [3, 3])
        self.assertAlmostEqual(result["max_qf_change_V"], 0.22)

    def test_boltzmann_density_transform_both_signs(self):
        self.run_helper()
        vt = 1.380649e-23 * 300 / 1.602176634e-19
        for i, row in enumerate(self.read_output()[:6]):
            target = 0.1 if i < 3 else -0.1
            self.assertAlmostEqual(float(row["electrons_m3"]) / (2e20 * math.exp((0.12-target)/vt)), 1)
            self.assertAlmostEqual(float(row["holes_m3"]) / (3e10 * math.exp((target+0.08)/vt)), 1)

    def test_sweep_contact_bias_overrides_static_bias(self):
        self.doc["sweep"]["contact"] = "gate"
        self.run_helper(bias=0.2)
        self.assertEqual(float(self.read_output()[0]["phin"]), 0.2)

    def test_split_increment_retained_in_correction_below_reference_ulp(self):
        self.doc["contacts"][0]["bias"] = 1.0
        for row in self.rows[:3]:
            row.update(phin="1", electron_qf_reference_V="1", electron_qf_increment_V="1e-17",
                       phip="1", hole_qf_reference_V="1", hole_qf_increment_V="0")
        self.doc["contacts"][1]["bias"] = 1.0
        for row in self.rows[3:6]:
            row.update(phin="1", electron_qf_reference_V="1", electron_qf_increment_V="0",
                       phip="1", hole_qf_reference_V="1", hole_qf_increment_V="0")
        result = self.run_helper()
        self.assertEqual(result["max_qf_change_V"], 1e-17)

    def test_hdf5_without_split_arrays(self):
        self.fields = self.fields[:7]
        self.rows = [{key: row[key] for key in self.fields} for row in self.rows]
        self.run_helper()
        self.assertEqual(float(self.read_output()[0]["phin"]), 0.1)

    def assert_rejected(self):
        self.output.write_text("existing candidate")
        with self.assertRaises(ValueError):
            self.run_helper()
        self.assertEqual(self.output.read_text(), "existing candidate")

    def test_rejects_unsupported_physics(self):
        original = copy.deepcopy(self.doc)
        variants = [{"carrier_statistics": "fermi_dirac"}, {"statistics": "fermi_dirac"},
                    {"temperature_K": 301}, {"temperature_K": float("nan")},
                    {"recombination_materials": []}, {"recombination_materials": ["PolySilicon"]},
                    {"band_to_band": {"model": "hurkx"}}, {"generation": {"rate": 1}},
                    {"electron_quantum_potential": True}, {"mobility": "constant"},
                    {"impact_ionization": "van_overstraeten"}]
        for change in variants:
            with self.subTest(change=change):
                self.doc = copy.deepcopy(original)
                self.doc["solver"].update(change)
                self.assert_rejected()

    def test_avalanche_requires_canonical_gradqf_without_fallback(self):
        av = {"model": "van_overstraeten", "generation": "current_density",
              "driving_force": "quasi_fermi_gradient", "contact_electric_field_fallback": False,
              "current_approximation": "element_edge_sg_gss_laux",
              "source_mapping_mode": "element_vertex_box_measure",
              "cell_reconstructed_midpoint_density": "bernoulli", "source_jacobian": "local_ad"}
        self.doc["solver"]["impact_ionization"] = av
        self.run_helper()
        for key, value in [("contact_electric_field_fallback", True), ("driving_force", "electric_field"),
                           ("generation", "carrier_density"), ("current_approximation", "cell_reconstructed"),
                           ("driving_force_interpolation", "electric_field"),
                           ("quasi_fermi_carrier_truncation", 1),
                           ("quasi_fermi_carrier_trucation", 1)]:
            with self.subTest(key=key):
                self.doc["solver"]["impact_ionization"] = dict(av, **{key: value})
                self.assert_rejected()

    def test_rejects_density_field_avalanche_interpolation(self):
        self.doc["solver"]["impact_ionization"] = {
            "model": "van_overstraeten", "generation": "current_density",
            "driving_force": "quasi_fermi_gradient", "contact_electric_field_fallback": False,
            "current_approximation": "element_edge_sg_gss_laux",
            "source_mapping_mode": "element_vertex_box_measure",
            "cell_reconstructed_midpoint_density": "bernoulli", "source_jacobian": "local_ad",
            "driving_force_interpolation": "none", "electron_driving_force_ref_density_m3": 1e18}
        self.assert_rejected()

    def test_rejects_zero_mobility_partial_split_and_bad_mesh(self):
        self.doc["solver"]["mobility"]["ialmob"]["constant_material_mobility_cm2_per_Vs"]["PolySilicon"]["hole"] = 0
        self.assert_rejected()
        self.doc["solver"]["mobility"]["ialmob"]["constant_material_mobility_cm2_per_Vs"]["PolySilicon"]["hole"] = 470.5
        self.fields.remove("hole_qf_reference_V")
        for row in self.rows:
            del row["hole_qf_reference_V"]
        self.assert_rejected()
        self.fields = self.fields[:7]
        self.rows = [{key: row[key] for key in self.fields} for row in self.rows]
        self.mesh["triangles"][0]["node_ids"] = [0, 1, 100]
        self.assert_rejected()

    def test_allows_multiple_equal_bias_contacts_on_component(self):
        self.doc["contacts"].append({"name": "gate2", "type": "ohmic", "bias": 0.1})
        self.mesh["contacts"].append({"name": "gate2", "node_ids": [1]})
        result = self.run_helper()
        self.assertEqual(result["components"][0]["contacts"], ["gate", "gate2"])

    def test_rejects_external_circuit(self):
        self.doc["sweep"]["external_circuit"] = {"enabled": True}
        self.assert_rejected()

    def test_rejects_non300_poly_material(self):
        (self.root / "materials.json").write_text(json.dumps({"materials": [
            {"name": "PolySilicon", "temperature_K": 350, "ni": 1e16, "mun": 0.1, "mup": 0.04}]}))
        self.assert_rejected()

    def test_rejects_floating_conflicting_and_nonohmic_contacts(self):
        original = copy.deepcopy(self.mesh)
        for contacts in [original["contacts"][1:], original["contacts"] + [{"name": "SLOT", "node_ids": [1]}]]:
            self.mesh = dict(original, contacts=contacts)
            self.assert_rejected()
        self.mesh = original
        self.doc["contacts"][0]["type"] = "schottky"
        self.assert_rejected()

    def test_rejects_poly_shared_with_semiconductor_or_unknown_material(self):
        self.mesh["triangles"][2]["node_ids"] = [2, 7, 8]
        self.assert_rejected()
        self.mesh["regions"][1]["material"] = "UnknownTransport"
        self.assert_rejected()

    def test_versioned_material_transport_keys_remove_insulating_exemption(self):
        # The oxide triangle shares nodes 2 and 5 with the two poly components.
        keys = ("intrinsic_carrier_density_cm3", "electron_mobility_cm2_per_V_s",
                "hole_mobility_cm2_per_V_s")
        oxide = {"name": "SiO2", "role": "dielectric", "eps_r": 3.9,
                 "provenance": "synthetic material guard regression", **dict.fromkeys(keys, 0)}
        materials = {
            "schema": "vela.templates_ldmos.materials.v1", "benchmark": "synthetic", "revision": 1,
            "unit_system": {"concentration": "cm^-3", "mobility": "cm^2/(V*s)", "energy": "eV",
                            "temperature": "K", "thermal_conductivity": "W/(m*K)",
                            "specific_heat": "J/(kg*K)", "mass_density": "kg/m^3"},
            "materials": [oxide, {"name": "PolySilicon", "role": "transport_semiconductor",
                                   "eps_r": 11.7, "provenance": "synthetic material guard regression",
                                   keys[0]: 1e10, keys[1]: 1417, keys[2]: 470.5}],
        }
        path = self.root / "materials.json"
        path.write_text(json.dumps(materials))
        self.run_helper()  # An explicitly zero-transport oxide remains supported.
        for key in keys:
            for value in (1.0, float("nan")):
                with self.subTest(key=key, value=value):
                    oxide.update(dict.fromkeys(keys, 0))
                    oxide[key] = value
                    path.write_text(json.dumps(materials))
                    self.assert_rejected()

    def test_rejects_malformed_nonfinite_or_incomplete_state(self):
        original = copy.deepcopy(self.rows)
        for field, value in [("psi", "nan"), ("electrons_m3", "-1"), ("holes_m3", "0"),
                             ("phin", "inf"), ("phin", "garbage"), ("node_id", "9"),
                             ("node_id", "1"), ("electron_qf_reference_V", "3")]:
            with self.subTest(field=field, value=value):
                self.rows = copy.deepcopy(original)
                self.rows[0][field] = value
                self.assert_rejected()
        self.rows = original[:-1]
        self.assert_rejected()

    def test_rejects_density_overflow_before_output(self):
        self.rows[0].update(phin="100", electron_qf_reference_V="100", electron_qf_increment_V="0")
        self.assert_rejected()

    def test_rejects_overwriting_input(self):
        self.persist()
        original = self.source.read_bytes()
        with self.assertRaises(ValueError):
            self.write_restart(self.doc, self.source, self.source, self.root, 2)
        self.assertEqual(self.source.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
