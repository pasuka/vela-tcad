from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO / "scripts" / "prepare_slot_ldmos_bvds.py"
SPEC = importlib.util.spec_from_file_location("prepare_slot_ldmos_bvds", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


class SlotLdmosPreparationTests(unittest.TestCase):
    def make_export(self, root: Path, *, obtuse: bool = False) -> Path:
        export = root / "export"
        export.mkdir()
        # Coordinates emulate the centimetre-valued SProcess neutral export.
        oxide_tip = [3.5e-4, 1.0e-6] if obtuse else [4.0e-4, 1.0e-4]
        nodes = [
            [0, 0.0, 0.0], [1, 1.0e-4, 0.0], [2, 0.0, 1.0e-4],
            [3, 3.0e-4, 0.0], [4, 4.0e-4, 0.0], [5, *oxide_tip],
            [6, 5.0e-4, 0.0], [7, 5.0e-4, 1.0e-4],
        ]
        write_csv(export / "nodes.csv", ["id", "x_um", "y_um"], nodes)
        write_csv(
            export / "elements.csv",
            ["id", "node0", "node1", "node2", "region", "material"],
            [
                [0, 0, 1, 2, "Silicon_1", "Si"],
                [1, 3, 4, 5, "Oxide_1", "SiO2"],
                [2, 4, 6, 7, "PolySilicon_1", "PolySilicon"],
            ],
        )
        write_csv(
            export / "contacts.csv",
            ["name", "node_ids", "region"],
            [
                ["source", "0;2", "Silicon_1"],
                ["drain", "1;2", "Silicon_1"],
                ["gate", "6;7", "PolySilicon_1"],
                ["SLOT", "3;5", "Oxide_1"],
                ["substrate", "0;1", "Silicon_1"],
            ],
        )
        write_csv(
            export / "doping.csv",
            ["node_id", "donors_cm3", "acceptors_cm3"],
            [[node[0], 1.0e15, 0.0] for node in nodes],
        )
        return export

    def test_legacy_profile_seals_slot_and_coordinate_policy(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_slot_ldmos_") as tmp:
            root = Path(tmp)
            export = self.make_export(root)
            output = root / "output"
            report = MODULE.prepare(
                export, output, "legacy_cell_reconstructed", 60.0, 1.0e4
            )

            mesh = json.loads((output / "mesh.json").read_text(encoding="utf-8"))
            boundary = json.loads(
                (output / "simulation_slot_boundary_check.json").read_text(encoding="utf-8")
            )
            bv = json.loads(
                (output / "simulation_bvds_legacy.json").read_text(encoding="utf-8")
            )
            manifest = json.loads(
                (output / "slot_ldmos_bvds_stages_manifest.json").read_text(
                    encoding="utf-8"
                )
            )
            materials = json.loads(
                (output / "materials.json").read_text(encoding="utf-8")
            )
            prepared_doping = {
                int(row["node_id"]): row
                for row in MODULE.read_csv(output / "doping.csv")
            }
            contacts = {item["name"]: item for item in boundary["contacts"]}

            self.assertEqual(mesh["nodes"][1]["x"], 1.0)
            self.assertEqual(report["coordinate_conversion"]["scale_to_um"], 1.0e4)
            self.assertEqual(contacts["SLOT"]["type"], "metal_gate")
            self.assertEqual(contacts["gate"]["type"], "metal_gate")
            self.assertEqual(contacts["drain"]["type"], "ohmic")
            self.assertFalse(report["contact_policy"]["SLOT"]["carrier_dirichlet"])
            self.assertEqual(
                set(report["contact_policy"]["SLOT"]["incident_material_cell_counts"]),
                {"SiO2"},
            )
            self.assertEqual(
                bv["mesh_geometry"],
                {"node_volume_policy": "barycentric", "require_non_obtuse": False},
            )
            self.assertEqual(
                bv["solver"]["impact_ionization"]["current_approximation"],
                "cell_reconstructed",
            )
            self.assertEqual(
                bv["sweep"]["external_circuit"]["resistance_ohm_um"], 1.0e12
            )
            self.assertEqual(
                manifest["execution_order"],
                [
                    "00_equilibrium",
                    "01_unit_resistor_1v",
                    "02_avalanche_off_60v",
                    "03_iic_postprocess_60v",
                    "04_avalanche_activation_1v",
                    "05_avalanche_on_60v",
                    "06_bvds_external_resistor_final",
                ],
            )
            self.assertEqual(
                manifest["external_circuit_contract"]["equation"],
                "Vouter_V = Vinner_V + resistance_ohm_um * Id_A_per_um",
            )
            poly = next(
                item for item in materials["materials"]
                if item["name"] == "PolySilicon"
            )
            for transport_key in ("ni", "mun", "mup", "Nc_m3", "Nv_m3"):
                self.assertNotIn(transport_key, poly)
            self.assertEqual(float(prepared_doping[0]["donors_cm3"]), 1.0e15)
            self.assertEqual(float(prepared_doping[2]["donors_cm3"]), 1.0e15)
            for node_id in range(3, 8):
                self.assertEqual(float(prepared_doping[node_id]["donors_cm3"]), 0.0)
                self.assertEqual(float(prepared_doping[node_id]["acceptors_cm3"]), 0.0)
            self.assertEqual(
                report["nontransport_region_policy"]["doping"][
                    "zeroed_node_count"
                ],
                5,
            )

    def test_staged_configs_separate_physics_and_chain_states(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_slot_ldmos_stages_") as tmp:
            root = Path(tmp)
            output = root / "output"
            MODULE.prepare(
                self.make_export(root),
                output,
                "legacy_cell_reconstructed",
                60.0,
                1.0e4,
            )

            def load(name: str) -> dict[str, object]:
                return json.loads((output / name).read_text(encoding="utf-8"))

            equilibrium = load("simulation_00_equilibrium.json")
            unit = load("simulation_01_unit_resistor_1v.json")
            avalanche_off = load("simulation_02_avalanche_off_60v.json")
            iic = load("simulation_03_iic_postprocess_60v.json")
            avalanche_activation = load(
                "simulation_04_avalanche_activation_1v.json"
            )
            avalanche_on = load("simulation_05_avalanche_on_60v.json")
            final = load("simulation_06_bvds_external_resistor_final.json")

            self.assertTrue(
                (output / "outputs/stages/00_equilibrium/states").is_dir()
            )
            self.assertTrue(
                (
                    output
                    / "outputs/stages/01_unit_resistor_1v/boundary_control_checkpoints"
                ).is_dir()
            )
            self.assertNotIn("impact_ionization", equilibrium["solver"])
            self.assertNotIn("external_circuit", equilibrium["sweep"])
            self.assertEqual(equilibrium["sweep"]["bias_points"], [0.0])

            for document in [
                unit,
                avalanche_off,
                iic,
                avalanche_activation,
                avalanche_on,
                final,
            ]:
                circuit = document["sweep"]["external_circuit"]
                self.assertEqual(circuit["mode"], "series_resistor")
                self.assertEqual(circuit["resistance_ohm_um"], 1.0e12)
                self.assertTrue(document["sweep"]["boundary_control"]["resume"])

            self.assertEqual(unit["sweep"]["bias_points"], [0.0, 0.1, 0.2, 0.5, 1.0])
            self.assertEqual(
                unit["sweep"]["external_circuit"]["max_inner_voltage_step_V"],
                0.01,
            )
            for document in [avalanche_off, iic, avalanche_on, final]:
                self.assertEqual(
                    document["sweep"]["external_circuit"][
                        "max_inner_voltage_step_V"
                    ],
                    0.01,
                )
                self.assertEqual(
                    document["sweep"]["boundary_control"][
                        "predictor_max_step_factor"
                    ],
                    2.0,
                )
            self.assertEqual(
                avalanche_activation["sweep"]["external_circuit"][
                    "max_inner_voltage_step_V"
                ],
                1.0e-4,
            )
            self.assertEqual(
                avalanche_activation["sweep"]["boundary_control"][
                    "predictor_max_step_factor"
                ],
                2.0,
            )
            self.assertEqual(
                unit["sweep"]["boundary_control"]["predictor_max_step_factor"],
                2.0,
            )
            self.assertNotIn("handoff", equilibrium["solver"])
            self.assertEqual(unit["solver"]["handoff"]["gummel_max_iter"], 0)
            self.assertFalse(
                unit["solver"]["handoff"]["require_gummel_convergence"]
            )
            self.assertEqual(
                unit["solver"]["residual_scales"],
                MODULE.SLOT_LDMOS_RESTART_RESIDUAL_SCALES,
            )

            self.assertNotIn("impact_ionization", avalanche_off["solver"])
            self.assertNotIn("diagnostics", avalanche_off["sweep"])
            self.assertEqual(
                iic["solver"]["impact_ionization"]["coupling_mode"],
                "postprocess_only",
            )
            for document in [avalanche_activation, avalanche_on]:
                self.assertEqual(
                    document["solver"]["impact_ionization"]["coupling_mode"],
                    "self_consistent",
                )
            self.assertEqual(
                avalanche_on["solver"]["carrier_regularization_scale"], 1.0e-8
            )
            self.assertEqual(
                avalanche_on["solver"]["handoff"]["gummel_max_iter"], 50
            )
            self.assertEqual(
                avalanche_on["solver"]["impact_ionization"]["source_jacobian"],
                "local_ad",
            )
            self.assertEqual(avalanche_activation["sweep"]["start"], 1.0)
            self.assertEqual(
                avalanche_activation["sweep"]["initial_state_file"],
                "outputs/stages/01_unit_resistor_1v/final_state.h5",
            )
            self.assertEqual(
                avalanche_activation["sweep"]["external_circuit"][
                    "initial_inner_voltage_V"
                ],
                MODULE.VELA_UNIT_RESISTOR_INNER_V,
            )
            self.assertEqual(avalanche_on["sweep"]["start"], 10.0)
            self.assertEqual(
                avalanche_on["sweep"]["initial_state_file"],
                "outputs/stages/04_avalanche_activation_1v/final_state.h5",
            )
            for document in [iic, avalanche_activation, avalanche_on, final]:
                audit = document["sweep"]["diagnostics"][
                    "release_bv_config_audit"
                ]
                self.assertTrue(audit["enabled"])
                self.assertTrue(audit["csv_file"].endswith("avalanche_summary.csv"))
            self.assertEqual(
                final["sweep"]["initial_state_file"],
                "outputs/stages/05_avalanche_on_60v/final_state.h5",
            )
            threshold_outer = (
                MODULE.SENTAURUS_REFERENCE_BVDS_V
                + MODULE.SENTAURUS_SERIES_RESISTANCE_OHM_UM
                * MODULE.SENTAURUS_BREAKDOWN_CURRENT_A_PER_UM
            )
            self.assertIn(threshold_outer, final["sweep"]["bias_points"])
            self.assertGreater(final["sweep"]["bias_points"][-1], threshold_outer)

    def test_sg_laux_profile_fails_closed_on_obtuse_mesh(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_slot_ldmos_obtuse_") as tmp:
            root = Path(tmp)
            export = self.make_export(root, obtuse=True)
            with self.assertRaisesRegex(
                MODULE.PreparationError, "element_edge_sg_gss_laux is forbidden"
            ):
                MODULE.prepare(
                    export,
                    root / "output",
                    "element_edge_sg_gss_laux",
                    60.0,
                    1.0e4,
                )


if __name__ == "__main__":
    unittest.main()
