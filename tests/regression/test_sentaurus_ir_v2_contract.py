"""Schema and loss-aware migration tests for M0 mesh contracts."""

from __future__ import annotations

import copy
import json
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from sentaurus_device_ir import parse_sde_device_ir  # noqa: E402
from sentaurus_ir_v2 import (  # noqa: E402
    DeviceIrMigrationError,
    migrate_device_ir_v1,
    require_known_provenance,
)

PN2D_SDE = (
    REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source" / "pn2d_sde.cmd"
)
IR_V2_SCHEMA = REPO / "schemas" / "vela.sentaurus_device_ir.v2.schema.json"
BUNDLE_SCHEMA = REPO / "schemas" / "vela.mesh_field_bundle.v1.schema.json"


def validate_schema(testcase: unittest.TestCase, document: dict, path: Path) -> None:
    schema = json.loads(path.read_text(encoding="utf-8"))
    try:
        import jsonschema  # type: ignore[import-not-found]
    except ImportError:
        for key in schema["required"]:
            testcase.assertIn(key, document)
        testcase.assertEqual(set(document), set(schema["properties"]))
        return
    jsonschema.validate(document, schema)


class DeviceIrV2MigrationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.migrated = migrate_device_ir_v1(parse_sde_device_ir(PN2D_SDE))

    def test_pn2d_v1_migrates_to_schema_valid_v2(self) -> None:
        validate_schema(self, self.migrated, IR_V2_SCHEMA)
        self.assertEqual("vela.sentaurus_device_ir.v2", self.migrated["schema"])
        body = self.migrated["geometry"]["bodies"][0]
        self.assertEqual("R.Si", body["region"])
        self.assertEqual({"sde": "Silicon", "vela": "Si"}, body["material"])
        self.assertEqual("ABA", body["boolean_mode"])
        self.assertEqual(0, body["creation_order"])

    def test_current_parser_preserves_geometry_and_contact_provenance(self) -> None:
        self.assertEqual("known", self.migrated["geometry"]["bodies"][0]["provenance"]["status"])
        self.assertEqual("known", self.migrated["contacts"][0]["provenance"]["status"])
        self.assertEqual(
            "known", self.migrated["contacts"][0]["assignments"][0]["provenance"]["status"])

    def test_legacy_v1_without_locations_is_explicit_and_fail_closed(self) -> None:
        legacy = copy.deepcopy(parse_sde_device_ir(PN2D_SDE))
        legacy["geometry"]["regions"][0].pop("location")
        migrated = migrate_device_ir_v1(legacy)
        manifest = migrated["capability_manifest"]
        self.assertFalse(manifest["provenance_complete"])
        required = ["/geometry/bodies/0/provenance"]
        with self.assertRaises(DeviceIrMigrationError) as context:
            require_known_provenance(migrated, required)
        self.assertIn(required[0], str(context.exception))

    def test_v1_contact_points_become_typed_assignments(self) -> None:
        contacts = {item["name"]: item for item in self.migrated["contacts"]}
        self.assertEqual([0.0, 0.25], contacts["Anode"]["assignments"][0]["position"])
        self.assertEqual("pick_point", contacts["Anode"]["assignments"][0]["kind"])


class MeshFieldBundleSchemaTest(unittest.TestCase):
    def test_minimal_bundle_contract_validates(self) -> None:
        field = {"association": "node", "unit": "cm^-3", "values": [1.0, 2.0, 3.0]}
        log_field = {
            "association": "node", "unit": "log10(cm^-3)",
            "values": [0.0, 0.3010299956639812, 0.47712125471966244],
        }
        document = {
            "schema": "vela.mesh_field_bundle.v1",
            "coordinate_unit": "um",
            "field_units": {"concentration": "cm^-3"},
            "provenance": {
                "device_ir_schema": "vela.sentaurus_device_ir.v2",
                "mesh_backend": "gmsh",
            },
            "gmsh": {"version": "4.15.2", "msh_version": "4.1", "options": {}},
            "nodes": [
                {"gmsh_tag": index + 10, "vela_id": index, "dfise_id": index,
                 "x": float(index == 1), "y": float(index == 2), "z": 0.0}
                for index in range(3)
            ],
            "triangles": [{
                "gmsh_tag": 20, "vela_id": 0, "dfise_id": 0,
                "node_ids": [0, 1, 2], "region_id": 0, "material_id": 0,
                "physical_surface_tags": [1001, 2001],
            }],
            "lines": [{
                "gmsh_tag": 30, "vela_id": 0, "dfise_id": 0,
                "node_ids": [0, 1], "kind": "contact",
                "physical_curve_tags": [3001],
            }],
            "regions": [{"id": 0, "name": "R.Si", "physical_name": "region:R.Si", "physical_tag": 1001}],
            "materials": [{"id": 0, "name": "Si", "physical_name": "material:Si", "physical_tag": 2001}],
            "contacts": [{"id": 0, "name": "Anode", "physical_name": "contact:Anode", "physical_tag": 3001}],
            "interfaces": [],
            "fields": {
                "DonorConcentration": field,
                "AcceptorConcentration": field,
                "DopingConcentration": field,
                "AbsDopingConcentration": field,
                "Log10AbsDopingConcentration": log_field,
                "species": {"BoronActiveConcentration": field},
            },
            "id_mapping": {
                "gmsh_node_to_vela": {"10": 0, "11": 1, "12": 2},
                "gmsh_triangle_to_vela": {"20": 0},
                "gmsh_line_to_vela": {"30": 0},
                "vela_node_to_dfise": {"0": 0, "1": 1, "2": 2},
                "vela_triangle_to_dfise": {"0": 0},
                "vela_line_to_dfise": {"0": 0},
            },
        }
        validate_schema(self, document, BUNDLE_SCHEMA)


if __name__ == "__main__":
    unittest.main()
