#!/usr/bin/env python3
"""Versioned Sentaurus Device IR v2 migration helpers.

The v1 parser did not retain source locations for every semantic record.  The
migration therefore preserves all known values and marks unrecoverable
provenance as ``unknown`` instead of inventing source evidence.
"""

from __future__ import annotations

import copy
from typing import Any, Iterable


DEVICE_IR_V1 = "vela.sentaurus_device_ir.v1"
DEVICE_IR_V2 = "vela.sentaurus_device_ir.v2"


class DeviceIrMigrationError(ValueError):
    """Raised when an IR cannot be migrated without guessing semantics."""


def _unknown(reason: str) -> dict[str, Any]:
    return {
        "status": "unknown",
        "source_location": None,
        "reason": reason,
    }


def _known(location: str) -> dict[str, Any]:
    return {
        "status": "known",
        "source_location": location,
        "reason": None,
    }


def _copy_with_provenance(
    records: Iterable[dict[str, Any]],
    unknown_paths: list[str],
    prefix: str,
) -> list[dict[str, Any]]:
    migrated: list[dict[str, Any]] = []
    for index, record in enumerate(records):
        item = copy.deepcopy(record)
        location = item.get("location")
        if isinstance(location, str) and location:
            item["provenance"] = _known(location)
        else:
            item["provenance"] = _unknown(
                "source location was not retained by Device IR v1")
            unknown_paths.append(f"{prefix}/{index}/provenance")
        migrated.append(item)
    return migrated


def migrate_device_ir_v1(document: dict[str, Any]) -> dict[str, Any]:
    """Migrate a v1 Device IR to the explicit v2 skeleton.

    The conversion is deliberately loss-aware.  Geometry and contact source
    locations cannot be recovered from v1 and are reported in the capability
    manifest so a provenance-sensitive Gmsh operation can fail closed.
    """
    if document.get("schema") != DEVICE_IR_V1:
        raise DeviceIrMigrationError(
            f"expected {DEVICE_IR_V1}, got {document.get('schema')!r}")

    materials = {
        item["region"]: {
            "sde": item["sde_material"],
            "vela": item["vela_material"],
        }
        for item in document.get("materials", [])
    }
    unknown_paths: list[str] = []
    bodies: list[dict[str, Any]] = []
    for index, region in enumerate(document["geometry"]["regions"]):
        name = region["name"]
        if name not in materials:
            raise DeviceIrMigrationError(
                f"v1 region {name!r} has no material assignment")
        location = region.get("location")
        provenance = _known(location) if isinstance(location, str) and location else _unknown(
            "body source location was not retained by Device IR v1")
        bodies.append({
            "id": f"body:{region['order']}",
            "region": name,
            "material": materials[name],
            "shape": region["shape"],
            "geometry": {
                "lower_left": copy.deepcopy(region["lower_left"]),
                "upper_right": copy.deepcopy(region["upper_right"]),
            },
            "creation_order": region["order"],
            "boolean_mode": region["boolean_mode"],
            "provenance": provenance,
        })
        if provenance["status"] == "unknown":
            unknown_paths.append(f"/geometry/bodies/{index}/provenance")

    inserted_vertices: list[dict[str, Any]] = []
    for index, vertex in enumerate(document["geometry"]["inserted_vertices"]):
        location = vertex.get("location")
        provenance = _known(location) if isinstance(location, str) and location else _unknown(
            "insert-vertex source location was not retained by Device IR v1")
        inserted_vertices.append({
            "position": [vertex["x"], vertex["y"]],
            "provenance": provenance,
        })
        if provenance["status"] == "unknown":
            unknown_paths.append(
                f"/geometry/inserted_vertices/{index}/provenance")

    contacts: list[dict[str, Any]] = []
    for contact_index, contact in enumerate(document["contacts"]):
        assignments = []
        for pick_index, point in enumerate(contact["pick_points"]):
            location = point.get("location")
            provenance = _known(location) if isinstance(location, str) and location else _unknown(
                "contact assignment source location was not retained by Device IR v1")
            assignments.append({
                "kind": "pick_point",
                "position": [point["x"], point["y"]],
                "entity_ids": [],
                "provenance": provenance,
            })
            if provenance["status"] == "unknown":
                unknown_paths.append(
                    f"/contacts/{contact_index}/assignments/{pick_index}/provenance")
        contact_location = contact.get("location")
        contact_provenance = (
            _known(contact_location)
            if isinstance(contact_location, str) and contact_location
            else _unknown("contact-set source location was not retained by Device IR v1")
        )
        contacts.append({
            "name": contact["name"],
            "display": copy.deepcopy(contact["display"]),
            "current_state": "active" if assignments else "inactive",
            "assignments": assignments,
            "provenance": contact_provenance,
        })
        if contact_provenance["status"] == "unknown":
            unknown_paths.append(f"/contacts/{contact_index}/provenance")

    doping = copy.deepcopy(document["doping"])
    doping["profiles"] = _copy_with_provenance(
        doping["profiles"], unknown_paths, "/doping/profiles")
    doping["placements"] = _copy_with_provenance(
        doping["placements"], unknown_paths, "/doping/placements")
    doping["windows"] = _copy_with_provenance(
        doping["windows"], unknown_paths, "/doping/windows")

    mesh_control = copy.deepcopy(document["mesh_control"])
    for key in ("refinement_sizes", "refinement_placements",
                "refinement_functions", "windows"):
        mesh_control[key] = _copy_with_provenance(
            mesh_control[key], unknown_paths, f"/mesh_control/{key}")

    return {
        "schema": DEVICE_IR_V2,
        "source": {
            "path": document["source"],
            "frontend": "unknown",
            "input_schema": DEVICE_IR_V1,
        },
        "defines": copy.deepcopy(document["defines"]),
        "geometry": {
            "bodies": bodies,
            "inserted_vertices": inserted_vertices,
        },
        "contacts": contacts,
        "doping": doping,
        "mesh_control": mesh_control,
        "capability_manifest": {
            "migration": "v1-to-v2",
            "provenance_complete": not unknown_paths,
            "unknown_provenance_paths": unknown_paths,
        },
        "metadata_only": copy.deepcopy(document["metadata_only"]),
    }


def require_known_provenance(
    document: dict[str, Any], required_paths: Iterable[str]
) -> None:
    """Fail closed if an operation requires provenance absent after migration."""
    unknown = set(document.get("capability_manifest", {}).get(
        "unknown_provenance_paths", []))
    blocked = sorted(set(required_paths) & unknown)
    if blocked:
        raise DeviceIrMigrationError(
            "operation requires provenance unavailable in migrated Device IR: "
            + ", ".join(blocked))
