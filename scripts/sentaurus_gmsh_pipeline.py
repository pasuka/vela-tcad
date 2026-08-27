#!/usr/bin/env python3
"""Public offline SDE -> Gmsh mesh and cross-format export pipeline."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Sequence

from compare_sentaurus_tdr_tdx import parse_dat, parse_grd
from gmsh_runtime import (
    EXPECTED_GMSH_VERSION,
    GmshRuntimeError,
    _load_staged_module,
    discover_runtime_sources,
    executable_version,
    sha256,
    stage_runtime,
)
from sentaurus_device_ir import SdeParseError, parse_sde_device_ir
from sentaurus_gmsh_backend import (
    GmshBackendError,
    build_m1_bundle,
    build_m2_bundle,
    build_m3_bundle,
    build_m5_bundle,
    read_msh_roundtrip,
    write_msh_with_node_data,
)
from sentaurus_ir_v2 import migrate_device_ir_v1
from sentaurus_mesh_exporters import (
    POINT_FIELD_NAMES,
    bundle_to_vela_mesh,
    compare_nodal_fields,
    read_doping_csv,
    read_legacy_vtk,
    read_vtu,
    write_dfise,
    write_doping_csv,
    write_json,
    write_legacy_vtk,
    write_vtu,
)


REPO = Path(__file__).resolve().parents[1]
RUN_SCHEMA = "vela.sentaurus_gmsh_preprocess.v1"
FAIL_SCHEMA = "vela.sentaurus_gmsh_preprocess.fail_report.v1"


class PipelineError(RuntimeError):
    def __init__(self, message: str, stage: str, details: dict[str, Any] | None = None):
        super().__init__(message)
        self.stage = stage
        self.details = details or {}


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _strict_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _load_runtime(args: argparse.Namespace) -> tuple[Any, dict[str, Any]]:
    sources = discover_runtime_sources(
        args.gmsh_python, args.gmsh_library, args.gmsh_executable)
    cli_version = executable_version(sources.executable)
    versions = {
        "expected": EXPECTED_GMSH_VERSION,
        "wrapper": sources.version,
        "executable": cli_version,
    }
    if len(set(versions.values())) != 1:
        raise PipelineError(
            "Gmsh wrapper, executable, and pinned versions do not match",
            "gmsh_runtime", {"versions": versions})
    staged = stage_runtime(sources, args.runtime_dir)
    try:
        gmsh = _load_staged_module(Path(staged["wrapper"]))
    except Exception as error:  # noqa: BLE001
        raise PipelineError(
            f"failed to import staged Gmsh Python API: {error}",
            "gmsh_runtime", {"staging": staged}) from error
    api_version = str(getattr(gmsh, "__version__", "unknown"))
    loaded = Path(str(getattr(gmsh, "libpath", ""))).resolve()
    if api_version != EXPECTED_GMSH_VERSION or loaded != Path(staged["library"]).resolve():
        raise PipelineError(
            "loaded Gmsh API or shared library does not match staged runtime",
            "gmsh_runtime",
            {"api_version": api_version, "loaded_library": str(loaded), "staging": staged},
        )
    manifest = {
        "versions": {**versions, "api": api_version},
        "python": {"executable": str(Path(sys.executable).resolve()),
                   "version": platform.python_version()},
        "sources": {
            "wrapper": str(sources.wrapper), "wrapper_sha256": sha256(sources.wrapper),
            "library": str(sources.library), "library_sha256": sha256(sources.library),
            "executable": str(sources.executable), "executable_sha256": sha256(sources.executable),
        },
        "staging": {**staged, "loaded_library": str(loaded)},
    }
    return gmsh, manifest


def _msh_fields(bundle: dict[str, Any], readback: dict[str, Any]) -> dict[str, list[float]]:
    result: dict[str, list[float]] = {}
    for name in POINT_FIELD_NAMES:
        values = readback["fields"].get(name)
        if values is None:
            continue
        result[name] = [values[item["gmsh_tag"]] for item in bundle["nodes"]]
    return result


def _csv_consistency(bundle: dict[str, Any], csv_fields: dict[str, list[float]]) -> list[str]:
    failures = []
    for name in ("DonorConcentration", "AcceptorConcentration", "DopingConcentration"):
        expected = bundle["fields"][name]["values"]
        actual = csv_fields.get(name, [])
        if len(actual) != len(expected) or any(a != b for a, b in zip(actual, expected)):
            failures.append(f"doping.csv {name} mismatch")
    return failures


def _expected_vtk(bundle: dict[str, Any]) -> dict[str, Any]:
    cells = [item["node_ids"] for item in bundle["triangles"]]
    cells.extend(item["node_ids"] for item in bundle["lines"])
    types = [5] * len(bundle["triangles"]) + [3] * len(bundle["lines"])
    metadata = {
        name: [] for name in (
            "RegionId", "MaterialId", "PhysicalSurfaceTag", "ElementId",
            "ContactId", "PhysicalCurveTag",
        )
    }
    for triangle in bundle["triangles"]:
        values = {
            "RegionId": triangle["region_id"],
            "MaterialId": triangle["material_id"],
            "PhysicalSurfaceTag": triangle["physical_surface_tags"][0],
            "ElementId": triangle["vela_id"],
            "ContactId": -1,
            "PhysicalCurveTag": -1,
        }
        for name, value in values.items():
            metadata[name].append(value)
    contact_id = {item["physical_tag"]: item["id"] for item in bundle["contacts"]}
    for line in bundle["lines"]:
        physical_tag = line["physical_curve_tags"][0]
        values = {
            "RegionId": -1,
            "MaterialId": -1,
            "PhysicalSurfaceTag": -1,
            "ElementId": line["vela_id"],
            "ContactId": contact_id.get(physical_tag, -1),
            "PhysicalCurveTag": physical_tag,
        }
        for name, value in values.items():
            metadata[name].append(value)
    return {
        "points": [[item["x"], item["y"], item["z"]] for item in bundle["nodes"]],
        "cells": cells,
        "cell_types": types,
        "cell_fields": metadata,
    }


def _compare_vtk_payload(
    label: str, actual: dict[str, Any], expected: dict[str, Any], failures: list[str]
) -> None:
    if actual["points"] != expected["points"]:
        failures.append(f"{label} point coordinates mismatch")
    if actual["cells"] != expected["cells"]:
        failures.append(f"{label} cell connectivity/order mismatch")
    if actual["cell_types"] != expected["cell_types"]:
        failures.append(f"{label} Tri3/Line2 cell types mismatch")
    if actual["cell_fields"] != expected["cell_fields"]:
        failures.append(f"{label} CellData ownership mismatch")


def _compare_msh_topology(
    bundle: dict[str, Any], msh: dict[str, Any], failures: list[str]
) -> None:
    expected_nodes = {
        item["gmsh_tag"]: [item["x"], item["y"], item["z"]]
        for item in bundle["nodes"]
    }
    actual_nodes = msh["node_coordinates"]
    if set(actual_nodes) != set(expected_nodes) or any(
        abs(actual - expected) > 1.0e-14 * max(1.0, abs(expected))
        for tag, expected_row in expected_nodes.items()
        for actual, expected in zip(actual_nodes.get(tag, []), expected_row)
    ):
        failures.append("MSH node tags/coordinates mismatch")
    expected_triangles = {
        item["gmsh_tag"]: [bundle["nodes"][node]["gmsh_tag"] for node in item["node_ids"]]
        for item in bundle["triangles"]
    }
    expected_lines = {
        item["gmsh_tag"]: [bundle["nodes"][node]["gmsh_tag"] for node in item["node_ids"]]
        for item in bundle["lines"]
    }
    if msh["elements"]["triangles"] != expected_triangles:
        failures.append("MSH Tri3 tags/connectivity mismatch")
    if msh["elements"]["lines"] != expected_lines:
        failures.append("MSH Line2 tags/connectivity mismatch")


def _decode_dfise_triangle(row: list[float], edges: list[list[float]]) -> list[int] | None:
    values = [int(value) for value in row]
    if len(values) != 4 or values[0] != 2:
        return None
    oriented = []
    for reference in values[1:]:
        index = reference if reference >= 0 else -reference - 1
        if index < 0 or index >= len(edges):
            return None
        edge = [int(value) for value in edges[index]]
        oriented.append(edge if reference >= 0 else list(reversed(edge)))
    if any(oriented[index][1] != oriented[(index + 1) % 3][0] for index in range(3)):
        return None
    return [edge[0] for edge in oriented]


def _compare_dfise(bundle: dict[str, Any], grd: dict[str, Any], failures: list[str]) -> None:
    expected_vertices = [[item["x"], item["y"]] for item in bundle["nodes"]]
    if len(grd["vertices"]) != len(expected_vertices) or any(
        abs(actual - expected) > 1.0e-13 * max(1.0, abs(expected))
        for actual_row, expected_row in zip(grd["vertices"], expected_vertices)
        for actual, expected in zip(actual_row, expected_row)
    ):
        failures.append("DF-ISE GRD vertex coordinates/order mismatch")
    expected_elements = len(bundle["triangles"]) + len(bundle["lines"])
    if grd["element_count"] != expected_elements or len(grd["elements"]) != expected_elements:
        failures.append("DF-ISE GRD element count mismatch")
        return
    for index, triangle in enumerate(bundle["triangles"]):
        if _decode_dfise_triangle(grd["elements"][index], grd["edges"]) != triangle["node_ids"]:
            failures.append(f"DF-ISE GRD triangle {index} orientation/connectivity mismatch")
            break
    offset = len(bundle["triangles"])
    for index, line in enumerate(bundle["lines"]):
        row = [int(value) for value in grd["elements"][offset + index]]
        if row != [1, *line["node_ids"]]:
            failures.append(f"DF-ISE GRD contact line {index} connectivity mismatch")
            break
    material_name = {item["id"]: item["name"] for item in bundle["materials"]}
    expected_regions: dict[str, list[int]] = {}
    expected_materials = []
    for region in bundle["regions"]:
        triangles = [item for item in bundle["triangles"] if item["region_id"] == region["id"]]
        expected_regions[region["name"]] = [item["dfise_id"] for item in triangles]
        expected_materials.append(material_name[triangles[0]["material_id"]])
    for contact in bundle["contacts"]:
        expected_regions[contact["name"]] = [
            offset + item["dfise_id"] for item in bundle["lines"]
            if contact["physical_tag"] in item["physical_curve_tags"]
        ]
        expected_materials.append("Contact")
    for interface in bundle["interfaces"]:
        expected_regions[interface["name"]] = [
            offset + item["dfise_id"] for item in bundle["lines"]
            if interface["physical_tag"] in item["physical_curve_tags"]
        ]
        expected_materials.append("Interface")
    if grd["region_elements"] != expected_regions:
        failures.append("DF-ISE GRD region/contact element ownership mismatch")
    if grd["materials"] != expected_materials:
        failures.append("DF-ISE GRD region material list mismatch")


def _cross_format_report(work: Path, bundle: dict[str, Any], msh: dict[str, Any]) -> dict[str, Any]:
    vtu = read_vtu(work / "device.vtu")
    vtk = read_legacy_vtk(work / "device.vtk")
    dat = parse_dat(work / "device.dat")
    grd = parse_grd(work / "device.grd")
    dat_fields = {name: entry["values"] for name, entry in dat.items()}
    fields = compare_nodal_fields(bundle, {
        "msh": _msh_fields(bundle, msh),
        "vtu": vtu["point_fields"],
        "vtk": vtk["point_fields"],
        "dat": dat_fields,
    })
    failures = list(fields["failures"])
    failures.extend(_csv_consistency(bundle, read_doping_csv(work / "doping.csv")))
    expected_physical = sorted(
        item["physical_name"] for key in ("regions", "materials", "contacts", "interfaces")
        for item in bundle[key])
    if msh["physical_names"] != expected_physical:
        failures.append("MSH Physical Names mismatch")
    _compare_msh_topology(bundle, msh, failures)
    expected_vtk = _expected_vtk(bundle)
    _compare_vtk_payload("VTU", vtu, expected_vtk, failures)
    _compare_vtk_payload("VTK", vtk, expected_vtk, failures)
    expected_cells = len(bundle["triangles"]) + len(bundle["lines"])
    if vtu["point_count"] != len(bundle["nodes"]) or vtu["cell_count"] != expected_cells:
        failures.append("VTU point/cell count mismatch")
    if grd["vertex_count"] != len(bundle["nodes"]):
        failures.append("DF-ISE GRD vertex count mismatch")
    if set(grd["regions"]) != {
            *(item["name"] for item in bundle["regions"]),
            *(item["name"] for item in bundle["contacts"]),
            *(item["name"] for item in bundle["interfaces"])}:
        failures.append("DF-ISE GRD region/contact names mismatch")
    _compare_dfise(bundle, grd, failures)
    return {
        "schema": "vela.sentaurus_gmsh_cross_format.v1",
        "status": "passed" if not failures else "failed",
        "failures": failures,
        "nodal_fields": fields,
        "msh_physical_names": msh["physical_names"],
        "vtu": {"point_count": vtu["point_count"], "cell_count": vtu["cell_count"],
                "cell_type_counts": {str(key): value for key, value in Counter(vtu["cell_types"]).items()}},
        "vtk": {"point_count": vtk["point_count"], "cell_count": vtk["cell_count"],
                "cell_type_counts": {str(key): value for key, value in Counter(vtk["cell_types"]).items()}},
        "dfise": {"vertex_count": grd["vertex_count"], "regions": grd["regions"],
                  "datasets": sorted(dat)},
    }


def run_pipeline(
    args: argparse.Namespace,
    *,
    device_ir_override: dict[str, Any] | None = None,
    frontend_manifest: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Run the Gmsh/export pipeline from a source deck or a frozen frontend IR.

    ``sentaurus_import.py device`` parses the SDE frontend exactly once and
    passes the resulting Device IR here.  This keeps evaluated-log semantics
    intact: the Gmsh backend must never silently reparse the original source
    with the offline frontend.
    """
    output = args.output_dir.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    for stale in ("run_manifest.json", "fail_report.json"):
        path = output / stale
        if path.exists():
            path.unlink()
    gmsh = None
    initialized = False
    try:
        try:
            gmsh, runtime = _load_runtime(args)
        except GmshRuntimeError as error:
            raise PipelineError(
                str(error), "gmsh_runtime", error.details) from error
        gmsh.initialize([f"vela-sentaurus-gmsh-{args.capability_level.lower()}", "-nopopup"], False)
        initialized = True
        library_version = gmsh.option.getString("General.Version")
        if library_version != EXPECTED_GMSH_VERSION:
            raise PipelineError(
                "loaded Gmsh shared library version does not match the pin",
                "gmsh_runtime", {"library_version": library_version})
        runtime["versions"]["library"] = library_version
        if device_ir_override is None:
            try:
                v1 = parse_sde_device_ir(args.sde.resolve())
                device_ir = migrate_device_ir_v1(v1)
            except (SdeParseError, ValueError) as error:
                raise PipelineError(str(error), "sde_parse") from error
        else:
            frozen = json.loads(json.dumps(device_ir_override, allow_nan=False))
            device_ir = (
                frozen if frozen.get("schema") == "vela.device_ir.v2"
                else migrate_device_ir_v1(frozen)
            )

        with tempfile.TemporaryDirectory(
                prefix=f".vela_gmsh_{args.capability_level.lower()}_",
                dir=output.parent) as tmp:
            work = Path(tmp)
            try:
                builder = {
                    "M1": build_m1_bundle,
                    "M2": build_m2_bundle,
                    "M3": build_m3_bundle,
                    "M5": build_m5_bundle,
                }[args.capability_level]
                if args.capability_level == "M5":
                    bundle, qualification = builder(
                        gmsh, device_ir, args.mesh_size,
                        args.mesh_strategy, args.transfinite_diagonal)
                else:
                    bundle, qualification = builder(gmsh, device_ir, args.mesh_size)
                write_json(work / "device_ir.json", device_ir)
                write_json(work / "mesh_field_bundle.json", bundle)
                write_json(work / "id_mapping.json", bundle["id_mapping"])
                write_json(work / "mesh.json", bundle_to_vela_mesh(bundle))
                write_doping_csv(work / "doping.csv", bundle)
                write_vtu(work / "device.vtu", bundle)
                write_legacy_vtk(work / "device.vtk", bundle)
                write_dfise(work / "device.grd", work / "device.dat", bundle)
                write_msh_with_node_data(gmsh, work / "device.msh", bundle)
                msh_readback = read_msh_roundtrip(gmsh, work / "device.msh")
                cross_format = _cross_format_report(work, bundle, msh_readback)
                if cross_format["status"] != "passed":
                    raise PipelineError(
                        "cross-format round-trip validation failed", "roundtrip", cross_format)
                write_json(work / "mesh_qualification.json", qualification)
                if args.capability_level == "M3":
                    write_json(
                        work / "refinement_field_dag.json",
                        bundle["provenance"]["refinement_field_dag"],
                    )
                if bundle["provenance"].get("transfinite_contract") is not None:
                    write_json(
                        work / "transfinite_mesh_contract.json",
                        bundle["provenance"]["transfinite_contract"],
                    )
                write_json(work / "cross_format_report.json", cross_format)
                write_json(work / "gmsh_runtime_manifest.json", runtime)
            except GmshBackendError as error:
                raise PipelineError(str(error), "gmsh_backend", error.details) from error
            except (OSError, ValueError) as error:
                raise PipelineError(str(error), "export") from error

            generated = sorted(path.name for path in work.iterdir() if path.is_file())
            hashes = {name: _file_sha256(work / name) for name in generated}
            manifest = {
                "schema": RUN_SCHEMA, "status": "generated",
                "milestone": args.capability_level,
                "sde": str(args.sde.resolve()), "sde_sha256": _file_sha256(args.sde.resolve()),
                "mesh_backend": "gmsh", "gmsh_version": EXPECTED_GMSH_VERSION,
                "mesh_strategy": bundle["provenance"].get("mesh_strategy"),
                "frontend": frontend_manifest,
                "generated": generated, "sha256": hashes,
                "mesh_qualification": qualification,
                "cross_format": cross_format,
                "capability": ({
                    "supported": ["strictly-separated Rectangle", "pick-point contact",
                                  "constant NoReplace region/material doping"],
                    "fail_closed": ["overlap/ABA/BAB", "touching interface", "insert-vertex",
                                    "Gaussian", "replacement", "refinement"],
                } if args.capability_level == "M1" else {
                    "supported": ["Rectangle ABA/BAB ownership", "conformal interface",
                                  "insert-vertex pick-point contact",
                                  "ordered Replace/LocalReplace/NoReplace",
                                  "line-baseline Gaussian"] + (
                                      ["Line/Rectangle/Polygon refinement window",
                                       "region/material/window anisotropic placement",
                                       "MaxLenInt DoubleSide/UseRegionNames",
                                       "auditable Gmsh field DAG"]
                                      if args.capability_level in {"M3", "M5"} else []) + (
                                      ["auto Transfinite Rectangle selection",
                                       "first-order right Tri3 structured mesh",
                                       "Left/Right/AlternateLeft/AlternateRight diagonals"]
                                      if args.capability_level == "M5" else []),
                    "fail_closed": ["AB/ABiA/ABiB/XX", "Polygon body",
                                    "contact entity/body/remove/imprint"] + (
                                        ["MaxTransDiff/MaxGradient phase 4B", "MaxInterval",
                                         "multibox", "offsetting", "MaxLenInt contact interface"]
                                        if args.capability_level in {"M3", "M5"}
                                        else ["refinement"]) + ([
                                        "Transfinite multi-body/insert-vertex/function mapping"]
                                        if args.capability_level == "M5" else []) + [
                                    "constant-profile decay",
                                    "analytical restriction decay"],
                }),
            }
            write_json(work / "run_manifest.json", manifest)
            for name in generated:
                os.replace(work / name, output / name)
            os.replace(work / "run_manifest.json", output / "run_manifest.json")
            return manifest
    finally:
        if initialized and gmsh is not None:
            gmsh.finalize()


def write_fail_report(output: Path, error: BaseException) -> Path:
    stage = error.stage if isinstance(error, PipelineError) else "unexpected"
    details = error.details if isinstance(error, PipelineError) else {}
    report = {
        "schema": FAIL_SCHEMA, "status": "failed", "stage": stage,
        "message": str(error), "details": details,
        "recommendation": (
            "install/configure Gmsh 4.15.2 or restrict the SDE deck to the selected "
            "capability subset; unsupported semantics require the named later milestone"),
    }
    path = output.resolve() / "fail_report.json"
    _strict_json(path, report)
    manifest = output.resolve() / "run_manifest.json"
    if manifest.exists():
        manifest.unlink()
    return path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sde", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--runtime-dir", type=Path,
                        default=REPO / "build" / "gmsh-python-runtime")
    parser.add_argument("--gmsh-python", type=Path)
    parser.add_argument("--gmsh-library", type=Path)
    parser.add_argument("--gmsh-executable", type=Path)
    parser.add_argument("--mesh-size", type=float, default=0.2)
    parser.add_argument("--capability-level", choices=("M1", "M2", "M3", "M5"), default="M5",
                        help="freeze the accepted semantic subset (default: M5)")
    parser.add_argument(
        "--mesh-strategy", choices=("auto", "transfinite-rect", "unstructured"),
        default="auto", help="M5 mesh selection policy (default: auto)")
    parser.add_argument(
        "--transfinite-diagonal",
        choices=("Left", "Right", "AlternateLeft", "AlternateRight"),
        default="AlternateRight", help="Gmsh Transfinite surface arrangement")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        manifest = run_pipeline(args)
    except Exception as error:  # noqa: BLE001 - CLI always leaves a strict report
        try:
            report = write_fail_report(args.output_dir, error)
        except Exception as report_error:  # noqa: BLE001
            print(f"pipeline failed and fail report could not be written: {report_error}",
                  file=sys.stderr)
            return 3
        print(f"Gmsh preprocessing failed; inspect {report}", file=sys.stderr)
        return 2
    print(json.dumps(manifest, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
