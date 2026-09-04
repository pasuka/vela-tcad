#!/usr/bin/env python3
"""Run a Sentaurus PN2D vector-current oracle on the exact Vela mesh."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Sequence


WORKTREE = Path(__file__).resolve().parents[1]
DEFAULT_DATA_REPO = Path(r"D:\code-repo\vela-tcad")
DEFAULT_OUTPUT = (
    WORKTREE
    / "build-release/reference_tcad/pn2d_sentaurus2022/sentaurus_vm_runs"
)
DEFAULT_IMPORTER = WORKTREE / "build-release/sentaurus_import.exe"
DEFAULT_REMOTE_ROOT = "~/sentaurus_runs/vela_oracle_2022"
RUN_SCHEMA = "vela.pn2d_same_mesh_sdevice_vector.v1"
COORDINATE_TOLERANCE_UM = 1.0e-12
RETURNED_ARTIFACTS = (
    "pn2d_same_mesh.grd",
    "pn2d_same_mesh.dat",
    "pn2d_same_mesh.tdr",
    "pn2d_same_mesh_m20_des.tdr",
    "pn2d_same_mesh_m20_des.log",
    "pn2d_same_mesh_m20.plt",
    "run_tdx.out",
    "run_sdevice.out",
)
_SAFE_RUN_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
_SAFE_REMOTE = re.compile(r"^[A-Za-z0-9_./~:-]+$")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8", newline="\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_edges(triangles: Iterable[tuple[int, int, int]]) -> list[tuple[int, int]]:
    return sorted(
        {
            tuple(sorted(edge))
            for triangle in triangles
            for edge in (
                (triangle[0], triangle[1]),
                (triangle[1], triangle[2]),
                (triangle[2], triangle[0]),
            )
        }
    )


def edge_reference(edge: tuple[int, int], edge_ids: dict[tuple[int, int], int]) -> int:
    canonical = tuple(sorted(edge))
    edge_id = edge_ids[canonical]
    return edge_id if edge == canonical else -(edge_id + 1)


def load_mesh(mesh_path: Path) -> tuple[dict[int, tuple[float, float]], list[tuple[int, int, int]], dict[str, list[int]]]:
    mesh = read_json(mesh_path)
    nodes = {int(row["id"]): (float(row["x"]), float(row["y"])) for row in mesh["nodes"]}
    if sorted(nodes) != list(range(len(nodes))):
        raise ValueError("Vela mesh node IDs must be contiguous from zero")
    triangles = [tuple(int(value) for value in row["node_ids"]) for row in mesh["triangles"]]
    if any(len(triangle) != 3 or any(node not in nodes for node in triangle) for triangle in triangles):
        raise ValueError("Vela mesh has an invalid triangle")
    contacts = {str(row["name"]).lower(): [int(value) for value in row["node_ids"]] for row in mesh["contacts"]}
    if set(contacts) != {"anode", "cathode"}:
        raise ValueError("PN mesh must contain only anode and cathode contacts")
    return nodes, triangles, contacts


def contact_segments(
    nodes: dict[int, tuple[float, float]],
    contact_nodes: list[int],
    edge_set: set[tuple[int, int]],
) -> list[tuple[int, int]]:
    ordered = sorted(contact_nodes, key=lambda node: (nodes[node][1], nodes[node][0]))
    segments = [tuple(sorted((left, right))) for left, right in zip(ordered, ordered[1:])]
    if not segments or any(segment not in edge_set for segment in segments):
        raise ValueError("contact nodes do not form consecutive boundary edges")
    return segments


def render_dfise_grid(
    nodes: dict[int, tuple[float, float]],
    triangles: list[tuple[int, int, int]],
    contacts: dict[str, list[int]],
) -> str:
    edges = canonical_edges(triangles)
    edge_ids = {edge: index for index, edge in enumerate(edges)}
    incidence = Counter(
        tuple(sorted(edge))
        for triangle in triangles
        for edge in (
            (triangle[0], triangle[1]),
            (triangle[1], triangle[2]),
            (triangle[2], triangle[0]),
        )
    )
    locations = "".join("e" if incidence[edge] == 1 else "i" for edge in edges)
    edge_set = set(edges)
    contact_order = ("cathode", "anode")
    contact_edges = {
        name: contact_segments(nodes, contacts[name], edge_set) for name in contact_order
    }
    vertices = "\n".join(
        f" {nodes[node][0]:.15e} {nodes[node][1]:.15e}" for node in sorted(nodes)
    )
    edge_rows = "\n".join(f" {start} {end}" for start, end in edges)
    triangle_rows = [
        " 2 "
        + " ".join(
            str(edge_reference(edge, edge_ids))
            for edge in (
                (triangle[0], triangle[1]),
                (triangle[1], triangle[2]),
                (triangle[2], triangle[0]),
            )
        )
        for triangle in triangles
    ]
    contact_rows: list[str] = []
    contact_element_ids: dict[str, list[int]] = {}
    next_element = len(triangles)
    for name in contact_order:
        ids = []
        for start, end in contact_edges[name]:
            contact_rows.append(f" 1 {start} {end}")
            ids.append(next_element)
            next_element += 1
        contact_element_ids[name] = ids
    element_rows = "\n".join(triangle_rows + contact_rows)

    def id_block(values: list[int]) -> str:
        return " ".join(str(value) for value in values)

    return f"""DF-ISE text

Info {{
  version = 1.1
  type = grid
  dimension = 2
  nb_vertices = {len(nodes)}
  nb_edges = {len(edges)}
  nb_faces = 0
  nb_elements = {next_element}
  nb_regions = 3
  regions = [ \"R.Si\" \"Cathode\" \"Anode\" ]
  materials = [ Silicon Contact Contact ]
}}

Data {{
  CoordSystem {{
    translate = [ 0 0 0 ]
    transform = [ 1 0 0 0 1 0 0 0 1 ]
  }}

  Vertices ({len(nodes)}) {{
{vertices}
  }}

  Edges ({len(edges)}) {{
{edge_rows}
  }}

  Locations ({len(edges)}) {{
{locations}
  }}

  Elements ({next_element}) {{
{element_rows}
  }}

  Region (\"R.Si\") {{
    material = Silicon
    Elements ({len(triangles)}) {{
 {id_block(list(range(len(triangles))))}
    }}
  }}

  Region (\"Cathode\") {{
    material = Contact
    Elements ({len(contact_element_ids['cathode'])}) {{
 {id_block(contact_element_ids['cathode'])}
    }}
  }}

  Region (\"Anode\") {{
    material = Contact
    Elements ({len(contact_element_ids['anode'])}) {{
 {id_block(contact_element_ids['anode'])}
    }}
  }}
}}
"""


def render_dfise_doping(
    nodes: dict[int, tuple[float, float]],
    triangles: list[tuple[int, int, int]],
    contacts: dict[str, list[int]],
    doping_path: Path,
) -> str:
    doping_rows = {int(row["node_id"]): row for row in read_csv(doping_path)}
    if set(doping_rows) != set(nodes):
        raise ValueError("doping node IDs do not exactly match the Vela mesh")
    datasets = {
        "DopingConcentration": [
            float(doping_rows[node]["donors_cm3"]) - float(doping_rows[node]["acceptors_cm3"])
            for node in sorted(nodes)
        ],
        "PhosphorusActiveConcentration": [
            float(doping_rows[node]["donors_cm3"]) for node in sorted(nodes)
        ],
        "BoronActiveConcentration": [
            float(doping_rows[node]["acceptors_cm3"]) for node in sorted(nodes)
        ],
    }
    blocks = []
    for name, values in datasets.items():
        value_rows = "\n".join(f" {value:.15e}" for value in values)
        blocks.append(
            f"""  Dataset (\"{name}\") {{
    function = {name}
    type = scalar
    dimension = 1
    location = vertex
    validity = [ \"R.Si\" ]
    Values ({len(nodes)}) {{
{value_rows}
    }}
  }}"""
        )
    names = " ".join(f'"{name}"' for name in datasets)
    functions = " ".join(datasets)
    edge_set = set(canonical_edges(triangles))
    contact_element_count = sum(
        len(contact_segments(nodes, contacts[name], edge_set))
        for name in ("cathode", "anode")
    )
    return f"""DF-ISE text

Info {{
  version = 1.0
  type = dataset
  dimension = 2
  nb_vertices = {len(nodes)}
  nb_edges = {len(edge_set)}
  nb_faces = 0
  nb_elements = {len(triangles) + contact_element_count}
  nb_regions = 3
  datasets = [ {names} ]
  functions = [ {functions} ]
}}

Data {{
{chr(10).join(blocks)}
}}
"""


def render_sdevice_deck(target_bias: float) -> str:
    return f"""File {{
  Grid = \"pn2d_same_mesh.tdr\"
  Doping = \"pn2d_same_mesh.tdr\"
  Parameter = \"models.par\"
  Plot = \"pn2d_same_mesh_m20_des.tdr\"
  Current = \"pn2d_same_mesh_m20.plt\"
  Output = \"pn2d_same_mesh_m20_des.log\"
}}

Electrode {{
  {{ Name=\"Anode\" Voltage=0.0 }}
  {{ Name=\"Cathode\" Voltage=0.0 }}
}}

Physics {{
  Mobility(DopingDependence HighFieldSaturation)
  Recombination(SRH Avalanche(VanOverstraeten))
  EffectiveIntrinsicDensity(OldSlotboom)
}}

Plot {{
  Potential eQuasiFermi hQuasiFermi eDensity hDensity
  ElectricField/Vector
  eCurrentDensity/Vector hCurrentDensity/Vector TotalCurrentDensity/Vector
  Doping DonorConcentration AcceptorConcentration SpaceCharge
  SRHRecombination eAlphaAvalanche hAlphaAvalanche AvalancheGeneration
  eMobility hMobility eVelocity hVelocity
}}

Math {{
  Extrapolate RelErrControl Digits=5 Iterations=80 NotDamped=100 Method=Blocked
}}

Solve {{
  Coupled(Iterations=100) {{ Poisson }}
  Coupled(Iterations=100) {{ Poisson Electron Hole }}
  Quasistationary(
    InitialStep=1e-4 MinStep=1e-10 MaxStep=0.02 Increment=1.2 Decrement=2.0
    Goal {{ Name=\"Anode\" Voltage={target_bias:.15g} }}
  ) {{ Coupled {{ Poisson Electron Hole }} }}
}}
"""


def build_bundle(mesh_path: Path, doping_path: Path, models_path: Path, bundle: Path, target_bias: float) -> dict[str, Any]:
    nodes, triangles, contacts = load_mesh(mesh_path)
    bundle.mkdir(parents=True, exist_ok=True)
    paths = {
        "grid": bundle / "pn2d_same_mesh.grd",
        "doping": bundle / "pn2d_same_mesh.dat",
        "deck": bundle / "pn2d_same_mesh_sdevice.cmd",
        "models": bundle / "models.par",
    }
    paths["grid"].write_text(render_dfise_grid(nodes, triangles, contacts), encoding="utf-8")
    paths["doping"].write_text(
        render_dfise_doping(nodes, triangles, contacts, doping_path),
        encoding="utf-8",
    )
    paths["deck"].write_text(render_sdevice_deck(target_bias), encoding="utf-8")
    shutil.copyfile(models_path, paths["models"])
    return {
        "nodes": len(nodes),
        "triangles": len(triangles),
        "contacts": contacts,
        "files": {name: str(path) for name, path in paths.items()},
        "sha256": {name: sha256(path) for name, path in paths.items()},
    }


def default_openssh(name: str) -> str:
    candidate = Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32/OpenSSH" / f"{name}.exe"
    return str(candidate) if candidate.is_file() else (shutil.which(name) or name)


def run_checked(argv: Sequence[str]) -> None:
    subprocess.run(list(argv), check=True)


def coordinate_mapping(
    canonical_nodes: dict[int, tuple[float, float]], export_dir: Path
) -> dict[int, int]:
    exported = read_csv(export_dir / "nodes.csv")
    if len(exported) != len(canonical_nodes):
        raise ValueError(
            f"same-mesh node-count mismatch: expected {len(canonical_nodes)}, got {len(exported)}"
        )
    mapping: dict[int, int] = {}
    used: set[int] = set()
    for row in exported:
        point = (float(row["x_um"]), float(row["y_um"]))
        matches = [
            node
            for node, expected in canonical_nodes.items()
            if abs(point[0] - expected[0]) <= COORDINATE_TOLERANCE_UM
            and abs(point[1] - expected[1]) <= COORDINATE_TOLERANCE_UM
        ]
        if len(matches) != 1 or matches[0] in used:
            raise ValueError(f"exported node {row['id']} has no unique canonical coordinate")
        mapping[int(row["id"])] = matches[0]
        used.add(matches[0])
    return mapping


def align_export(export_dir: Path, aligned_dir: Path, mapping: dict[int, int]) -> None:
    if aligned_dir.exists():
        shutil.rmtree(aligned_dir)
    shutil.copytree(export_dir, aligned_dir)
    nodes = read_csv(aligned_dir / "nodes.csv")
    aligned_nodes = [
        {**row, "id": mapping[int(row["id"])]} for row in nodes
    ]
    write_csv(aligned_dir / "nodes.csv", sorted(aligned_nodes, key=lambda row: int(row["id"])), list(nodes[0]))
    for name, id_columns in (
        ("elements.csv", ("node0", "node1", "node2")),
        ("doping.csv", ("node_id",)),
    ):
        path = aligned_dir / name
        rows = read_csv(path)
        for row in rows:
            for column in id_columns:
                row[column] = str(mapping[int(row[column])])
        if "node_id" in rows[0]:
            rows.sort(key=lambda row: int(row["node_id"]))
        write_csv(path, rows, list(rows[0]))
    contacts_path = aligned_dir / "contacts.csv"
    contacts = read_csv(contacts_path)
    for row in contacts:
        row["node_ids"] = ";".join(str(mapping[int(value)]) for value in row["node_ids"].split(";") if value)
    write_csv(contacts_path, contacts, list(contacts[0]))
    for path in (aligned_dir / "fields").glob("*.csv"):
        rows = read_csv(path)
        if not rows or "node_id" not in rows[0]:
            continue
        for row in rows:
            row["node_id"] = str(mapping[int(row["node_id"])] )
        rows.sort(key=lambda row: int(row["node_id"]))
        write_csv(path, rows, list(rows[0]))


def validate_topology(
    canonical_triangles: list[tuple[int, int, int]], aligned_dir: Path
) -> dict[str, Any]:
    rows = read_csv(aligned_dir / "elements.csv")
    actual = {
        tuple(sorted(int(row[name]) for name in ("node0", "node1", "node2")))
        for row in rows
        if row["region"] == "R.Si"
    }
    expected = {tuple(sorted(triangle)) for triangle in canonical_triangles}
    if actual != expected:
        raise ValueError("same-mesh triangle connectivity mismatch after TDX round trip")
    vector_fields = {}
    for field in ("eCurrentDensity", "hCurrentDensity", "TotalCurrentDensity"):
        path = aligned_dir / "fields" / f"{field}_region0.csv"
        rows = read_csv(path)
        if not rows or not {"node_id", "component0", "component1"} <= set(rows[0]):
            raise ValueError(f"{field} is not a two-component node vector")
        if len(rows) != len({node for triangle in expected for node in triangle}):
            raise ValueError(f"{field} does not cover every mesh node")
        vector_fields[field] = {"path": str(path), "sha256": sha256(path), "nodes": len(rows)}
    return {
        "passed": True,
        "node_count": len({node for triangle in expected for node in triangle}),
        "triangle_count": len(actual),
        "vector_fields": vector_fields,
    }


def write_classical_dd_state(aligned_dir: Path, state_path: Path) -> None:
    """Write the six-column Vela state used by classical DD diagnostics."""
    field_names = {
        "psi": "ElectrostaticPotential",
        "phin": "eQuasiFermiPotential",
        "phip": "hQuasiFermiPotential",
        "electrons_m3": "eDensity",
        "holes_m3": "hDensity",
    }
    columns: dict[str, dict[int, float]] = {}
    for output_name, field_name in field_names.items():
        rows = read_csv(aligned_dir / "fields" / f"{field_name}_region0.csv")
        columns[output_name] = {
            int(row["node_id"]): float(row["component0"]) for row in rows
        }
    node_ids = set(columns["psi"])
    if any(set(values) != node_ids for values in columns.values()):
        raise ValueError("classical DD restart fields do not share identical node support")
    rows = []
    for node in sorted(node_ids):
        rows.append(
            {
                "node_id": node,
                "psi": columns["psi"][node],
                "phin": columns["phin"][node],
                "phip": columns["phip"][node],
                "electrons_m3": columns["electrons_m3"][node] * 1.0e6,
                "holes_m3": columns["holes_m3"][node] * 1.0e6,
            }
        )
    write_csv(
        state_path,
        rows,
        ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"],
    )


def finalize_local_artifacts(manifest: dict[str, Any]) -> None:
    artifacts = Path(manifest["artifacts_dir"])
    export_dir = artifacts / "neutral_export"
    nodes, triangles, _ = load_mesh(Path(manifest["mesh_file"]))
    mapping = coordinate_mapping(nodes, export_dir)
    aligned = artifacts / "aligned_export"
    align_export(export_dir, aligned, mapping)
    topology = validate_topology(triangles, aligned)
    state = artifacts / "sdevice_state.csv"
    write_classical_dd_state(aligned, state)
    missing = [name for name in RETURNED_ARTIFACTS if not (artifacts / name).is_file()]
    if missing:
        raise ValueError("missing returned artifacts: " + ", ".join(missing))
    manifest.update(
        {
            "neutral_export": str(export_dir),
            "aligned_export": str(aligned),
            "sdevice_state": str(state),
            "coordinate_mapping": {str(key): value for key, value in mapping.items()},
            "topology_gate": topology,
            "returned_sha256": {
                name: sha256(artifacts / name) for name in RETURNED_ARTIFACTS
            },
            "aligned_nodes_sha256": sha256(aligned / "nodes.csv"),
            "aligned_elements_sha256": sha256(aligned / "elements.csv"),
            "sdevice_state_sha256": sha256(state),
            "passed": True,
        }
    )


def run_live(
    manifest: dict[str, Any], ssh_target: str, ssh_bin: str, scp_bin: str, importer: Path
) -> None:
    remote_dir = manifest["remote_dir"]
    bundle = Path(manifest["bundle_dir"])
    artifacts = Path(manifest["artifacts_dir"])
    artifacts.mkdir(parents=True, exist_ok=True)
    run_checked([ssh_bin, ssh_target, f"mkdir -p {remote_dir}"])
    for name in ("pn2d_same_mesh.grd", "pn2d_same_mesh.dat", "pn2d_same_mesh_sdevice.cmd", "models.par"):
        run_checked([scp_bin, str(bundle / name), f"{ssh_target}:{remote_dir}/"])
    run_checked([
        ssh_bin,
        ssh_target,
        f"cd {remote_dir} && test -s pn2d_same_mesh.tdr || "
        "tdx -d pn2d_same_mesh.grd pn2d_same_mesh.dat pn2d_same_mesh.tdr > run_tdx.out 2>&1",
    ])
    run_checked([
        ssh_bin,
        ssh_target,
        f"cd {remote_dir} && sdevice pn2d_same_mesh_sdevice.cmd > run_sdevice.out 2>&1",
    ])
    for name in RETURNED_ARTIFACTS:
        run_checked([scp_bin, f"{ssh_target}:{remote_dir}/{name}", str(artifacts) + os.sep])
    export_dir = artifacts / "neutral_export"
    run_checked([
        str(importer.resolve()), "--tdr", str(artifacts / "pn2d_same_mesh_m20_des.tdr"),
        "--export-dir", str(export_dir), "--compensated-doping-policy", "reported",
    ])
    finalize_local_artifacts(manifest)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-repo", type=Path, default=DEFAULT_DATA_REPO)
    parser.add_argument("--mesh-file", type=Path, default=None)
    parser.add_argument("--doping-file", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--target-bias", type=float, default=-20.0)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=None)
    parser.add_argument("--scp-bin", default=None)
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    parser.add_argument("--importer", type=Path, default=DEFAULT_IMPORTER)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--reuse-local-export", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    run_id = args.run_id or datetime.now().strftime("pn2d_same_mesh_vector_%Y%m%d_%H%M%S")
    if not _SAFE_RUN_ID.fullmatch(run_id):
        raise ValueError("run ID contains unsupported characters")
    remote_root = args.remote_root.rstrip("/")
    if not _SAFE_REMOTE.fullmatch(remote_root):
        raise ValueError("remote root contains unsupported characters")
    source = args.data_repo / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/imported_reference/vela"
    mesh = args.mesh_file or (source / "mesh.json")
    doping = args.doping_file or (source / "doping.csv")
    models = WORKTREE / "reference_tcad/pn2d_sentaurus2018_coarse7x3/source/models.par"
    run_root = (args.output_dir / run_id).resolve()
    bundle = run_root / "source"
    artifacts = run_root / "artifacts"
    build = build_bundle(mesh, doping, models, bundle, args.target_bias)
    manifest: dict[str, Any] = {
        "schema": RUN_SCHEMA,
        "run_id": run_id,
        "dry_run": args.dry_run,
        "sentaurus_version": "T-2022.03-SP2",
        "sdevice_input_format": "tdr_converted_from_explicit_dfise",
        "tdx_resume_policy": "reuse_nonempty_tdr_for_same_run_id",
        "reused_local_export": args.reuse_local_export,
        "target_bias_V": args.target_bias,
        "mesh_file": str(mesh.resolve()),
        "mesh_sha256": sha256(mesh),
        "doping_file": str(doping.resolve()),
        "doping_sha256": sha256(doping),
        "bundle_dir": str(bundle),
        "artifacts_dir": str(artifacts),
        "remote_dir": f"{remote_root}/{run_id}",
        "build": build,
        "passed": None,
    }
    manifest_path = run_root / "manifest.json"
    try:
        if not args.dry_run:
            if args.reuse_local_export:
                finalize_local_artifacts(manifest)
            else:
                run_live(
                    manifest,
                    args.ssh_target,
                    args.ssh_bin or default_openssh("ssh"),
                    args.scp_bin or default_openssh("scp"),
                    args.importer,
                )
        write_json(manifest_path, manifest)
        print(json.dumps({"manifest": str(manifest_path), "passed": manifest["passed"]}))
        return 0 if args.dry_run or manifest["passed"] else 1
    except Exception as error:
        manifest["passed"] = False
        manifest["error"] = str(error)
        write_json(manifest_path, manifest)
        print(str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
