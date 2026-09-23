#!/usr/bin/env python3
"""Run the local-refinement sensitivity workflow for the Genius NPN BJT."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
import state_archive


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
DEFAULT_ROOT = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "genius_bjt_sentaurus2022"
    / "mesh_sensitivity"
    / "local_refined"
)
DEFAULT_RUNNER = REPO / "build-release" / "vela_example_runner.exe"
DEFAULT_IMPORTER = REPO / "build-release" / "sentaurus_import.exe"
DEFAULT_PYTHON = Path(sys.executable)
TARGET_BIASES = [index / 10.0 for index in range(31)]


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(command: list[str], *, log: Path | None = None) -> str:
    process = subprocess.run(
        command,
        cwd=REPO,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    if log is not None:
        log.parent.mkdir(parents=True, exist_ok=True)
        log.with_suffix(".stdout.log").write_text(process.stdout, encoding="utf-8")
        log.with_suffix(".stderr.log").write_text(process.stderr, encoding="utf-8")
    if process.returncode != 0:
        detail = process.stderr.strip() or process.stdout.strip()
        raise RuntimeError(f"command failed ({process.returncode}): {' '.join(command)}\n{detail}")
    return process.stdout


def find_field(export: Path, name: str) -> Path:
    matches = sorted((export / "fields").glob(f"{name}_region*.csv"))
    matches = [path for path in matches if not path.stem.endswith("_cells")]
    if len(matches) != 1:
        raise ValueError(f"expected one node field for {name} in {export}, found {matches}")
    return matches[0]


def field_map(path: Path) -> dict[int, float]:
    result: dict[int, float] = {}
    for row in read_csv(path):
        node = int(row.get("id", row.get("node_id", "")))
        result[node] = float(row["component0"])
    return result


def make_seed_state(export: Path, output: Path, mesh_file: Path) -> int:
    names = {
        "psi": "ElectrostaticPotential",
        "phin": "eQuasiFermiPotential",
        "phip": "hQuasiFermiPotential",
        "electrons_cm3": "eDensity",
        "holes_cm3": "hDensity",
    }
    fields = {key: field_map(find_field(export, value)) for key, value in names.items()}
    node_ids = sorted(fields["psi"])
    if any(sorted(field) != node_ids for field in fields.values()):
        raise ValueError("SDevice seed fields do not share one node-id set")
    rows = []
    for node in node_ids:
        rows.append(
            {
                "node_id": node,
                "psi": fields["psi"][node],
                "phin": fields["phin"][node],
                "phip": fields["phip"][node],
                "electrons_m3": fields["electrons_cm3"][node] * 1.0e6,
                "holes_m3": fields["holes_cm3"][node] * 1.0e6,
            }
        )
    mesh = read_json(mesh_file)
    if len(rows) != len(mesh["nodes"]):
        raise ValueError("SDevice seed node count differs from the converted mesh")
    state_archive.write(output, state_archive.rows_to_fields(rows), dict(
        mode="dd", mesh_sha256=state_archive.mesh_identity(mesh, 1e-6),
        potential_origin_V=0., source_fields_sha256={
            name: sha256(find_field(export, name)) for name in names.values()}))
    return len(rows)


def prepare(args: argparse.Namespace) -> dict:
    from audit_genius_bjt_sentaurus_reference import audit_structure

    root = args.artifact_root.resolve()
    raw = root / "raw"
    required = {
        "mesh": raw / "bjt_msh.tdr",
        "vbe070": raw / "bjt_m1_vbe070_des.tdr",
        "vce3": raw / "bjt_m1_vce_0003_des.tdr",
        "plt": raw / "bjt_m1_collector_bjt_m1_des.plt",
    }
    missing = [str(path) for path in required.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing returned Sentaurus artifacts: " + ", ".join(missing))

    exports = root / "sentaurus_exports"
    for name in ("mesh", "vbe070", "vce3"):
        run(
            [str(args.importer.resolve()), "--tdr", str(required[name]), "--export-dir", str(exports / name)],
            log=root / "logs" / f"import_{name}",
        )
    input_root = root / "vela" / "input"
    run(
        [
            str(args.python.resolve()),
            str(REPO / "scripts" / "convert_tcad_export.py"),
            "--input-dir",
            str(exports / "mesh"),
            "--output-dir",
            str(input_root),
            "--device",
            "bjt2d",
            "--simulation-types",
            "iv",
        ],
        log=root / "logs" / "convert_mesh",
    )
    count = make_seed_state(exports / "vbe070", root / "vela" / "sdevice_vbe070_seed.h5",
                            input_root / "mesh.json")
    metadata = read_json(exports / "mesh" / "metadata.json")
    structure_audit = audit_structure(exports / "mesh")
    if not structure_audit["pass"]:
        raise RuntimeError(f"refined structure audit failed: {structure_audit}")
    result = {
        "node_count": count,
        "mesh_metadata": metadata,
        "structure_audit": structure_audit,
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in required.items()},
        "mesh_sha256": sha256(input_root / "mesh.json"),
        "doping_sha256": sha256(input_root / "doping.csv"),
    }
    write_json(root / "prepare_manifest.json", result)
    return result


def patch_paths(cfg: dict, root: Path) -> None:
    cfg["mesh_file"] = str((root / "vela" / "input" / "mesh.json").resolve())
    cfg["node_doping_file"] = str((root / "vela" / "input" / "doping.csv").resolve())
    cfg["materials_file"] = str((FIXTURE / "vela" / "materials_sentaurus2022.json").resolve())


def set_contact(cfg: dict, name: str, bias: float) -> None:
    for contact in cfg["contacts"]:
        if contact["name"] == name:
            contact["bias"] = bias
            return
    raise KeyError(name)


def run_vela_config(args: argparse.Namespace, cfg: dict, name: str) -> dict:
    root = args.artifact_root.resolve()
    config = root / "vela" / "configs" / f"{name}.json"
    write_json(config, cfg)
    stdout = run(
        [str(args.runner.resolve()), "--config", str(config)],
        log=root / "logs" / name,
    )
    lines = [line for line in stdout.splitlines() if line.strip()]
    return json.loads(lines[-1]) if lines else {}


def spatial(args: argparse.Namespace) -> dict:
    root = args.artifact_root.resolve()
    final_state = root / "vela" / "m1_vce300_state.h5"
    if not final_state.is_file():
        raise FileNotFoundError(f"run the Vela collector sweep first: {final_state}")
    cfg = read_json(FIXTURE / "vela" / "configs" / "m1_spatial_vce3.json")
    patch_paths(cfg, root)
    cfg.pop("sweep", None)
    set_contact(cfg, "base", 0.70)
    set_contact(cfg, "collector", 3.0)
    cfg.update(
        {
            "simulation_type": "newton_solve_from_state",
            "state_file": str(final_state.resolve()),
            "output_state_file": str((root / "vela" / "m1_spatial_vce3_state.h5").resolve()),
            "output_vtk": str((root / "vela" / "m1_spatial_vce3.vtk").resolve()),
        }
    )
    status = run_vela_config(args, cfg, "m1_spatial_vce3_fixed_state")
    if not status.get("converged", False):
        raise RuntimeError(f"Vela 3 V fixed-state export failed: {status}")
    return status


def vela(args: argparse.Namespace) -> dict:
    root = args.artifact_root.resolve()
    seed = root / "vela" / "sdevice_vbe070_seed.h5"
    if not seed.is_file():
        raise FileNotFoundError(f"run prepare first: {seed}")

    base = read_json(FIXTURE / "vela" / "configs" / "m1_spatial_vce3.json")
    patch_paths(base, root)
    base.pop("sweep", None)
    set_contact(base, "base", 0.70)
    set_contact(base, "collector", 0.0)
    base.update(
        {
            "simulation_type": "newton_solve_from_state",
            "state_file": str(seed.resolve()),
            "output_state_file": str((root / "vela" / "m1_vbe070_state.h5").resolve()),
        }
    )
    seed_status = run_vela_config(args, base, "m1_seed_relaxation")
    if not seed_status.get("converged", False):
        raise RuntimeError(f"Vela seed relaxation failed: {seed_status}")

    sweep = read_json(FIXTURE / "vela" / "configs" / "m1_collector_sweep.json")
    patch_paths(sweep, root)
    sweep["output_csv"] = str((root / "vela" / "m1_collector.csv").resolve())
    sweep["sweep"]["initial_state_file"] = str((root / "vela" / "m1_vbe070_state.h5").resolve())
    sweep["sweep"]["write_state_file"] = str((root / "vela" / "m1_vce300_state.h5").resolve())
    sweep["sweep"]["diagnostics"]["terminal_balance"]["csv_file"] = str(
        (root / "vela" / "m1_collector_terminal_balance.csv").resolve()
    )
    sweep_status = run_vela_config(args, sweep, "m1_collector_sweep")

    spatial_status = spatial(args)

    result = {"seed": seed_status, "sweep": sweep_status, "spatial": spatial_status}
    write_json(root / "vela_run_manifest.json", result)
    return result


def select_terminal_rows(path: Path) -> list[dict[str, object]]:
    grouped: dict[int, list[dict[str, str]]] = {}
    for row in read_csv(path):
        grouped.setdefault(int(row["point_index"]), []).append(row)
    points = []
    for rows in grouped.values():
        by_contact = {row["contact"]: row for row in rows}
        if not {"collector", "base", "emitter"}.issubset(by_contact):
            continue
        points.append((float(rows[0]["bias_V"]), by_contact))
    selected = []
    for target in TARGET_BIASES:
        bias, by_contact = min(points, key=lambda item: abs(item[0] - target))
        if abs(bias - target) > 1.0e-8:
            raise ValueError(f"Vela terminal ledger has no exact VCE={target:.1f} V point")
        currents = {name: float(by_contact[name]["current_total_A_per_um"]) for name in by_contact}
        selected.append(
            {
                "VCE_V": target,
                "Ic_A_per_um": currents["collector"],
                "Ib_A_per_um": currents["base"],
                "Ie_A_per_um": currents["emitter"],
                "beta_abs": abs(currents["collector"] / currents["base"]),
                "kcl_abs_A_per_um": abs(sum(currents.values())),
            }
        )
    return selected


def maximum_log_error(lhs: list[dict], rhs: list[dict], key: str, start: float = 0.5) -> float:
    values = []
    for left, right in zip(lhs, rhs, strict=True):
        if float(left["VCE_V"]) < start:
            continue
        a, b = abs(float(left[key])), abs(float(right[key]))
        values.append(abs(math.log10(max(a, 1.0e-300) / max(b, 1.0e-300))))
    return max(values)


def analyze(args: argparse.Namespace) -> dict:
    from audit_genius_bjt_sentaurus_reference import normalize_plt

    root = args.artifact_root.resolve()
    refined_sentaurus, plt_audit = normalize_plt(
        root / "raw" / "bjt_m1_collector_bjt_m1_des.plt", "M1-refined"
    )
    refined_vela = select_terminal_rows(root / "vela" / "m1_collector_terminal_balance.csv")
    baseline_sentaurus = read_csv(FIXTURE / "reference_curves" / "bjt_m1_output.csv")
    write_csv(root / "comparison" / "sentaurus_refined_curve.csv", refined_sentaurus)
    write_csv(root / "comparison" / "vela_refined_curve.csv", refined_vela)

    vtk = root / "vela" / "m1_spatial_vce3.vtk"
    transport_dir = root / "comparison" / "transport"
    run(
        [
            str(args.python.resolve()),
            str(REPO / "scripts" / "compare_genius_bjt_transport_fields.py"),
            "--reference-root",
            str(FIXTURE),
            "--sentaurus-fields-root",
            str(root / "sentaurus_exports" / "vce3"),
            "--vela-vtk",
            str(vtk),
            "--output-dir",
            str(transport_dir),
        ],
        log=root / "logs" / "transport_compare",
    )
    transport = read_json(transport_dir / "transport_source_comparison.json")
    refined_spatial = read_json(root / "comparison" / "spatial" / "spatial_comparison_summary.json")
    baseline_transport = read_json(FIXTURE / "comparison" / "transport_source_comparison.json")
    baseline_spatial = read_json(FIXTURE / "comparison" / "spatial_comparison_summary.json")

    report = {
        "schema_version": 1,
        "bias_contract": "VBE=0.70 V; VCE=0.0..3.0 V in 0.1 V increments",
        "mesh": {
            "baseline_nodes": 5611,
            "baseline_triangles": 10940,
            "refined_nodes": len(read_csv(root / "sentaurus_exports" / "mesh" / "nodes.csv")),
            "refined_triangles": len(read_csv(root / "sentaurus_exports" / "mesh" / "elements.csv")),
        },
        "sdevice_plt_audit": plt_audit,
        "sdevice_mesh_sensitivity_max_log_error_decade": {
            key: maximum_log_error(refined_sentaurus, baseline_sentaurus, key)
            for key in ("Ic_A_per_um", "Ib_A_per_um", "Ie_A_per_um", "beta_abs")
        },
        "refined_vela_vs_sdevice_max_log_error_decade": {
            key: maximum_log_error(refined_vela, refined_sentaurus, key)
            for key in ("Ic_A_per_um", "Ib_A_per_um", "Ie_A_per_um", "beta_abs")
        },
        "vce3": {
            "sdevice": refined_sentaurus[-1],
            "vela": refined_vela[-1],
            "hole_current_p95_log_error_decade": transport["current_density"]["hole"]
            ["log10_magnitude_error"]["p95_absolute_error"],
            "hole_current_normalized_vector_rmse": transport["current_density"]["hole"]
            ["normalized_vector_rmse"],
            "hole_current_cosine_similarity": transport["current_density"]["hole"]
            ["global_vector_cosine_similarity"],
        },
        "coarse_to_refined_vce3": {
            "hole_current_p95_log_error_decade": {
                "coarse": baseline_transport["current_density"]["hole"]["log10_magnitude_error"]
                ["p95_absolute_error"],
                "refined": transport["current_density"]["hole"]["log10_magnitude_error"]
                ["p95_absolute_error"],
            },
            "hole_current_normalized_vector_rmse": {
                "coarse": baseline_transport["current_density"]["hole"]["normalized_vector_rmse"],
                "refined": transport["current_density"]["hole"]["normalized_vector_rmse"],
            },
            "potential_rmse_V": {
                "coarse": baseline_spatial["fields"]["electrostatic_potential"]["gate"]["observed"]["rmse"],
                "refined": refined_spatial["fields"]["electrostatic_potential"]["gate"]["observed"]["rmse"],
            },
            "electron_density_rmse_decade": {
                "coarse": baseline_spatial["fields"]["electron_density"]["gate"]["observed"]["rmse"],
                "refined": refined_spatial["fields"]["electron_density"]["gate"]["observed"]["rmse"],
            },
            "hole_density_rmse_decade": {
                "coarse": baseline_spatial["fields"]["hole_density"]["gate"]["observed"]["rmse"],
                "refined": refined_spatial["fields"]["hole_density"]["gate"]["observed"]["rmse"],
            },
            "hole_density_full_domain_max_decade": {
                "coarse": baseline_spatial["fields"]["hole_density"]["full_domain_characterization"]
                ["maximum_absolute_error"],
                "refined": refined_spatial["fields"]["hole_density"]["full_domain_characterization"]
                ["maximum_absolute_error"],
            },
        },
        "gates": {
            "sdevice_curve_complete_and_kcl": bool(plt_audit["pass"]),
            "refined_spatial_state": bool(refined_spatial["overall_pass"]),
            "refined_transport_sources": bool(transport["overall_pass"]),
        },
        "inputs_sha256": {
            "refined_sde": sha256(FIXTURE / "source" / "bjt_sde_local_refined.cmd"),
            "sdevice_deck": sha256(FIXTURE / "source" / "bjt_m1_des.cmd"),
            "mesh": sha256(root / "vela" / "input" / "mesh.json"),
            "doping": sha256(root / "vela" / "input" / "doping.csv"),
            "materials": sha256(FIXTURE / "vela" / "materials_sentaurus2022.json"),
            "runner": sha256(args.runner.resolve()),
        },
    }
    write_json(root / "comparison" / "mesh_sensitivity_report.json", report)
    metrics = report["sdevice_mesh_sensitivity_max_log_error_decade"]
    cross = report["refined_vela_vs_sdevice_max_log_error_decade"]
    coarse_refined = report["coarse_to_refined_vce3"]
    lines = [
        "# Genius BJT local mesh-refinement sensitivity",
        "",
        f"- Mesh: 5611/10940 -> {report['mesh']['refined_nodes']}/{report['mesh']['refined_triangles']} nodes/triangles.",
        f"- SDevice mesh sensitivity (VCE=0.5-3 V): Ic {metrics['Ic_A_per_um']:.6g}, Ib {metrics['Ib_A_per_um']:.6g}, Ie {metrics['Ie_A_per_um']:.6g}, beta {metrics['beta_abs']:.6g} decade.",
        f"- Refined Vela vs refined SDevice: Ic {cross['Ic_A_per_um']:.6g}, Ib {cross['Ib_A_per_um']:.6g}, Ie {cross['Ie_A_per_um']:.6g}, beta {cross['beta_abs']:.6g} decade.",
        f"- Refined 3 V hole-current P95: {report['vce3']['hole_current_p95_log_error_decade']:.6g} decade.",
        f"- Refined 3 V hole-current normalized vector RMSE/cosine: {report['vce3']['hole_current_normalized_vector_rmse']:.6g} / {report['vce3']['hole_current_cosine_similarity']:.6g}.",
        f"- Hole-current P95 coarse -> refined: {coarse_refined['hole_current_p95_log_error_decade']['coarse']:.6g} -> {coarse_refined['hole_current_p95_log_error_decade']['refined']:.6g} decade.",
        f"- Hole-density masked RMSE coarse -> refined: {coarse_refined['hole_density_rmse_decade']['coarse']:.6g} -> {coarse_refined['hole_density_rmse_decade']['refined']:.6g} decade.",
        f"- Gates: SDevice curve/KCL={report['gates']['sdevice_curve_complete_and_kcl']}, spatial={report['gates']['refined_spatial_state']}, transport/source={report['gates']['refined_transport_sources']}.",
    ]
    (root / "comparison" / "mesh_sensitivity_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage", choices=("prepare", "vela", "spatial", "analyze", "all"), default="all"
    )
    parser.add_argument("--artifact-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--importer", type=Path, default=DEFAULT_IMPORTER)
    parser.add_argument("--python", type=Path, default=DEFAULT_PYTHON)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results = {}
    if args.stage in ("prepare", "all"):
        results["prepare"] = prepare(args)
    if args.stage in ("vela", "all"):
        results["vela"] = vela(args)
    elif args.stage == "spatial":
        results["spatial"] = spatial(args)
    if args.stage in ("analyze", "all"):
        results["analyze"] = analyze(args)
    print(json.dumps(results, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
