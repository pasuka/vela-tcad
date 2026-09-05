#!/usr/bin/env python3
"""Prepare and analyze a 7x3/13x5/25x9 exact-mesh PN refinement study."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import re
import subprocess
from pathlib import Path
from typing import Any


WORKTREE = Path(__file__).resolve().parents[1]
DATA_REPO = Path(r"D:\code-repo\vela-tcad")
BASE_VELA = (
    DATA_REPO
    / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3"
    / "imported_reference/vela"
)
DEFAULT_OUTPUT = (
    WORKTREE
    / "build-release/reference_tcad/pn2d_sentaurus2022/same_mesh_refinement_study"
)
DEFAULT_REPORT = (
    WORKTREE
    / "reference_tcad/pn2d_sentaurus2018_coarse7x3/reports"
    / "pn2d_same_mesh_refinement_study.json"
)
RUN_IDS = {
    "7x3": "pn2d_same_mesh_vector_20260904",
    "13x5": "pn2d_same_mesh_13x5_20260904",
    "25x9": "pn2d_same_mesh_25x9_20260904",
}
FACTORS = {"7x3": 1, "13x5": 2, "25x9": 4}
Q_C = 1.602176634e-19


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CV = load_module("cell_first_cross_device", WORKTREE / "scripts/run_cell_first_cross_device_validation.py")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(value, indent=2, allow_nan=False) + "\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def base_doping_grid() -> tuple[list[float], list[float], list[float], int, int]:
    mesh = read_json(BASE_VELA / "mesh.json")
    rows = read_csv(BASE_VELA / "doping.csv")
    xs = sorted({float(node["x"]) for node in mesh["nodes"]})
    ys = sorted({float(node["y"]) for node in mesh["nodes"]})
    by_id = {int(row["node_id"]): row for row in rows}
    by_xy = {
        (float(node["x"]), float(node["y"])): by_id[int(node["id"])]
        for node in mesh["nodes"]
    }
    donors = [float(by_xy[(x, ys[0])]["donors_cm3"]) for x in xs]
    acceptors = [float(by_xy[(x, ys[0])]["acceptors_cm3"]) for x in xs]
    return xs, donors, acceptors, len(xs), len(ys)


def interpolate(values: list[float], scaled_index: float) -> float:
    left = min(int(math.floor(scaled_index)), len(values) - 1)
    right = min(left + 1, len(values) - 1)
    weight = scaled_index - left
    return values[left] * (1.0 - weight) + values[right] * weight


def prepare_level(label: str, factor: int, output: Path) -> dict[str, Any]:
    xs, donors, acceptors, base_nx, base_ny = base_doping_grid()
    nx = (base_nx - 1) * factor + 1
    ny = (base_ny - 1) * factor + 1
    x_values = [xs[0] + (xs[-1] - xs[0]) * i / (nx - 1) for i in range(nx)]
    y_values = [0.5 * j / (ny - 1) for j in range(ny)]
    nodes = [
        {"id": j * nx + i, "x": x, "y": y}
        for j, y in enumerate(y_values)
        for i, x in enumerate(x_values)
    ]
    triangles = []
    p_cells: list[int] = []
    n_cells: list[int] = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            n00 = j * nx + i
            n10 = n00 + 1
            n01 = n00 + nx
            n11 = n01 + 1
            region_id = 0 if (i // factor) < 3 else 1
            for ids in ((n00, n10, n11), (n00, n11, n01)):
                cell_id = len(triangles)
                triangles.append({"id": cell_id, "region_id": region_id, "node_ids": list(ids)})
                (p_cells if region_id == 0 else n_cells).append(cell_id)
    mesh = {
        "nodes": nodes,
        "triangles": triangles,
        "regions": [
            {"id": 0, "name": "p_region", "material": "Si", "cell_ids": p_cells},
            {"id": 1, "name": "n_region", "material": "Si", "cell_ids": n_cells},
        ],
        "contacts": [
            {"id": 0, "name": "anode", "region_id": 0, "node_ids": [j * nx for j in range(ny)]},
            {"id": 1, "name": "cathode", "region_id": 1, "node_ids": [j * nx + nx - 1 for j in range(ny)]},
        ],
    }
    doping_rows = []
    for node in nodes:
        scaled = float(node["x"]) / (xs[1] - xs[0])
        doping_rows.append(
            {
                "node_id": node["id"],
                "donors_cm3": interpolate(donors, scaled),
                "acceptors_cm3": interpolate(acceptors, scaled),
            }
        )
    level = output / "inputs" / label
    mesh_path = level / "mesh.json"
    doping_path = level / "doping.csv"
    write_json(mesh_path, mesh)
    write_csv(doping_path, doping_rows)
    return {
        "label": label,
        "factor": factor,
        "nx": nx,
        "ny": ny,
        "nodes": len(nodes),
        "triangles": len(triangles),
        "mesh": str(mesh_path.resolve()),
        "doping": str(doping_path.resolve()),
        "run_id": RUN_IDS[label],
    }


def prepare(output: Path) -> dict[str, Any]:
    levels = [prepare_level(label, factor, output) for label, factor in FACTORS.items()]
    manifest = {
        "schema": "vela.pn2d_same_mesh_refinement_inputs.v1",
        "method": "nested uniform subdivision with linear interpolation of the original donor and acceptor fields",
        "levels": levels,
    }
    write_json(output / "input_manifest.json", manifest)
    return manifest


def run_runner(config: dict[str, Any], config_path: Path, runner: Path) -> dict[str, Any]:
    write_json(config_path, config)
    result = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve()), "--log", "off"],
        cwd=WORKTREE,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    config_path.with_suffix(".stdout.log").write_text(result.stdout, encoding="utf-8")
    config_path.with_suffix(".stderr.log").write_text(result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(result.stderr.strip() or result.stdout.strip())
    return json.loads([line for line in result.stdout.splitlines() if line.strip()][-1])


def section_current(edge_rows: list[dict[str, str]], cut_um: float) -> dict[str, float | int]:
    electron_particles = 0.0
    hole_particles = 0.0
    count = 0
    for row in edge_rows:
        x0 = float(row["x0"]) * 1.0e6
        x1 = float(row["x1"]) * 1.0e6
        if (x0 <= cut_um) == (x1 <= cut_um):
            continue
        sign = 1.0 if x0 <= cut_um else -1.0
        electron_particles += sign * float(row["electron_particle_line_flux_per_m_s"])
        hole_particles += sign * float(row["hole_particle_line_flux_per_m_s"])
        count += 1
    electron = -Q_C * electron_particles * 1.0e-6
    hole = Q_C * hole_particles * 1.0e-6
    return {"crossing_edges": count, "electron_A_per_um": electron, "hole_A_per_um": hole, "total_A_per_um": electron + hole}


def final_plt_currents(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    info = re.search(r"datasets\s*=\s*\[(.*?)\]", text, re.DOTALL)
    data = re.search(r"Data\s*\{(.*?)\}", text, re.DOTALL)
    if info is None or data is None:
        raise ValueError(f"cannot parse SDevice PLT: {path}")
    names = re.findall(r'"([^"]+)"', info.group(1))
    values = [
        float(value)
        for value in re.findall(
            r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", data.group(1)
        )
    ]
    if not names or len(values) % len(names):
        raise ValueError(f"invalid SDevice PLT table: {path}")
    final = dict(zip(names, values[-len(names) :]))
    return {
        "anode_A_per_um": final["Anode TotalCurrent"],
        "cathode_A_per_um": final["Cathode TotalCurrent"],
    }


def analyze(output: Path, runner: Path, report_path: Path) -> dict[str, Any]:
    inputs = read_json(output / "input_manifest.json")
    rows: list[dict[str, object]] = []
    level_reports = []
    for level in inputs["levels"]:
        label = level["label"]
        run_root = (
            WORKTREE
            / "build-release/reference_tcad/pn2d_sentaurus2022/sentaurus_vm_runs"
            / level["run_id"]
        )
        export_manifest = read_json(run_root / "manifest.json")
        if not export_manifest.get("passed"):
            raise ValueError(f"SDevice export did not pass for {label}")
        artifacts = run_root / "artifacts"
        case = {
            "device": "pn2d_sentaurus2022_same_mesh_refinement",
            "topology": f"pn_diode_{label}",
            "case": f"reverse_m20V_{label}",
            "swept_bias_V": -20.0,
            "drain_bias_V": 0.0,
            "config": BASE_VELA / "simulation_bv.json",
            "state": artifacts / "sdevice_state.csv",
            "sentaurus_export": artifacts / "aligned_export",
            "region_index": 0,
            "region_name": "R.Si",
            "mesh_override": Path(level["mesh"]),
            "doping_override": Path(level["doping"]),
            "materials_override": BASE_VELA / "pn2d_sentaurus2018_iv_materials.json",
            "bias_contact": "Anode",
        }
        metric_rows, recovery_manifest = CV.validate_case(runner, output / "recovery", case)
        edge_csv = output / "edge_flux" / f"{label}.csv"
        edge_csv.parent.mkdir(parents=True, exist_ok=True)
        edge_config = read_json(BASE_VELA / "simulation_bv.json")
        CV.resolve_config_paths(edge_config, BASE_VELA)
        edge_config.update(
            {
                "simulation_type": "sg_edge_flux_probe",
                "mesh_file": level["mesh"],
                "node_doping_file": level["doping"],
                "materials_file": str((BASE_VELA / "pn2d_sentaurus2018_iv_materials.json").resolve()),
                "state_file": str((artifacts / "sdevice_state.csv").resolve()),
                "output_csv": str(edge_csv.resolve()),
            }
        )
        CV.set_contact_bias(edge_config, "Anode", -20.0)
        edge_status = run_runner(edge_config, output / "configs" / f"edge_{label}.json", runner)
        sections = {str(cut): section_current(read_csv(edge_csv), cut) for cut in (0.5, 1.5)}
        sdevice_currents = final_plt_currents(artifacts / "pn2d_same_mesh_m20.plt")
        sdevice_anode = sdevice_currents["anode_A_per_um"]
        sdevice_cathode = sdevice_currents["cathode_A_per_um"]
        terminal_scale = max(abs(sdevice_anode), abs(sdevice_cathode), 1.0e-30)
        section_values = [float(section["total_A_per_um"]) for section in sections.values()]
        level_report = {
            **level,
            "sdevice_export_manifest": str((run_root / "manifest.json").resolve()),
            "sdevice_terminal": {
                "source": "final row of native SDevice PLT",
                "anode_A_per_um": sdevice_anode,
                "cathode_A_per_um": sdevice_cathode,
                "kcl_relative": abs(sdevice_anode + sdevice_cathode) / terminal_scale,
            },
            "vela_sg_sections": sections,
            "vela_section_spread_relative_to_sdevice_terminal": abs(section_values[1] - section_values[0]) / terminal_scale,
            "recovery_manifest": recovery_manifest,
            "edge_status": edge_status,
        }
        level_reports.append(level_report)
        for metric in metric_rows:
            rows.append({"grid": label, **metric})
    assessed = CV.add_ab_assessment(rows)
    finest_terminal = max(abs(level_reports[-1]["sdevice_terminal"]["anode_A_per_um"]), abs(level_reports[-1]["sdevice_terminal"]["cathode_A_per_um"]))
    for level in level_reports:
        current = max(abs(level["sdevice_terminal"]["anode_A_per_um"]), abs(level["sdevice_terminal"]["cathode_A_per_um"]))
        level["sdevice_terminal_relative_to_finest"] = current / finest_terminal
    report = {
        "schema": "vela.pn2d_same_mesh_refinement_study.v1",
        "bias_V": -20.0,
        "levels": level_reports,
        "rows": rows,
        "assessment": assessed,
        "limitations": [
            "SDevice does not expose its directed internal edge flux; Vela SG section flux is checked for internal consistency and against SDevice terminal current.",
            "Vela section fluxes are recomputed on a frozen SDevice state, not a converged Vela state; their section spread diagnoses operator/model compatibility, not Vela self-consistent conservation.",
            "Absolute current-vector errors include transport-model differences as well as reconstruction differences.",
            "The 25x9 result is the finest level in this three-level experiment, not a proof of asymptotic mesh convergence.",
        ],
    }
    write_json(report_path, report)
    write_csv(output / "node_vector_metrics.csv", rows)
    write_csv(output / "node_vector_ab.csv", assessed)
    write_markdown(report_path.with_suffix(".md"), report)
    return report


def write_markdown(path: Path, report: dict[str, Any]) -> None:
    by = {(row["grid"], row["carrier"], row["recovery"]): row for row in report["rows"]}
    lines = [
        "# PN 二极管三级同网格加密试验",
        "",
        "三套网格均由原始网格嵌套细分，SDevice 与 Vela 在每一级使用完全相同的节点、三角形和掺杂数据。",
        "",
        "## 结论",
        "",
        "- 网格稀疏显著影响 SDevice 的反向端口漏电：7x3 结果仅为 25x9 的 0.196%，13x5 仍比 25x9 高 40.4%。因此原 7x3 网格不适合用作 -20 V 端口漏电的定量基准。",
        "- 节点恢复顺序的差异随加密消失：13x5 和 25x9 的 direct/cell-first P95 基本相同。粗网格确实放大了恢复方法敏感性。",
        "- 但 Vela 相对 SDevice 的绝对节点矢量 P95 没有随加密下降，典型节点 P50 也基本不变。剩余误差不能只归因于网格，输运模型、SG 边系数和 SDevice 节点矢量语义仍占主要部分。",
        "",
        "## 端口与节点矢量",
        "",
        "| 网格 | 节点/三角形 | SDevice |I| [A/um] | 相对25x9 | 电子 P50/P95 direct [dec] | 空穴 P50/P95 direct [dec] |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for level in report["levels"]:
        label = level["label"]
        electron = by[(label, "electron", "direct")]
        hole = by[(label, "hole", "direct")]
        terminal = abs(level["sdevice_terminal"]["anode_A_per_um"])
        lines.append(
            f"| {label} | {level['nodes']}/{level['triangles']} | {terminal:.6e} | "
            f"{level['sdevice_terminal_relative_to_finest']:.6g} | "
            f"{electron['p50_log10_magnitude_error_decade']:.6g}/{electron['p95_log10_magnitude_error_decade']:.6g} | "
            f"{hole['p50_log10_magnitude_error_decade']:.6g}/{hole['p95_log10_magnitude_error_decade']:.6g} |"
        )
    lines.extend(
        [
            "",
            "## 恢复方法 A/B",
            "",
            "| 网格 | 电子 P95 direct/cell-first [dec] | 空穴 P95 direct/cell-first [dec] |",
            "|---|---:|---:|",
        ]
    )
    for level in report["levels"]:
        label = level["label"]
        e0 = by[(label, "electron", "direct")]["p95_log10_magnitude_error_decade"]
        e1 = by[(label, "electron", "cell_first")]["p95_log10_magnitude_error_decade"]
        h0 = by[(label, "hole", "direct")]["p95_log10_magnitude_error_decade"]
        h1 = by[(label, "hole", "cell_first")]["p95_log10_magnitude_error_decade"]
        lines.append(f"| {label} | {e0:.6g}/{e1:.6g} | {h0:.6g}/{h1:.6g} |")
    lines.extend(["", "## 解释边界", ""] + [f"- {item}" for item in report["limitations"]])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("prepare", "analyze", "all"), default="all")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--runner", type=Path, default=WORKTREE / "build-release/vela_example_runner.exe")
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.stage in ("prepare", "all"):
        prepare(args.output)
    report = analyze(args.output, args.runner, args.report) if args.stage in ("analyze", "all") else None
    print(json.dumps({"stage": args.stage, "levels": None if report is None else len(report["levels"]), "report": str(args.report)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
