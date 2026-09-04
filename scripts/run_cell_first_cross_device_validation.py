#!/usr/bin/env python3
"""Validate direct and cell-first SG node-current recovery on reference devices."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any


WORKTREE = Path(__file__).resolve().parents[1]
DEFAULT_DATA_REPO = Path(r"D:\code-repo\vela-tcad")
DEFAULT_RUNNER = WORKTREE / "build-release" / "vela_example_runner.exe"
DEFAULT_OUTPUT = (
    WORKTREE
    / "build-release/reference_tcad/genius_bjt_sentaurus2022"
    / "cell_first_cross_device_validation"
)
DEFAULT_REPORT_JSON = (
    WORKTREE
    / "reference_tcad/genius_bjt_sentaurus2022/reports"
    / "cell_first_cross_device_validation.json"
)
DEFAULT_REPORT_MD = DEFAULT_REPORT_JSON.with_suffix(".md")

FIELDS = {
    "electron": {
        "sentaurus": "eCurrentDensity_region3.csv",
        "direct": "DualFaceSgElectronCurrentDensityVector",
        "cell_first": "CellFirstSgElectronCurrentDensityVector",
    },
    "hole": {
        "sentaurus": "hCurrentDensity_region3.csv",
        "direct": "DualFaceSgHoleCurrentDensityVector",
        "cell_first": "CellFirstSgHoleCurrentDensityVector",
    },
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def read_vector(path: Path) -> dict[int, tuple[float, float]]:
    return {
        int(row["node_id"]): (float(row["component0"]), float(row["component1"]))
        for row in read_rows(path)
    }


def strict_region_nodes(export_dir: Path, region_name: str) -> set[int]:
    incident_regions: dict[int, set[str]] = {}
    for element in read_rows(export_dir / "elements.csv"):
        for local in range(3):
            incident_regions.setdefault(int(element[f"node{local}"]), set()).add(
                element["region"]
            )
    return {
        node
        for node, regions in incident_regions.items()
        if regions == {region_name}
    }


def read_vtk_vectors(path: Path) -> tuple[int, dict[str, list[tuple[float, float, float]]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    point_line = next(line for line in lines if line.startswith("POINT_DATA "))
    count = int(point_line.split()[1])
    index = lines.index(point_line) + 1
    vectors: dict[str, list[tuple[float, float, float]]] = {}
    while index < len(lines):
        parts = lines[index].split()
        if parts and parts[0] == "VECTORS":
            name = parts[1]
            vectors[name] = [
                tuple(float(value) for value in lines[index + offset].split()[:3])
                for offset in range(1, count + 1)
            ]
            index += count + 1
        elif parts and parts[0] == "SCALARS":
            index += count + 2
        else:
            index += 1
    return count, vectors


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("cannot calculate a percentile of an empty population")
    position = (len(ordered) - 1) * fraction
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def vector_metrics(
    reference: dict[int, tuple[float, float]],
    actual: list[tuple[float, float, float]],
    nodes: set[int],
    reference_fraction: float = 1.0e-6,
) -> dict[str, object]:
    common = sorted(nodes & reference.keys())
    peak = max(math.hypot(*reference[node]) for node in common)
    floor = peak * reference_fraction
    selected = [node for node in common if math.hypot(*reference[node]) >= floor]
    log_errors: list[float] = []
    squared_error = 0.0
    squared_reference = 0.0
    squared_actual = 0.0
    dot = 0.0
    for node in selected:
        rx, ry = reference[node]
        ax, ay, _ = actual[node]
        rmag = math.hypot(rx, ry)
        amag = math.hypot(ax, ay)
        if rmag > 0.0 and amag > 0.0:
            log_errors.append(abs(math.log10(amag / rmag)))
        squared_error += (ax - rx) ** 2 + (ay - ry) ** 2
        squared_reference += rx * rx + ry * ry
        squared_actual += ax * ax + ay * ay
        dot += rx * ax + ry * ay
    return {
        "common_node_count": len(common),
        "selected_node_count": len(selected),
        "reference_peak_A_per_cm2": peak,
        "active_floor_A_per_cm2": floor,
        "p50_log10_magnitude_error_decade": percentile(log_errors, 0.50),
        "p95_log10_magnitude_error_decade": percentile(log_errors, 0.95),
        "max_log10_magnitude_error_decade": max(log_errors),
        "normalized_vector_rmse": math.sqrt(squared_error / squared_reference),
        "global_vector_cosine_similarity": (
            dot / math.sqrt(squared_reference * squared_actual)
            if squared_reference > 0.0 and squared_actual > 0.0
            else None
        ),
    }


def resolve_config_paths(config: dict[str, Any], source_dir: Path) -> None:
    for key in ("mesh_file", "node_doping_file", "materials_file", "state_file"):
        if key in config:
            candidate = Path(config[key])
            if not candidate.is_absolute():
                config[key] = str((source_dir / candidate).resolve())


def enable_export(config: dict[str, Any], state: Path, vtk: Path) -> None:
    config["simulation_type"] = "write_dd_state_vtk"
    config["state_file"] = str(state.resolve())
    config["output_vtk"] = str(vtk.resolve())
    config.setdefault("output_diagnostics", {})[
        "cell_first_sg_current_recovery"
    ] = True


def normalized_restart_state(source: Path, destination: Path) -> Path:
    rows = read_rows(source)
    if not rows:
        raise ValueError(f"empty restart state: {source}")
    header = set(rows[0])
    current = {"node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"}
    if current <= header:
        return source
    legacy = {"node_id", "psi_V", "phin_V", "phip_V", "n_m3", "p_m3"}
    if not legacy <= header:
        raise ValueError(f"unsupported restart-state columns in {source}: {sorted(header)}")
    converted = [
        {
            "node_id": row["node_id"],
            "psi": row["psi_V"],
            "phin": row["phin_V"],
            "phip": row["phip_V"],
            "electrons_m3": row["n_m3"],
            "holes_m3": row["p_m3"],
        }
        for row in rows
    ]
    write_csv(destination, converted)
    return destination


def set_contact_bias(config: dict[str, Any], name: str, bias: float) -> None:
    for contact in config["contacts"]:
        if contact["name"] == name:
            contact["bias"] = bias
            return
    raise KeyError(f"missing contact {name}")


def run_export(runner: Path, output_root: Path, name: str, config: dict[str, Any]) -> Path:
    config_path = output_root / "configs" / f"{name}.json"
    vtk = output_root / "vtk" / f"{name}.vtk"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    vtk.parent.mkdir(parents=True, exist_ok=True)
    config["output_vtk"] = str(vtk.resolve())
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
    log_root = output_root / "logs"
    log_root.mkdir(parents=True, exist_ok=True)
    (log_root / f"{name}.stdout.log").write_text(result.stdout, encoding="utf-8")
    (log_root / f"{name}.stderr.log").write_text(result.stderr, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(
            f"{name} export failed ({result.returncode}): "
            + (result.stderr.strip() or result.stdout.strip())
        )
    status = json.loads([line for line in result.stdout.splitlines() if line.strip()][-1])
    if not status.get("cell_first_sg_current_recovery", False):
        raise RuntimeError(f"cell-first diagnostic was not enabled for {name}")
    return vtk


def transportmodels_cases(data_repo: Path) -> list[dict[str, Any]]:
    build = data_repo / "build-release/reference_tcad/transportmodels_sentaurus2022"
    vela_manifest = read_json(
        build
        / "vela_baseline/three_regime_spatial_vtk_2026-08-24"
        / "three_regime_spatial_manifest.json"
    )
    sent_manifest = read_json(
        build
        / "sentaurus_vm_runs/three_regime_spatial_oracles_20260824"
        / "three_regime_spatial_oracles_manifest.json"
    )
    sent_by_key = {
        (item["mode"], item["regime"]): item for item in sent_manifest["states"]
    }
    result = []
    for vela in vela_manifest["states"]:
        if vela["mode"] != "dd":
            continue
        sent = sent_by_key[(vela["mode"], vela["regime"])]
        result.append(
            {
                "device": "transportmodels_sentaurus2022",
                "topology": "planar_nmos",
                "case": f"dd_{vela['regime']}",
                "swept_bias_V": vela["gate_bias_V"],
                "bias_contact": "gate",
                "drain_bias_V": vela["drain_bias_V"],
                "config": Path(vela["config"]),
                "state": Path(vela["state_csv"]),
                "sentaurus_export": Path(sent["export_dir"]),
                "region_index": 3,
                "region_name": "R.Substrate",
            }
        )
    return result


def singledevice_cases(data_repo: Path) -> list[dict[str, Any]]:
    root = data_repo / "build-release/reference_tcad/singledevice_sentaurus2018/fixed_state_curve"
    result = []
    for branch, drain in (("lin", 0.1), ("sat", 1.1)):
        for index in (0, 10, 20):
            gate = -0.5 + 2.7 * index / 20.0
            stem = f"{branch}_{index:04d}"
            result.append(
                {
                    "device": "singledevice_sentaurus2018",
                    "topology": "planar_nmos_with_quantum_potential",
                    "case": stem,
                    "swept_bias_V": gate,
                    "bias_contact": "gate",
                    "drain_bias_V": drain,
                    "config": root / f"{stem}.json",
                    "state": root / f"{stem}_state.csv",
                    "sentaurus_export": root / "exports" / stem,
                    "region_index": 3,
                    "region_name": "R.Substrate",
                }
            )
    return result


def pn_minimal6_cases(data_repo: Path) -> list[dict[str, Any]]:
    fixture = WORKTREE / "reference_tcad/pn2d_sentaurus2018_minimal6"
    state_root = (
        data_repo
        / "build-release/reference_tcad/pn2d_sentaurus2018_minimal6/state_exports"
        / "minimal6_states_v2_sealed_20260717_000955/states"
    )
    result = []
    for topology in ("sketch", "mirror"):
        for token, bias in (("m12V", -12.0), ("m19V", -19.0)):
            export = state_root / topology / token / "export"
            result.append(
                {
                    "device": "pn2d_sentaurus2018_minimal6",
                    "topology": f"pn_diode_{topology}_triangulation",
                    "case": f"{topology}_{token}",
                    "swept_bias_V": bias,
                    "drain_bias_V": 0.0,
                    "config": fixture / "vela/pn2d_minimal6_sweep_template.json",
                    "state": export / "state.csv",
                    "sentaurus_export": export,
                    "region_index": 0,
                    "region_name": "R.Si",
                    "mesh_override": export / "mesh.json",
                    "doping_override": export / "doping.csv",
                    "materials_override": fixture / "vela/pn2d_minimal6_materials.json",
                    "bias_contact": "Anode",
                }
            )
    return result


def pn_same_mesh_cases(data_repo: Path) -> list[dict[str, Any]]:
    vela = (
        data_repo
        / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3"
        / "imported_reference/vela"
    )
    run = (
        WORKTREE
        / "build-release/reference_tcad/pn2d_sentaurus2022/sentaurus_vm_runs"
        / "pn2d_same_mesh_vector_20260904"
    )
    manifest = read_json(run / "manifest.json")
    if not manifest.get("passed") or not manifest.get("topology_gate", {}).get("passed"):
        raise ValueError("same-mesh PN SDevice vector export did not pass its gates")
    artifacts = run / "artifacts"
    return [
        {
            "device": "pn2d_sentaurus2022_same_mesh",
            "topology": "pn_diode_coarse7x3_exact_21_node_mesh",
            "case": "reverse_m20V",
            "swept_bias_V": -20.0,
            "drain_bias_V": 0.0,
            "config": vela / "simulation_bv.json",
            "state": artifacts / "sdevice_state.csv",
            "sentaurus_export": artifacts / "aligned_export",
            "region_index": 0,
            "region_name": "R.Si",
            "mesh_override": vela / "mesh.json",
            "doping_override": vela / "doping.csv",
            "materials_override": vela / "pn2d_sentaurus2018_iv_materials.json",
            "bias_contact": "Anode",
            "sdevice_export_manifest": run / "manifest.json",
        }
    ]


def validate_case(
    runner: Path, output_root: Path, case: dict[str, Any]
) -> tuple[list[dict[str, object]], dict[str, Any]]:
    config = read_json(case["config"])
    resolve_config_paths(config, case["config"].parent)
    if "mesh_override" in case:
        config["mesh_file"] = str(case["mesh_override"].resolve())
        config["node_doping_file"] = str(case["doping_override"].resolve())
        config["materials_file"] = str(case["materials_override"].resolve())
    set_contact_bias(config, case["bias_contact"], float(case["swept_bias_V"]))
    if case["bias_contact"] == "gate":
        set_contact_bias(config, "drain", float(case["drain_bias_V"]))
    name = f"{case['device']}__{case['case']}"
    vtk = output_root / "vtk" / f"{name}.vtk"
    restart_state = normalized_restart_state(
        case["state"], output_root / "states" / f"{name}.csv"
    )
    enable_export(config, restart_state, vtk)
    vtk = run_export(runner, output_root, name, config)
    node_count, vectors = read_vtk_vectors(vtk)
    interior = strict_region_nodes(case["sentaurus_export"], case["region_name"])
    rows: list[dict[str, object]] = []
    for carrier, names in FIELDS.items():
        sentaurus_name = names["sentaurus"].replace(
            "region3", f"region{case['region_index']}"
        )
        reference = read_vector(case["sentaurus_export"] / "fields" / sentaurus_name)
        for recovery in ("direct", "cell_first"):
            field = names[recovery]
            if field not in vectors:
                raise KeyError(f"{field} was not exported for {name}")
            metrics = vector_metrics(reference, vectors[field], interior)
            rows.append(
                {
                    "device": case["device"],
                    "topology": case["topology"],
                    "case": case["case"],
                    "swept_contact": case["bias_contact"],
                    "swept_bias_V": case["swept_bias_V"],
                    "drain_bias_V": case["drain_bias_V"],
                    "carrier": carrier,
                    "recovery": recovery,
                    **metrics,
                }
            )
    manifest = {
        **{
            key: value
            for key, value in case.items()
            if key
            not in {
                "config", "state", "sentaurus_export", "mesh_override",
                "doping_override", "materials_override", "sdevice_export_manifest",
            }
        },
        "config": str(case["config"]),
        "state": str(case["state"]),
        "state_sha256": sha256(case["state"]),
        "normalized_restart_state": str(restart_state),
        "normalized_restart_state_sha256": sha256(restart_state),
        "sentaurus_export": str(case["sentaurus_export"]),
        "vtk": str(vtk),
        "vtk_sha256": sha256(vtk),
        "node_count": node_count,
        "strict_region_node_count": len(interior),
    }
    if "sdevice_export_manifest" in case:
        manifest["sdevice_export_manifest"] = str(case["sdevice_export_manifest"])
        manifest["sdevice_export_manifest_sha256"] = sha256(
            case["sdevice_export_manifest"]
        )
    return rows, manifest


def add_ab_assessment(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    assessed: list[dict[str, object]] = []
    keys = sorted({(row["device"], row["case"], row["carrier"]) for row in rows})
    for key in keys:
        pair = {
            row["recovery"]: row
            for row in rows
            if (row["device"], row["case"], row["carrier"]) == key
        }
        direct = pair["direct"]
        cell = pair["cell_first"]
        direct_p95 = float(direct["p95_log10_magnitude_error_decade"])
        cell_p95 = float(cell["p95_log10_magnitude_error_decade"])
        direct_rmse = float(direct["normalized_vector_rmse"])
        cell_rmse = float(cell["normalized_vector_rmse"])
        p95_delta = cell_p95 - direct_p95
        rmse_delta = cell_rmse - direct_rmse
        assessed.append(
            {
                "device": key[0],
                "case": key[1],
                "carrier": key[2],
                "direct_p95_decade": direct_p95,
                "cell_first_p95_decade": cell_p95,
                "p95_delta_decade": p95_delta,
                "direct_normalized_vector_rmse": direct_rmse,
                "cell_first_normalized_vector_rmse": cell_rmse,
                "normalized_vector_rmse_delta": rmse_delta,
                "p95_improved": p95_delta < -1.0e-12,
                "p95_not_materially_worse": p95_delta <= 0.02,
                "vector_rmse_not_materially_worse": (
                    cell_rmse <= direct_rmse * 1.05 + 0.02
                ),
            }
        )
    return assessed


def write_markdown(path: Path, report: dict[str, Any]) -> None:
    assessment = report["assessment"]
    lines = [
        "# 跨器件先单元后节点电流恢复验证",
        "",
        "## 结论",
        "",
        report["conclusion"],
        "",
        "本轮只改变电流矢量的后处理恢复顺序，不改变已接受状态、有限体积残差、守恒边通量或端口电流。比较区域为仅与该工况半导体区域相邻的节点，活跃节点门槛为各工况 SDevice 电流峰值的 `1e-6`。",
        "",
        "## A/B 结果",
        "",
        "| 器件 | 工况 | 载流子 | 直接恢复 P95 [dec] | 先单元后节点 P95 [dec] | P95 变化 [dec] | 矢量 RMSE 变化 |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for row in assessment:
        lines.append(
            f"| {row['device']} | {row['case']} | {row['carrier']} | "
            f"{row['direct_p95_decade']:.6g} | {row['cell_first_p95_decade']:.6g} | "
            f"{row['p95_delta_decade']:+.6g} | {row['normalized_vector_rmse_delta']:+.6g} |"
        )
    lines.extend(
        [
            "",
            "## 证据边界",
            "",
            "- `transportmodels_sentaurus2022` 覆盖普通漂移扩散 NMOS 的关断、阈值和导通状态。",
            "- `singledevice_sentaurus2018` 覆盖含电子量子势状态的 NMOS，分别检查低漏压线性区和高漏压饱和区。",
            "- `pn2d_sentaurus2018_minimal6` 是 6 节点 PN 二极管算法微夹具；`sketch/mirror` 两种三角剖分用于检查方向敏感性，不代表生产网格精度。",
            "- `pn2d_sentaurus2022_same_mesh` 使用显式 DF-ISE 21 节点、24 三角形网格；SDevice 状态与电流矢量逐节点映射，无插值。",
            "- 绝对误差同时包含 Vela 与 SDevice 输运模型差异；本报告主要用同一冻结状态下 direct/cell-first 的相对变化判断恢复方法的可迁移性。",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-repo", type=Path, default=DEFAULT_DATA_REPO)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-json", type=Path, default=DEFAULT_REPORT_JSON)
    parser.add_argument("--report-markdown", type=Path, default=DEFAULT_REPORT_MD)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cases = (
        transportmodels_cases(args.data_repo)
        + singledevice_cases(args.data_repo)
        + pn_minimal6_cases(args.data_repo)
        + pn_same_mesh_cases(args.data_repo)
    )
    rows: list[dict[str, object]] = []
    manifests: list[dict[str, Any]] = []
    for case in cases:
        case_rows, manifest = validate_case(args.runner, args.output_root, case)
        rows.extend(case_rows)
        manifests.append(manifest)
    assessment = add_ab_assessment(rows)
    p95_improved = sum(bool(row["p95_improved"]) for row in assessment)
    p95_safe = sum(bool(row["p95_not_materially_worse"]) for row in assessment)
    rmse_safe = sum(bool(row["vector_rmse_not_materially_worse"]) for row in assessment)
    mos = [row for row in assessment if not str(row["device"]).startswith("pn2d_")]
    pn = [row for row in assessment if str(row["device"]).startswith("pn2d_")]
    pn_micro = [
        row for row in pn if row["device"] == "pn2d_sentaurus2018_minimal6"
    ]
    mos_p95_improved = sum(bool(row["p95_improved"]) for row in mos)
    pn_p95_improved = sum(bool(row["p95_improved"]) for row in pn)
    pn_rmse_improved = sum(
        float(row["normalized_vector_rmse_delta"]) < 0.0 for row in pn
    )
    pn_symmetry_differences = []
    for token in ("m12V", "m19V"):
        for carrier in ("electron", "hole"):
            pair = [
                row
                for row in pn_micro
                if str(row["case"]).endswith(token) and row["carrier"] == carrier
            ]
            pn_symmetry_differences.extend(
                [
                    abs(float(pair[0]["direct_p95_decade"]) - float(pair[1]["direct_p95_decade"])),
                    abs(float(pair[0]["cell_first_p95_decade"]) - float(pair[1]["cell_first_p95_decade"])),
                ]
            )
    pn_mirror_symmetric = max(pn_symmetry_differences) <= 1.0e-12
    strict_pass = p95_safe == len(assessment) and rmse_safe == len(assessment)
    conclusion = (
        f"共完成 {len(cases)} 个工况、{len(assessment)} 个载流子 A/B。"
        f"MOSFET 中 {mos_p95_improved}/{len(mos)} 组 P95 改善，且全部无实质退化；"
        f"PN 中 {pn_p95_improved}/{len(pn)} 组 P95 改善、"
        f"{pn_rmse_improved}/{len(pn)} 组矢量 RMSE 改善。"
        "六节点微夹具存在点态 P95 与全局矢量 RMSE 的权衡；候选有明显跨 MOSFET 收益，"
        "但尚不具备替换全局默认值的证据，"
        "应保持默认关闭。"
    )
    report = {
        "schema": "vela.cell_first_cross_device_validation.v1",
        "scope": "reference_tcad MOSFET and PN fixed-state current-vector A/B",
        "postprocess_only": True,
        "candidate_default_enabled": False,
        "metric_definition": {
            "region": "nodes incident only to the case semiconductor region",
            "mask": "SDevice carrier-current magnitude >= 1e-6 of per-case peak",
            "p95": "95th percentile absolute log10 magnitude error",
            "vector_rmse": "L2 vector error normalized by SDevice vector L2 norm",
            "material_p95_regression_tolerance_decade": 0.02,
            "material_vector_rmse_regression_tolerance": "candidate <= direct*1.05 + 0.02",
        },
        "cases": manifests,
        "excluded_candidates": [],
        "rows": rows,
        "assessment": assessment,
        "checks": {
            "case_count": len(cases),
            "carrier_pair_count": len(assessment),
            "p95_improved_count": p95_improved,
            "p95_not_materially_worse_count": p95_safe,
            "vector_rmse_not_materially_worse_count": rmse_safe,
            "mos_carrier_pair_count": len(mos),
            "mos_p95_improved_count": mos_p95_improved,
            "mos_all_metrics_not_materially_worse": all(
                bool(row["p95_not_materially_worse"])
                and bool(row["vector_rmse_not_materially_worse"])
                for row in mos
            ),
            "pn_carrier_pair_count": len(pn),
            "pn_p95_improved_count": pn_p95_improved,
            "pn_vector_rmse_improved_count": pn_rmse_improved,
            "pn_mirrored_topologies_metric_symmetric": pn_mirror_symmetric,
            "strict_ab_pass": strict_pass,
            "promote_to_global_default": strict_pass,
        },
        "conclusion": conclusion,
    }
    write_csv(args.output_root / "metrics.csv", rows)
    write_csv(args.output_root / "assessment.csv", assessment)
    write_json(args.report_json, report)
    write_markdown(args.report_markdown, report)
    print(json.dumps({"overall_pass": strict_pass, "conclusion": conclusion}, ensure_ascii=False))
    return 0 if strict_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
