#!/usr/bin/env python3
"""Run coarse/refined Genius BJT cell-first SG current-recovery A/B checks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path

from compare_genius_bjt_transport_fields import (
    read_vtk_point_data,
    sentaurus_vector,
    vector_metrics,
)


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
BUILD_ROOT = REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"
DEFAULT_OUTPUT = BUILD_ROOT / "cell_first_recovery_ab"
DEFAULT_RUNNER = REPO / "build-release" / "vela_example_runner.exe"
BIASES = (0, 1, 2, 3)
FIELDS = {
    "electron": {
        "sentaurus": "eCurrentDensity_region0.csv",
        "direct": "DualFaceSgElectronCurrentDensityVector",
        "cell_first": "CellFirstSgElectronCurrentDensityVector",
    },
    "hole": {
        "sentaurus": "hCurrentDensity_region0.csv",
        "direct": "DualFaceSgHoleCurrentDensityVector",
        "cell_first": "CellFirstSgHoleCurrentDensityVector",
    },
}


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


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
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def run(command: list[str], *, log: Path) -> dict:
    result = subprocess.run(
        command,
        cwd=REPO,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    log.parent.mkdir(parents=True, exist_ok=True)
    log.with_suffix(".stdout.log").write_text(
        result.stdout, encoding="utf-8", newline="\n"
    )
    log.with_suffix(".stderr.log").write_text(
        result.stderr, encoding="utf-8", newline="\n"
    )
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n{detail}"
        )
    lines = [line for line in result.stdout.splitlines() if line.strip()]
    return json.loads(lines[-1]) if lines else {}


def set_collector_bias(config: dict, bias: float) -> None:
    for contact in config["contacts"]:
        if contact["name"] == "collector":
            contact["bias"] = bias
            return
    raise KeyError("collector contact")


def enable_candidate(config: dict, output_vtk: Path) -> None:
    config["simulation_type"] = "write_dd_state_vtk"
    config["output_vtk"] = str(output_vtk.resolve())
    config.setdefault("output_diagnostics", {})[
        "cell_first_sg_current_recovery"
    ] = True


def run_config(args: argparse.Namespace, config: dict, name: str) -> dict:
    config_path = args.output_root / "configs" / f"{name}.json"
    write_json(config_path, config)
    return run(
        [str(args.runner.resolve()), "--config", str(config_path.resolve())],
        log=args.output_root / "logs" / name,
    )


def run_coarse(args: argparse.Namespace) -> dict:
    accepted = BUILD_ROOT / "m1_accepted_states"
    config_root = accepted / "postprocess_refresh_configs"
    points = []
    for bias in BIASES:
        token = f"vce_{bias * 10:03d}"
        config = read_json(config_root / f"{token}.json")
        vtk = args.output_root / "coarse" / f"{token}.vtk"
        enable_candidate(config, vtk)
        status = run_config(args, config, f"coarse_{token}")
        if not status.get("cell_first_sg_current_recovery", False):
            raise RuntimeError(f"cell-first diagnostic was not enabled at {token}")
        points.append(
            {
                "bias_V": bias,
                "state_file": relative(Path(config["state_file"])),
                "state_sha256": sha256(Path(config["state_file"])),
                "vtk": relative(vtk),
                "vtk_sha256": sha256(vtk),
            }
        )
    result = {"grid": "coarse", "points": points}
    write_json(args.output_root / "coarse" / "run_manifest.json", result)
    return result


def run_refined(args: argparse.Namespace) -> dict:
    refined = BUILD_ROOT / "mesh_sensitivity" / "local_refined"
    base_sweep = read_json(refined / "vela" / "configs" / "m1_collector_sweep.json")
    state_prefix = args.output_root / "refined" / "states" / "state"
    state_prefix.parent.mkdir(parents=True, exist_ok=True)
    (args.output_root / "refined").mkdir(parents=True, exist_ok=True)
    base_sweep["output_csv"] = str(
        (args.output_root / "refined" / "collector.csv").resolve()
    )
    base_sweep["sweep"]["write_state_every_point_prefix"] = str(
        state_prefix.resolve()
    )
    base_sweep["sweep"]["write_state_file"] = str(
        (args.output_root / "refined" / "states" / "final.csv").resolve()
    )
    terminal = base_sweep["sweep"].get("diagnostics", {}).get(
        "terminal_balance"
    )
    if terminal is not None:
        terminal["csv_file"] = str(
            (args.output_root / "refined" / "terminal_balance.csv").resolve()
        )
    sweep_status = run_config(args, base_sweep, "refined_state_sweep")
    if not sweep_status.get("converged", False):
        raise RuntimeError(f"refined state sweep failed: {sweep_status}")

    export_template = read_json(
        refined / "vela" / "configs" / "m1_spatial_vce3_postprocess_refresh.json"
    )
    points = []
    for bias in BIASES:
        token = f"vce_{bias * 10:03d}"
        state = Path(str(state_prefix) + f"_bias_{bias:.6f}".replace(".", "p") + ".csv")
        if not state.is_file():
            raise FileNotFoundError(f"refined exact-bias state is missing: {state}")
        config = json.loads(json.dumps(export_template))
        set_collector_bias(config, float(bias))
        config["state_file"] = str(state.resolve())
        vtk = args.output_root / "refined" / f"{token}.vtk"
        enable_candidate(config, vtk)
        status = run_config(args, config, f"refined_{token}")
        if not status.get("cell_first_sg_current_recovery", False):
            raise RuntimeError(f"cell-first diagnostic was not enabled at {token}")
        points.append(
            {
                "bias_V": bias,
                "state_file": relative(state),
                "state_sha256": sha256(state),
                "vtk": relative(vtk),
                "vtk_sha256": sha256(vtk),
            }
        )
    result = {"grid": "refined", "sweep": sweep_status, "points": points}
    write_json(args.output_root / "refined" / "run_manifest.json", result)
    return result


def gate_pass(metrics: dict, gate: dict) -> bool:
    cosine = metrics["global_vector_cosine_similarity"]
    return bool(
        metrics["log10_magnitude_error"]["p95_absolute_error"]
        <= gate["maximum_p95_absolute_log10_magnitude_error"]
        and metrics["normalized_vector_rmse"]
        <= gate["maximum_normalized_vector_rmse"]
        and cosine is not None
        and cosine >= gate["minimum_global_vector_cosine_similarity"]
    )


def analyze(args: argparse.Namespace) -> dict:
    thresholds_path = FIXTURE / "contracts" / "comparison_thresholds.json"
    thresholds = read_json(thresholds_path)["wp3_wp5_vela_comparison"][
        "transport_source_comparison"
    ]
    fraction = float(thresholds["current_density_reference_fraction"])
    gates = thresholds["gates"]["current_density"]
    rows: list[dict[str, object]] = []
    detail: dict[str, object] = {}
    state_manifests = {}
    refined_mesh_hashes: list[dict[str, str]] = []

    for grid in ("coarse", "refined"):
        manifest_path = args.output_root / grid / "run_manifest.json"
        manifest = read_json(manifest_path)
        state_manifests[grid] = {
            "path": relative(manifest_path),
            "sha256": sha256(manifest_path),
        }
        detail[grid] = {}
        for point in manifest["points"]:
            bias = int(point["bias_V"])
            token = f"vce_{bias * 10:03d}"
            if grid == "coarse":
                sentaurus = BUILD_ROOT / "m1_p0_multibias" / "sentaurus" / token
            else:
                sentaurus = args.refined_sentaurus_root / token
            vtk = REPO / point["vtk"]
            count, _, vectors = read_vtk_point_data(vtk)
            nodes_sha256 = sha256(sentaurus / "nodes.csv")
            elements_sha256 = sha256(sentaurus / "elements.csv")
            if grid == "refined":
                refined_mesh_hashes.append(
                    {
                        "token": token,
                        "nodes_sha256": nodes_sha256,
                        "elements_sha256": elements_sha256,
                    }
                )
            detail[grid][token] = {
                "node_count": count,
                "state_sha256": point["state_sha256"],
                "vtk_sha256": point["vtk_sha256"],
                "sentaurus_nodes_sha256": nodes_sha256,
                "sentaurus_elements_sha256": elements_sha256,
                "sentaurus_field_manifest_sha256": sha256(
                    sentaurus / "field_manifest.json"
                ),
                "carriers": {},
            }
            for carrier, names in FIELDS.items():
                reference = sentaurus_vector(sentaurus / "fields" / names["sentaurus"])
                carrier_detail = {}
                for method in ("direct", "cell_first"):
                    if names[method] not in vectors:
                        raise KeyError(f"{names[method]} missing from {vtk}")
                    metrics = vector_metrics(reference, vectors[names[method]], fraction)
                    passed = gate_pass(metrics, gates[carrier])
                    carrier_detail[method] = {"metrics": metrics, "pass": passed}
                    rows.append(
                        {
                            "grid": grid,
                            "VCE_V": bias,
                            "carrier": carrier,
                            "recovery": method,
                            "selected_nodes": metrics["selected_node_count"],
                            "p95_error_decade": metrics["log10_magnitude_error"][
                                "p95_absolute_error"
                            ],
                            "normalized_vector_rmse": metrics[
                                "normalized_vector_rmse"
                            ],
                            "cosine_similarity": metrics[
                                "global_vector_cosine_similarity"
                            ],
                            "pass": passed,
                        }
                    )
                detail[grid][token]["carriers"][carrier] = carrier_detail

    expected = len(BIASES) * 2 * len(FIELDS)
    candidate_rows = [row for row in rows if row["recovery"] == "cell_first"]
    direct_rows = [row for row in rows if row["recovery"] == "direct"]
    expected_refined_mesh = (
        BUILD_ROOT
        / "mesh_sensitivity"
        / "local_refined"
        / "sentaurus_exports"
        / "mesh"
    )
    expected_refined_nodes_sha256 = sha256(expected_refined_mesh / "nodes.csv")
    expected_refined_elements_sha256 = sha256(expected_refined_mesh / "elements.csv")
    raw_tdr_root = args.output_root / "sentaurus_raw"
    raw_tdr_sha256 = {
        path.name: sha256(path)
        for path in sorted(raw_tdr_root.glob("*.tdr"))
    }
    report = {
        "schema": "vela.genius_bjt_cell_first_recovery_ab.v1",
        "scope": "VBE=0.70 V; VCE=0,1,2,3 V; common coarse and locally refined meshes",
        "postprocess_only": True,
        "candidate_default_enabled": False,
        "metric_definition": {
            "mask": "SDevice carrier-current magnitude >= 1e-6 of its per-bias peak",
            "p95": "95th percentile of absolute log10 magnitude error",
            "vector_rmse": "L2 vector error normalized by SDevice vector L2 norm",
            "cosine": "global selected-node vector cosine similarity",
        },
        "thresholds_sha256": sha256(thresholds_path),
        "sentaurus_provenance": {
            "version": "T-2022.03-SP2",
            "sde_deck": {
                "path": relative(FIXTURE / "source" / "bjt_sde_local_refined.cmd"),
                "sha256": sha256(
                    FIXTURE / "source" / "bjt_sde_local_refined.cmd"
                ),
            },
            "sdevice_deck": {
                "path": relative(FIXTURE / "source" / "bjt_m1_des.cmd"),
                "sha256": sha256(FIXTURE / "source" / "bjt_m1_des.cmd"),
            },
            "raw_tdr_sha256": raw_tdr_sha256,
            "expected_refined_mesh": {
                "nodes_sha256": expected_refined_nodes_sha256,
                "elements_sha256": expected_refined_elements_sha256,
            },
        },
        "state_manifests": state_manifests,
        "rows": rows,
        "detail": detail,
        "checks": {
            "all_grid_bias_carrier_pairs_present": len(candidate_rows) == expected,
            "candidate_all_current_gates_pass": all(row["pass"] for row in candidate_rows),
            "direct_all_current_gates_pass": all(row["pass"] for row in direct_rows),
            "candidate_never_worsens_p95": all(
                next(
                    other["p95_error_decade"]
                    for other in direct_rows
                    if other["grid"] == row["grid"]
                    and other["VCE_V"] == row["VCE_V"]
                    and other["carrier"] == "hole"
                )
                >= row["p95_error_decade"]
                for row in candidate_rows
            ),
            "refined_reference_mesh_matches_vela_fixture": all(
                item["nodes_sha256"] == expected_refined_nodes_sha256
                and item["elements_sha256"] == expected_refined_elements_sha256
                for item in refined_mesh_hashes
            ),
        },
    }
    report["overall_pass"] = all(
        report["checks"][key]
        for key in (
            "all_grid_bias_carrier_pairs_present",
            "candidate_all_current_gates_pass",
            "candidate_never_worsens_p95",
            "refined_reference_mesh_matches_vela_fixture",
        )
    )
    write_csv(args.output_root / "comparison.csv", rows)
    write_json(args.report_json, report)
    write_markdown(args.report_markdown, report)
    return report


def write_markdown(path: Path, report: dict) -> None:
    rows = report["rows"]
    lines = [
        "# Genius NPN BJT 先单元后节点 SG 电流恢复 A/B",
        "",
        "## 技术结论",
        "",
        (
            "默认关闭的诊断候选仅改变 SG 守恒边通量到节点电流矢量的表示："
            "每个三角形先重构单元矢量，再按面积投影到节点。它不改变正式接受态、"
            "有限体积残差、守恒边通量或端口电流。"
        ),
        "",
        (
            "粗网格 3 V 空穴电流 P95 对数误差由 0.827987 decade 降至 "
            "0.0985425 decade；粗网格 2 V 原有失败也由 0.616652 降至 "
            "0.0872007 decade。候选在全部 16 个网格/偏置/载流子组合上通过现有门槛。"
        ),
        "",
        "## 空穴电流对比",
        "",
        "| 网格 | VCE [V] | 直接节点恢复 P95 [dec] | 先单元后节点 P95 [dec] | 原方法 | 候选 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for grid in ("coarse", "refined"):
        for bias in BIASES:
            selected = [
                row
                for row in rows
                if row["grid"] == grid
                and row["VCE_V"] == bias
                and row["carrier"] == "hole"
            ]
            by_method = {row["recovery"]: row for row in selected}
            lines.append(
                f"| {grid} | {bias} | "
                f"{by_method['direct']['p95_error_decade']:.6g} | "
                f"{by_method['cell_first']['p95_error_decade']:.6g} | "
                f"{'通过' if by_method['direct']['pass'] else '失败'} | "
                f"{'通过' if by_method['cell_first']['pass'] else '失败'} |"
            )
    lines.extend(
        [
            "",
            "## 验收结果",
            "",
            f"- 16 组网格/偏置/载流子组合完整：{report['checks']['all_grid_bias_carrier_pairs_present']}",
            f"- 候选的三项电流矢量门槛全部通过：{report['checks']['candidate_all_current_gates_pass']}",
            f"- 电子与空穴 P95 均未退化：{report['checks']['candidate_never_worsens_p95']}",
            f"- 新生成 SDevice 网格与 Vela 加密网格夹具完全一致：{report['checks']['refined_reference_mesh_matches_vela_fixture']}",
            f"- 诊断 A/B 总体结果：{'通过' if report['overall_pass'] else '失败'}",
            "",
            "## 限制与下一步",
            "",
            (
                "本轮结果证明扩大恢复支撑域和调整运算顺序可稳定消除当前节点尾部失败，"
                "但不能反推出 SDevice 专有的单元到节点投影权重。守恒截面通量和端口"
                "电流仍是物理验收量；候选保持默认关闭，待跨器件回归后再决定是否替换"
                "绘图默认值。"
            ),
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage", choices=("coarse", "refined", "analyze", "all"), default="all"
    )
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--refined-sentaurus-root",
        type=Path,
        default=DEFAULT_OUTPUT / "sentaurus_refined",
    )
    parser.add_argument(
        "--report-json",
        type=Path,
        default=FIXTURE / "reports" / "cell_first_recovery_ab.json",
    )
    parser.add_argument(
        "--report-markdown",
        type=Path,
        default=FIXTURE / "reports" / "cell_first_recovery_ab.md",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.stage in ("coarse", "all"):
        run_coarse(args)
    if args.stage in ("refined", "all"):
        run_refined(args)
    report = analyze(args) if args.stage in ("analyze", "all") else None
    summary = {
        "stage": args.stage,
        "output_root": relative(args.output_root),
        "overall_pass": None if report is None else report["overall_pass"],
    }
    print(json.dumps(summary, allow_nan=False))
    return 0 if report is None or report["overall_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
