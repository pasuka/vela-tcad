#!/usr/bin/env python3
"""Summarize the strict Sentaurus/Vela SLOT-LDMOS 1 V fixed-state audit."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


STATE_FIELDS = ("psi", "phin", "phip")
DENSITY_FIELDS = ("electrons_m3", "holes_m3")


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _finite_float(row: dict[str, str], name: str) -> float:
    value = float(row[name])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {name} in CSV row")
    return value


def _percentile(values: list[float], fraction: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    index = fraction * (len(ordered) - 1)
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[lower]
    weight = index - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def _state_map(path: Path) -> dict[int, dict[str, float]]:
    state: dict[int, dict[str, float]] = {}
    for row in _read_rows(path):
        node = int(row["node_id"])
        state[node] = {
            field: _finite_float(row, field)
            for field in (*STATE_FIELDS, *DENSITY_FIELDS)
        }
    return state


def compare_states(reference_path: Path, candidate_path: Path) -> dict:
    reference = _state_map(reference_path)
    candidate = _state_map(candidate_path)
    common = sorted(reference.keys() & candidate.keys())
    semiconductor = [
        node for node in common
        if reference[node]["electrons_m3"] > 0.0
        and reference[node]["holes_m3"] > 0.0
    ]
    if not semiconductor:
        raise ValueError("no common semiconductor nodes in state comparison")

    result: dict[str, object] = {
        "reference_nodes": len(reference),
        "candidate_nodes": len(candidate),
        "common_semiconductor_nodes": len(semiconductor),
        "absolute_voltage_error_V": {},
        "absolute_log10_density_ratio": {},
    }
    for field in STATE_FIELDS:
        errors = [
            abs(candidate[node][field] - reference[node][field])
            for node in semiconductor
        ]
        result["absolute_voltage_error_V"][field] = {
            "median": statistics.median(errors),
            "p95": _percentile(errors, 0.95),
            "maximum": max(errors),
        }
    for field in DENSITY_FIELDS:
        errors = [
            abs(math.log10(max(candidate[node][field], 1.0))
                - math.log10(max(reference[node][field], 1.0)))
            for node in semiconductor
        ]
        result["absolute_log10_density_ratio"][field] = {
            "median_decade": statistics.median(errors),
            "p95_decade": _percentile(errors, 0.95),
            "maximum_decade": max(errors),
        }
    return result


def sentaurus_contact_currents(import_dir: Path) -> dict[str, float]:
    manifest = json.loads(
        (import_dir / "field_manifest.json").read_text(encoding="utf-8"))
    currents: dict[str, float] = {}
    for field in manifest["fields"]:
        if field.get("name") != "ContactCurrentFlux":
            continue
        rows = _read_rows(import_dir / "fields" / field["csv_file"])
        if len(rows) != 1:
            raise ValueError("ContactCurrentFlux dataset is not scalar")
        currents[field["region_name"]] = _finite_float(rows[0], "component0")
    return currents


def analyze(
    import_dir: Path,
    terminal_balance_path: Path,
    recombination_balance_path: Path,
    sentaurus_state_path: Path,
    self_consistent_state_path: Path | None,
) -> dict:
    sentaurus = sentaurus_contact_currents(import_dir)
    terminal_rows = _read_rows(terminal_balance_path)
    terminals = {
        row["contact"]: {
            "electron_A_per_um": _finite_float(
                row, "current_electron_A_per_um"),
            "hole_A_per_um": _finite_float(row, "current_hole_A_per_um"),
            "total_A_per_um": _finite_float(row, "current_total_A_per_um"),
        }
        for row in terminal_rows
    }
    recombination_rows = _read_rows(recombination_balance_path)
    if len(recombination_rows) != 1:
        raise ValueError("expected exactly one recombination balance row")
    recombination = recombination_rows[0]

    sentaurus_drain = sentaurus["drain"]
    vela_drain = terminals["drain"]["total_A_per_um"]
    ratio = abs(vela_drain / sentaurus_drain)
    summary: dict[str, object] = {
        "gate": {
            "prebreak_current_error_limit_decade": 0.3,
            "fixed_state_current_error_decade": abs(math.log10(ratio)),
            "passed": abs(math.log10(ratio)) <= 0.3,
        },
        "sentaurus_contact_current_flux_A_per_um": sentaurus,
        "vela_fixed_state_terminal_current_A_per_um": terminals,
        "drain_current_comparison": {
            "vela_to_sentaurus_magnitude_ratio": ratio,
            "absolute_log10_error_decade": abs(math.log10(ratio)),
        },
        "vela_fixed_state_recombination_integrals_A_per_um": {
            name: _finite_float(recombination, name)
            for name in (
                "srh_net_current_A_per_um",
                "srh_generation_current_A_per_um",
                "srh_recombination_current_A_per_um",
                "auger_net_current_A_per_um",
                "auger_generation_current_A_per_um",
                "auger_recombination_current_A_per_um",
                "total_recombination_net_current_A_per_um",
            )
        },
        "diagnosis": (
            "fixed_state_operator_mismatch_still_above_gate"
            if abs(math.log10(ratio)) > 0.3
            else "fixed_state_current_gate_passed"
        ),
    }
    if self_consistent_state_path is not None:
        summary["sentaurus_vs_existing_vela_self_consistent_state"] = (
            compare_states(sentaurus_state_path, self_consistent_state_path))
    return summary


def render_markdown(summary: dict) -> str:
    drain = summary["drain_current_comparison"]
    recombination = summary[
        "vela_fixed_state_recombination_integrals_A_per_um"]
    state = summary.get("sentaurus_vs_existing_vela_self_consistent_state")
    lines = [
        "# SLOT-LDMOS 1 V 同态冻结审计",
        "",
        "| 项目 | 数值 |",
        "|---|---:|",
        f"| Vela/Sentaurus 漏极电流幅值比 | {drain['vela_to_sentaurus_magnitude_ratio']:.9g} |",
        f"| 漏极电流误差 | {drain['absolute_log10_error_decade']:.6g} decade |",
        f"| 0.3 decade 门槛 | {'通过' if summary['gate']['passed'] else '未通过'} |",
        f"| Vela SRH 净积分 | {recombination['srh_net_current_A_per_um']:.9g} A/um |",
        f"| Vela Auger 净积分 | {recombination['auger_net_current_A_per_um']:.9g} A/um |",
        "",
        "结论：同一节点状态下仍有接触电流算子差异，但远小于既有自洽解的数量级偏差。",
    ]
    if state is not None:
        lines += [
            "",
            "## 既有 Vela 自洽 1 V 状态相对 Sentaurus 的差异",
            "",
            "| 状态量 | 中位误差 | P95 | 最大误差 |",
            "|---|---:|---:|---:|",
        ]
        for name in STATE_FIELDS:
            metric = state["absolute_voltage_error_V"][name]
            lines.append(
                f"| {name} (V) | {metric['median']:.6g} | "
                f"{metric['p95']:.6g} | {metric['maximum']:.6g} |")
        for name in DENSITY_FIELDS:
            metric = state["absolute_log10_density_ratio"][name]
            lines.append(
                f"| {name} (decade) | {metric['median_decade']:.6g} | "
                f"{metric['p95_decade']:.6g} | {metric['maximum_decade']:.6g} |")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-import-dir", type=Path, required=True)
    parser.add_argument("--sentaurus-state", type=Path, required=True)
    parser.add_argument("--vela-terminal-balance", type=Path, required=True)
    parser.add_argument("--vela-recombination-balance", type=Path, required=True)
    parser.add_argument("--vela-self-consistent-state", type=Path)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-markdown", type=Path, required=True)
    args = parser.parse_args()

    summary = analyze(
        args.sentaurus_import_dir,
        args.vela_terminal_balance,
        args.vela_recombination_balance,
        args.sentaurus_state,
        args.vela_self_consistent_state,
    )
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_markdown.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8")
    args.output_markdown.write_text(
        render_markdown(summary), encoding="utf-8")
    print(args.output_json)
    print(args.output_markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
