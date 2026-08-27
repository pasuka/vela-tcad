#!/usr/bin/env python3
"""Summarize Templates/LDMOS phase-2 evidence and the phase-3 stop gate."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any, Sequence


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def state(path: Path) -> dict[int, dict[str, float]]:
    return {
        int(row["node_id"]): {
            key: float(value) for key, value in row.items() if key != "node_id"
        }
        for row in rows(path)
    }


def optional_state(path: Path) -> dict[int, dict[str, float]] | None:
    return state(path) if path.is_file() else None


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        raise ValueError("cannot summarize an empty sample")
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def distribution(values: list[float]) -> dict[str, float | int]:
    return {
        "count": len(values),
        "median": statistics.median(values),
        "p95": percentile(values, 0.95),
        "max": max(values),
    }


def field_metrics(reference: dict[int, dict[str, float]],
                  candidate: dict[int, dict[str, float]]) -> dict[str, Any]:
    common = sorted(set(reference) & set(candidate))
    semiconductor = [
        node for node in common
        if reference[node].get("electrons_m3", 0.0) > 0.0 or
        reference[node].get("holes_m3", 0.0) > 0.0
    ]
    result: dict[str, Any] = {
        "semiconductor_node_count": len(semiconductor),
        "psi_abs_error_V": distribution([
            abs(candidate[node]["psi"] - reference[node]["psi"])
            for node in semiconductor
        ]),
    }
    for column, label in (("electrons_m3", "electron"), ("holes_m3", "hole")):
        peak = max(reference[node][column] for node in semiconductor)
        support = [node for node in semiconductor if reference[node][column] > peak * 1.0e-12]
        result[f"{label}_abs_error_dex"] = distribution([
            abs(math.log10(max(candidate[node][column], 1.0)) -
                math.log10(max(reference[node][column], 1.0)))
            for node in support
        ])
        result[f"{label}_support_floor_m3"] = peak * 1.0e-12
    return result


def sweep_status(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"points": 0, "all_converged": False, "last": None,
                "missing": str(path)}
    data = rows(path)
    return {
        "points": len(data),
        "all_converged": bool(data) and all(int(row["converged"]) == 1 for row in data),
        "last": ({
            "bias_V": float(data[-1]["bias_V"]),
            "converged": int(data[-1]["converged"]) == 1,
            "iterations": int(data[-1]["iterations"]),
            "reason": data[-1]["newton_convergence_reason"],
            "failure_reason": data[-1]["failure_reason"],
            "psi_residual": float(data[-1]["final_psi_residual_norm"]),
        } if data else None),
    }


def formula_density_metrics(reference: dict[int, dict[str, float]],
                            sg_path: Path) -> dict[str, Any]:
    accumulated: dict[int, dict[str, list[float]]] = {}
    for row in rows(sg_path):
        for endpoint in (0, 1):
            node = int(row[f"node{endpoint}"])
            entry = accumulated.setdefault(node, {"n": [], "p": []})
            entry["n"].append(float(row[f"electron_density{endpoint}_m3"]))
            entry["p"].append(float(row[f"hole_density{endpoint}_m3"]))
    result: dict[str, Any] = {}
    for source, key, label in (
        ("electrons_m3", "n", "electron"), ("holes_m3", "p", "hole")):
        common = [node for node in accumulated if node in reference and reference[node][source] > 0.0]
        peak = max(reference[node][source] for node in common)
        support = [node for node in common if reference[node][source] > peak * 1.0e-12]
        result[f"{label}_formula_abs_error_dex"] = distribution([
            abs(math.log10(max(statistics.median(accumulated[node][key]), 1.0)) -
                math.log10(max(reference[node][source], 1.0)))
            for node in support
        ])
    return result


def gate(status: bool, identifier: str, summary: str, metrics: dict[str, Any] | None = None) -> dict[str, Any]:
    return {
        "id": identifier,
        "status": "pass" if status else "fail",
        "summary": summary,
        "metrics": metrics or {},
    }


def summarize(stage1_dir: Path, output_dir: Path) -> dict[str, Any]:
    qualification = stage1_dir / "qualification"
    manifest = json.loads((output_dir / "deck_manifest.json").read_text(encoding="utf-8"))
    oracle = manifest["classical_equilibrium_oracle"]
    reference = state(Path(oracle["path"]))
    mapped = optional_state(output_dir / "g_contact_polysi_eq_repeat_state.csv")
    provisional = optional_state(output_dir / "g_contact_provisional_eq_state.csv")
    mapped_metrics = field_metrics(reference, mapped) if mapped else {
        "status": "unavailable_no_converged_repeat_state"}
    provisional_metrics = field_metrics(reference, provisional) if provisional else {
        "status": "unavailable_no_converged_state"}
    mapping = manifest["gate_mapping"]
    mapped_status = sweep_status(output_dir / "g_contact_polysi_eq.csv")
    repeat_status = sweep_status(output_dir / "g_contact_polysi_eq_repeat.csv")
    drain_ramp = sweep_status(output_dir / "g3_drain_prebias.csv")
    idvg_seed = sweep_status(output_dir / "g3_idvg_seed.csv")
    idvg_reclose = sweep_status(output_dir / "g3_idvg_seed_repeat.csv")
    sg_path = output_dir / "sentaurus_eq_sg_edge_flux.csv"
    formula = formula_density_metrics(reference, sg_path) if sg_path.is_file() else {
        "status": "unavailable_missing_fixed_state_probe"}

    field_comparable = False
    if mapped:
        psi = mapped_metrics["psi_abs_error_V"]
        electron = mapped_metrics["electron_abs_error_dex"]
        hole = mapped_metrics["hole_abs_error_dex"]
        field_comparable = (
            psi["median"] <= 2.0e-3 and psi["p95"] <= 10.0e-3 and
            electron["median"] <= 0.03 and electron["p95"] <= 0.10 and
            hole["median"] <= 0.03 and hole["p95"] <= 0.10
        )
    gates = [
        gate(mapping["spread_V"] <= 1.0e-12, "polysi_gate_mapping",
             "PolySi(N) offset is derived from a constant sealed gate boundary.", mapping),
        gate(mapped_status["all_converged"] and repeat_status["all_converged"],
             "equilibrium_same_bias_reclose",
             "Mapped equilibrium and its explicit same-bias repeat must converge.",
             {"first": mapped_status, "repeat": repeat_status}),
        gate(field_comparable and oracle["qualified_for_spatial_scoring"],
             "equilibrium_spatial_comparable",
             "Mapped equilibrium is checked against the sealed state on common semiconductor nodes.",
             mapped_metrics),
        gate(oracle["qualified_for_spatial_scoring"], "sentaurus_g4_equilibrium_oracle",
             "The G4/equivalent classical Sentaurus equilibrium ablation must be normalized before spatial scoring.",
             oracle),
        gate(drain_ramp["all_converged"], "g3_equilibrium_to_low_drain_ramp",
             "The classical mapped equilibrium must continue to Vd=0.1 V.", drain_ramp),
        gate(idvg_reclose["all_converged"], "g3_low_drain_same_bias_reclose",
             "The alternate Sentaurus-state G3 seed must reclose before an Id-Vg sweep.",
             {"seed": idvg_seed, "repeat": idvg_reclose}),
    ]
    hard_failures = [item["id"] for item in gates if item["status"] == "fail"]
    summary = {
        "schema": "vela.templates_ldmos.phase23_summary.v1",
        "benchmark": "sentaurus_t2022_03_sp2_templates_ldmos",
        "status": "fail" if hard_failures else "pass",
        "highest_level": "L1",
        "stop_rule_triggered": True,
        "gates": gates,
        "comparisons": {
            "mapped_polysi_vs_selected_classical_oracle": mapped_metrics,
            "provisional_zero_flatband_vs_selected_classical_oracle": provisional_metrics,
            "fixed_state_formula_replay": formula,
        },
        "hard_failures": hard_failures,
        "interpretation": [
            "The PolySi boundary mapping is closed independently of the transport solve.",
            "The original G0 state contains a different physics layer and cannot substitute for a G4 classical equilibrium oracle.",
            "The correct mapped boundary exposes a low-drain continuation/reclose failure; stage 3 is stopped before curve scoring.",
            "No G0-to-G3 amplitude score is reported because that would mix QP and IALMob deltas.",
        ],
    }
    (output_dir / "phase23_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# Templates/LDMOS phase 2/3 execution summary", "",
        f"Status: **{summary['status']}**; highest qualified level: **L1**.", "",
        "| Gate | Status | Summary |", "| --- | --- | --- |",
    ]
    lines.extend(
        f"| `{item['id']}` | `{item['status']}` | {item['summary']} |"
        for item in gates
    )
    lines.extend(["", "## Stop decision", "", *(
        f"- `{identifier}`" for identifier in hard_failures), "",
        "The Id-Vg curve was not scored because the G3 low-drain fixed point did not pass same-bias reclose.",
    ])
    (output_dir / "phase23_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage1-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    summary = summarize(args.stage1_dir.resolve(), args.output_dir.resolve())
    print(json.dumps(summary, indent=2))
    return 0 if summary["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
