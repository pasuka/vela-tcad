#!/usr/bin/env python3
"""Summarize Templates/LDMOS exact-mesh restart and reclose qualification."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

from templates_ldmos_contracts import validate_document


BASE_STATE_FIELDS = ("psi", "phin", "phip", "electrons_m3", "holes_m3")
PSI_RECLOSE_LIMIT_V = 10.0e-6
RESOLVED_CURRENT_RELATIVE_LIMIT = 1.0e-3
CURRENT_RESOLUTION_FLOOR_A_PER_UM = 1.0e-15


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"{path}: expected at least one data row")
    return rows


def _finite_float(text: str, *, path: Path, column: str) -> float:
    try:
        value = float(text)
    except ValueError as exc:
        raise ValueError(f"{path}: non-numeric {column}={text!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"{path}: non-finite {column}={text!r}")
    return value


def state_max_abs_deltas(first: Path,
                         second: Path,
                         fields: Sequence[str] = BASE_STATE_FIELDS
                         ) -> dict[str, float]:
    first_rows = read_csv(first)
    second_rows = read_csv(second)
    if len(first_rows) != len(second_rows):
        raise ValueError(
            f"state row count differs: {first}={len(first_rows)}, "
            f"{second}={len(second_rows)}")
    result = {field: 0.0 for field in fields}
    for index, (left, right) in enumerate(zip(first_rows, second_rows)):
        if left.get("node_id") != right.get("node_id"):
            raise ValueError(
                f"state node order differs at row {index}: "
                f"{left.get('node_id')!r} != {right.get('node_id')!r}")
        for field in fields:
            if field not in left or field not in right:
                raise ValueError(f"state comparison requires column {field!r}")
            delta = abs(
                _finite_float(left[field], path=first, column=field) -
                _finite_float(right[field], path=second, column=field))
            result[field] = max(result[field], delta)
    return result


def sweep_point(path: Path) -> dict[str, Any]:
    row = read_csv(path)[-1]
    required = (
        "converged", "iterations", "newton_convergence_reason",
        "final_psi_residual_norm", "current_total_A_per_um",
    )
    missing = [field for field in required if field not in row]
    if missing:
        raise ValueError(f"{path}: missing sweep columns {missing}")
    document = {
        "converged": row["converged"] == "1",
        "iterations": int(row["iterations"]),
        "convergence_reason": row["newton_convergence_reason"],
        "final_psi_residual_norm": _finite_float(
            row["final_psi_residual_norm"], path=path,
            column="final_psi_residual_norm"),
        "current_total_A_per_um": _finite_float(
            row["current_total_A_per_um"], path=path,
            column="current_total_A_per_um"),
    }
    return document


def relative_delta(first: float, second: float) -> float:
    return abs(second - first) / max(abs(first), abs(second), 1.0e-300)


def terminal_kcl(path: Path) -> dict[str, float]:
    rows = read_csv(path)
    currents = [
        _finite_float(row["current_total_A_per_um"], path=path,
                      column="current_total_A_per_um")
        for row in rows
    ]
    current_sum = math.fsum(currents)
    scale = max(abs(value) for value in currents)
    document = {
        "terminal_current_sum_A_per_um": current_sum,
        "terminal_current_scale_A_per_um": scale,
        "normalized_imbalance": abs(current_sum) / max(scale, 1.0e-300),
    }
    return document


def _gate(gate_id: str,
          passed: bool,
          metric: float | None,
          limit: float | None,
          unit: str,
          summary: str) -> dict[str, Any]:
    return {
        "id": gate_id,
        "status": "pass" if passed else "fail",
        "metric": metric,
        "limit": limit,
        "unit": unit,
        "summary": summary,
    }


def summarize(root: Path) -> dict[str, Any]:
    files = {
        "frozen_eq_config": root / "frozen_eq_0v.json",
        "sentaurus_eq": root / "sentaurus_eq_0v_state.csv",
        "frozen_eq": root / "frozen_eq_0v_roundtrip.csv",
        "eq_first_config": root / "reclose_eq_0v.json",
        "eq_repeat_config": root / "reclose_eq_0v_repeat.json",
        "eq_first_state": root / "reclose_eq_0v_state.csv",
        "eq_repeat_state": root / "reclose_eq_0v_repeat_state.csv",
        "eq_first_sweep": root / "reclose_eq_0v.csv",
        "eq_repeat_sweep": root / "reclose_eq_0v_repeat.csv",
        "idvg_imported_state": root / "reclose_idvg_vg0_vd0p1_state.csv",
        "idvg_first_config": root / "reclose_idvg_vg0_vd0p1.json",
        "idvg_settled_config": root / "reclose_idvg_vg0_vd0p1_repeat.json",
        "idvg_repeat_config": root / "reclose_idvg_vg0_vd0p1_repeat2.json",
        "idvg_settled_state": root / "reclose_idvg_vg0_vd0p1_repeat_state.csv",
        "idvg_repeat_state": root / "reclose_idvg_vg0_vd0p1_repeat2_state.csv",
        "idvg_settled_sweep": root / "reclose_idvg_vg0_vd0p1_repeat.csv",
        "idvg_repeat_sweep": root / "reclose_idvg_vg0_vd0p1_repeat2.csv",
        "idvd_first_state": root / "reclose_idvd_vg4_vd0p1_state.csv",
        "idvd_first_config": root / "reclose_idvd_vg4_vd0p1.json",
        "idvd_repeat_config": root / "reclose_idvd_vg4_vd0p1_repeat.json",
        "idvd_repeat_state": root / "reclose_idvd_vg4_vd0p1_repeat_state.csv",
        "idvd_first_sweep": root / "reclose_idvd_vg4_vd0p1.csv",
        "idvd_repeat_sweep": root / "reclose_idvd_vg4_vd0p1_repeat.csv",
        "idvd_first_kcl": root / "reclose_idvd_vg4_vd0p1_terminal_balance.csv",
        "idvd_repeat_kcl": root / "reclose_idvd_vg4_vd0p1_repeat_terminal_balance.csv",
    }
    missing = [str(path) for path in files.values() if not path.is_file()]
    if missing:
        raise ValueError(f"restart qualification is missing files: {missing}")

    frozen = state_max_abs_deltas(files["sentaurus_eq"], files["frozen_eq"])
    eq = state_max_abs_deltas(files["eq_first_state"], files["eq_repeat_state"])
    idvg_imported = state_max_abs_deltas(
        files["idvg_imported_state"], files["idvg_settled_state"])
    idvg = state_max_abs_deltas(files["idvg_settled_state"], files["idvg_repeat_state"])
    idvd = state_max_abs_deltas(files["idvd_first_state"], files["idvd_repeat_state"])

    eq_points = [sweep_point(files[name]) for name in ("eq_first_sweep", "eq_repeat_sweep")]
    idvg_points = [
        sweep_point(files[name])
        for name in ("idvg_settled_sweep", "idvg_repeat_sweep")
    ]
    idvd_points = [
        sweep_point(files[name])
        for name in ("idvd_first_sweep", "idvd_repeat_sweep")
    ]
    all_points = eq_points + idvg_points + idvd_points
    idvd_current_delta = relative_delta(
        idvd_points[0]["current_total_A_per_um"],
        idvd_points[1]["current_total_A_per_um"])
    idvg_current_delta = relative_delta(
        idvg_points[0]["current_total_A_per_um"],
        idvg_points[1]["current_total_A_per_um"])
    idvg_current_resolved = max(
        abs(idvg_points[0]["current_total_A_per_um"]),
        abs(idvg_points[1]["current_total_A_per_um"]),
    ) >= CURRENT_RESOLUTION_FLOOR_A_PER_UM
    kcl_first = terminal_kcl(files["idvd_first_kcl"])
    kcl_repeat = terminal_kcl(files["idvd_repeat_kcl"])

    gates = [
        _gate(
            "binary64_state_roundtrip", all(value == 0.0 for value in frozen.values()),
            max(frozen.values()), 0.0, "stored_field_unit",
            "The 17-digit CSV read/write replay preserves every classical persistent field."),
        _gate(
            "frozen_state_no_change", all(value == 0.0 for value in frozen.values()),
            max(frozen.values()), 0.0, "stored_field_unit",
            "frozen_state returns the supplied exact-mesh state without modification."),
        _gate(
            "all_reclose_runs_converged", all(point["converged"] for point in all_points),
            float(sum(not point["converged"] for point in all_points)), 0.0, "count",
            "All equilibrium, low-drain IdVg and IdVd reclose runs converged."),
        _gate(
            "equilibrium_same_bias_psi", eq["psi"] <= PSI_RECLOSE_LIMIT_V,
            eq["psi"], PSI_RECLOSE_LIMIT_V, "V",
            "The equilibrium Vela fixed point recloses on the same branch."),
        _gate(
            "idvg_low_drain_same_bias_psi", idvg["psi"] <= PSI_RECLOSE_LIMIT_V,
            idvg["psi"], PSI_RECLOSE_LIMIT_V, "V",
            "The settled low-drain IdVg Vela fixed point recloses on the same branch."),
        _gate(
            "idvd_prebias_same_bias_psi", idvd["psi"] <= PSI_RECLOSE_LIMIT_V,
            idvd["psi"], PSI_RECLOSE_LIMIT_V, "V",
            "The resolved-current IdVd prebias state recloses on the same branch."),
        _gate(
            "idvd_prebias_current", idvd_current_delta <= RESOLVED_CURRENT_RELATIVE_LIMIT,
            idvd_current_delta, RESOLVED_CURRENT_RELATIVE_LIMIT, "fraction",
            "The resolved drain current is stable across the same-bias reclose."),
        _gate(
            "idvd_prebias_kcl_no_regression",
            kcl_repeat["normalized_imbalance"] <= kcl_first["normalized_imbalance"],
            kcl_repeat["normalized_imbalance"], kcl_first["normalized_imbalance"],
            "fraction",
            "Normalized terminal-current imbalance does not regress."),
    ]
    status = "pass" if all(gate["status"] == "pass" for gate in gates) else "fail"
    document = {
        "schema": "vela.templates_ldmos.restart_qualification.v1",
        "benchmark": "sentaurus_t2022_03_sp2_templates_ldmos",
        "generated_at": utc_now(),
        "status": status,
        "thresholds": {
            "psi_reclose_limit_V": PSI_RECLOSE_LIMIT_V,
            "resolved_current_relative_limit": RESOLVED_CURRENT_RELATIVE_LIMIT,
            "current_resolution_floor_A_per_um": CURRENT_RESOLUTION_FLOOR_A_PER_UM,
        },
        "gates": gates,
        "cases": {
            "frozen_equilibrium": {"max_abs_delta": frozen},
            "equilibrium_reclose": {"max_abs_delta": eq, "sweep_points": eq_points},
            "idvg_low_drain_reclose": {
                "max_abs_imported_to_settled_delta": idvg_imported,
                "max_abs_settled_reclose_delta": idvg,
                "sweep_points": idvg_points,
                "current_relative_delta": idvg_current_delta,
                "current_resolved": idvg_current_resolved,
                "classification": (
                    "resolved_current" if idvg_current_resolved
                    else "below_absolute_current_resolution_floor_non_gating"),
            },
            "idvd_prebias_reclose": {
                "max_abs_delta": idvd,
                "sweep_points": idvd_points,
                "current_relative_delta": idvd_current_delta,
                "kcl_before": kcl_first,
                "kcl_after": kcl_repeat,
            },
        },
        "interpretation": [
            "Sentaurus-to-Vela state closure is diagnostic only; it includes boundary, physics and discretization differences.",
            "The low-drain IdVg gate compares settled Vela fixed points. Its terminal current is below the declared absolute resolution floor, so relative current and KCL are report-only.",
            "The stage-1.5 gate qualifies serialization, frozen replay and same-bias closure; it does not approve the provisional metal-gate flat-band mapping or stage-2 physics.",
        ],
        "evidence": [
            {"role": role, "path": path.name, "sha256": sha256_file(path)}
            for role, path in sorted(files.items())
        ],
    }
    validate_document(document)
    return document


def render_markdown(document: dict[str, Any]) -> str:
    lines = [
        "# Templates/LDMOS stage 1.5 restart qualification",
        "",
        f"- Status: `{document['status']}`",
        f"- Generated: `{document['generated_at']}`",
        "- Scope: exact-mesh serialization, frozen replay and same-bias reclose",
        "",
        "## Gates",
        "",
        "| Gate | Status | Metric | Limit | Unit |",
        "| --- | --- | ---: | ---: | --- |",
    ]
    for gate in document["gates"]:
        lines.append(
            f"| `{gate['id']}` | `{gate['status']}` | {gate['metric']:.12g} | "
            f"{gate['limit']:.12g} | {gate['unit']} |")
    idvg = document["cases"]["idvg_low_drain_reclose"]
    idvd = document["cases"]["idvd_prebias_reclose"]
    lines.extend([
        "",
        "## Key diagnostics",
        "",
        f"- IdVg imported-to-settled psi change: `{idvg['max_abs_imported_to_settled_delta']['psi']:.12g} V` (non-gating).",
        f"- IdVg settled reclose psi change: `{idvg['max_abs_settled_reclose_delta']['psi']:.12g} V`.",
        f"- IdVg current classification: `{idvg['classification']}`.",
        f"- IdVd reclose current relative delta: `{idvd['current_relative_delta']:.12g}`.",
        f"- IdVd normalized KCL imbalance: `{idvd['kcl_before']['normalized_imbalance']:.12g}` -> `{idvd['kcl_after']['normalized_imbalance']:.12g}`.",
        "",
        "## Interpretation",
        "",
    ])
    lines.extend(f"- {item}" for item in document["interpretation"])
    return "\n".join(lines) + "\n"


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qualification-dir", type=Path, required=True)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--output-md", type=Path)
    args = parser.parse_args(argv)
    document = summarize(args.qualification_dir)
    output_json = args.output_json or args.qualification_dir / "qualification_summary.json"
    output_md = args.output_md or args.qualification_dir / "qualification_summary.md"
    output_json.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    output_md.write_text(render_markdown(document), encoding="utf-8")
    print(json.dumps({"status": document["status"], "output": str(output_json)}))
    return 0 if document["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
