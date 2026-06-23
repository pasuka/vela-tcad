#!/usr/bin/env python3
"""Compare a restart initial state against its converged state on BV support nodes."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from pathlib import Path
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import diagnose_pn2d_bv_sg_flux_forms as fluxforms
import diagnose_pn2d_bv_active_edge_flux_factors as fluxfac


STATE_FIELDS = ["psi", "phin", "phip", "n", "p"]
POTENTIAL_FIELDS = ["psi", "phin", "phip"]
DENSITY_FIELDS = ["n", "p"]


def optional_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_state_csv(path: Path, node_count: int | None = None) -> dict[str, list[float]]:
    rows = read_rows(path)
    if node_count is None:
        node_count = max(int(row["node_id"]) for row in rows) + 1
    values = {field: [0.0] * node_count for field in STATE_FIELDS}
    for row in rows:
        node_id = int(row["node_id"])
        values["psi"][node_id] = float(row["psi"])
        values["phin"][node_id] = float(row["phin"])
        values["phip"][node_id] = float(row["phip"])
        values["n"][node_id] = float(row["electrons_m3"])
        values["p"][node_id] = float(row["holes_m3"])
    return values


def load_baseline(args: argparse.Namespace, node_count: int) -> dict[str, list[float]]:
    if args.baseline_state_csv is not None:
        return load_state_csv(args.baseline_state_csv, node_count)
    if args.baseline_vtk_root is None or args.mesh is None or args.bias is None:
        raise SystemExit("baseline requires --baseline-state-csv or --mesh + --baseline-vtk-root + --bias")
    return fluxforms.load_vela_state(fluxfac.discover_vtk(args.baseline_vtk_root, args.bias), node_count)


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def ratio(num: float, den: float) -> float | None:
    if den == 0.0 or not math.isfinite(num) or not math.isfinite(den):
        return None
    return num / den


def support_nodes(path: Path) -> dict[str, list[int]]:
    groups: dict[str, list[int]] = {}
    for row in read_rows(path):
        support_class = row.get("support_class", "")
        if support_class in ("", "inactive"):
            continue
        groups.setdefault(support_class, []).append(int(row["node_id"]))
    return groups


def finite(items: list[float | None]) -> list[float]:
    return [item for item in items if item is not None and math.isfinite(item)]


def summarize_class(
    node_ids: list[int],
    initial: dict[str, list[float]],
    final: dict[str, list[float]],
    baseline: dict[str, list[float]],
) -> dict[str, Any]:
    out: dict[str, Any] = {"count": len(node_ids)}
    for field in POTENTIAL_FIELDS:
        initial_shift = [initial[field][node_id] - baseline[field][node_id] for node_id in node_ids]
        final_shift = [final[field][node_id] - baseline[field][node_id] for node_id in node_ids]
        relaxation = [final[field][node_id] - initial[field][node_id] for node_id in node_ids]
        retained = [ratio(abs(fin), abs(init)) for fin, init in zip(final_shift, initial_shift)]
        out[f"initial_minus_baseline_{field}_median_V"] = median(initial_shift)
        out[f"final_minus_baseline_{field}_median_V"] = median(final_shift)
        out[f"final_minus_initial_{field}_median_V"] = median(relaxation)
        out[f"abs_shift_retained_{field}_median"] = median(finite(retained))
    for field in DENSITY_FIELDS:
        out[f"initial_over_baseline_{field}_median"] = median(finite([
            ratio(initial[field][node_id], baseline[field][node_id]) for node_id in node_ids
        ]))
        out[f"final_over_baseline_{field}_median"] = median(finite([
            ratio(final[field][node_id], baseline[field][node_id]) for node_id in node_ids
        ]))
        out[f"final_over_initial_{field}_median"] = median(finite([
            ratio(final[field][node_id], initial[field][node_id]) for node_id in node_ids
        ]))
    return out


def node_rows(
    groups: dict[str, list[int]],
    initial: dict[str, list[float]],
    final: dict[str, list[float]],
    baseline: dict[str, list[float]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for support_class, node_ids in groups.items():
        for node_id in node_ids:
            row: dict[str, Any] = {"support_class": support_class, "node_id": node_id}
            for field in POTENTIAL_FIELDS:
                row[f"initial_minus_baseline_{field}_V"] = initial[field][node_id] - baseline[field][node_id]
                row[f"final_minus_baseline_{field}_V"] = final[field][node_id] - baseline[field][node_id]
                row[f"final_minus_initial_{field}_V"] = final[field][node_id] - initial[field][node_id]
            for field in DENSITY_FIELDS:
                row[f"initial_over_baseline_{field}"] = ratio(initial[field][node_id], baseline[field][node_id])
                row[f"final_over_baseline_{field}"] = ratio(final[field][node_id], baseline[field][node_id])
                row[f"final_over_initial_{field}"] = ratio(final[field][node_id], initial[field][node_id])
            rows.append(row)
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def clean_json(value: Any) -> Any:
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    return value


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--initial-state-csv", type=Path, required=True)
    parser.add_argument("--final-state-csv", type=Path, required=True)
    parser.add_argument("--baseline-state-csv", type=Path)
    parser.add_argument("--mesh", type=Path)
    parser.add_argument("--baseline-vtk-root", type=Path)
    parser.add_argument("--bias", type=float)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    if args.mesh is not None:
        nodes, _, _ = fluxforms.sgdiag.read_mesh(args.mesh)
        node_count = len(nodes)
    else:
        node_count = None
    initial = load_state_csv(args.initial_state_csv, node_count)
    node_count = len(initial["psi"])
    final = load_state_csv(args.final_state_csv, node_count)
    baseline = load_baseline(args, node_count)
    groups = support_nodes(args.support_csv)

    summary = {
        "support_classes": {
            support_class: summarize_class(node_ids, initial, final, baseline)
            for support_class, node_ids in groups.items()
        },
        "support_class_count": len(groups),
    }
    rows = node_rows(groups, initial, final, baseline)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "restart_state_relaxation_nodes.csv", rows)
    (args.out_dir / "restart_state_relaxation_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
