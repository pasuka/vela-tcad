#!/usr/bin/env python3
"""Run PN2D BV no-impact variants and compare high-bias carrier branches."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_branch_profiles as branch_profiles


REPO = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Variant:
    tag: str
    description: str
    solver_updates: dict[str, Any]


VARIANTS = [
    Variant("baseline", "no-impact baseline: SRH + Masetti high-field QF mobility", {}),
    Variant(
        "srh_tau_equal_1e_m5",
        "SRH with equal electron/hole lifetimes of 1e-5 s",
        {"taun": 1.0e-5, "taup": 1.0e-5},
    ),
    Variant(
        "srh_tau_equal_1e_m6",
        "SRH with equal electron/hole lifetimes of 1e-6 s",
        {"taun": 1.0e-6, "taup": 1.0e-6},
    ),
    Variant(
        "srh_tau_equal_1e_m7",
        "SRH with equal electron/hole lifetimes of 1e-7 s",
        {"taun": 1.0e-7, "taup": 1.0e-7},
    ),
    Variant(
        "srh_tau_equal_1e_m8",
        "SRH with equal electron/hole lifetimes of 1e-8 s",
        {"taun": 1.0e-8, "taup": 1.0e-8},
    ),
    Variant(
        "srh_tau_asym_1e_m5_3e_m6",
        "explicit baseline SRH lifetimes: taun=1e-5 s, taup=3e-6 s",
        {"taun": 1.0e-5, "taup": 3.0e-6},
    ),
    Variant("no_srh", "disable SRH recombination", {"recombination": ["none"]}),
    Variant(
        "low_field_masetti",
        "disable high-field mobility limiting, keep Masetti doping dependence",
        {"mobility": {"model": "masetti"}},
    ),
    Variant(
        "electric_field_drive",
        "keep high-field mobility but drive it with electric field",
        {"mobility": {"high_field_driving_force": "electric_field"}},
    ),
    Variant(
        "contact_relax_n_0p10",
        "enable n-contact minority-electron relaxation above 0.10 V",
        {
            "contact_boundary_minority_electron_relaxation": True,
            "contact_boundary_minority_electron_relaxation_contact_side": "n_contact_only",
            "contact_boundary_minority_electron_relaxation_bias_threshold_V": 0.10,
        },
    ),
    Variant(
        "no_srh_low_field_masetti",
        "disable SRH and high-field mobility limiting together",
        {"recombination": ["none"], "mobility": {"model": "masetti"}},
    ),
]


SUMMARY_FIELDS = [
    "variant",
    "description",
    "bias_V",
    "band",
    "node_count",
    "candidate_minus_sentaurus_psi_median_V",
    "candidate_minus_baseline_psi_median_V",
    "candidate_minus_sentaurus_psi_minus_phin_median_V",
    "candidate_minus_baseline_psi_minus_phin_median_V",
    "candidate_minus_sentaurus_phip_minus_psi_median_V",
    "candidate_minus_baseline_phip_minus_psi_median_V",
    "candidate_log10_electron_over_sentaurus_median",
    "candidate_log10_electron_over_baseline_median",
    "candidate_log10_hole_over_sentaurus_median",
    "candidate_log10_hole_over_baseline_median",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", default="-10,-13.2")
    parser.add_argument(
        "--bands",
        default=(
            "pre_junction_p:0.7:0.9,junction:0.95:1.05,"
            "post_junction_n:1.1:1.3"
        ),
    )
    parser.add_argument(
        "--variants",
        default=",".join(variant.tag for variant in VARIANTS),
        help="Comma-separated variant tags to run.",
    )
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--skip-run", action="store_true")
    return parser.parse_args()


def parse_float_list(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def parse_tags(raw: str) -> list[str]:
    return [item.strip() for item in raw.split(",") if item.strip()]


def selected_variants(tags: list[str]) -> list[Variant]:
    by_tag = {variant.tag: variant for variant in VARIANTS}
    missing = [tag for tag in tags if tag not in by_tag]
    if missing:
        raise SystemExit(f"unknown variants: {', '.join(missing)}")
    return [by_tag[tag] for tag in tags]


def deep_update(base: dict[str, Any], updates: dict[str, Any]) -> dict[str, Any]:
    result = json.loads(json.dumps(base))
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = deep_update(result[key], value)
        else:
            result[key] = value
    return result


def prepare_config(base: dict[str, Any],
                   base_dir: Path,
                   case_dir: Path,
                   variant: Variant,
                   biases: list[float]) -> Path:
    cfg = json.loads(json.dumps(base))
    cfg["mesh_file"] = str((base_dir / cfg["mesh_file"]).resolve())
    if "node_doping_file" in cfg:
        cfg["node_doping_file"] = str((base_dir / cfg["node_doping_file"]).resolve())
    cfg["output_csv"] = f"{variant.tag}.csv"
    solver = cfg.setdefault("solver", {})
    solver["impact_ionization"] = {"model": "none"}
    solver.update(deep_update(solver, variant.solver_updates))
    sweep = cfg.setdefault("sweep", {})
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = variant.tag
    run_biases = list(biases)
    if run_biases and abs(run_biases[0]) > 1.0e-15:
        run_biases.insert(0, 0.0)
    sweep["bias_points"] = run_biases
    if run_biases:
        sweep["start"] = run_biases[0]
        sweep["stop"] = run_biases[-1]
    case_dir.mkdir(parents=True, exist_ok=True)
    config_path = case_dir / f"{variant.tag}.json"
    config_path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
    return config_path.resolve()


def run_variant(runner: Path, config: Path, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(runner), "--config", str(config)],
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def read_nodes(mesh: Path) -> list[dict[str, float]]:
    return branch_profiles.read_nodes(mesh)


def load_variant_state(vtk_root: Path, bias: float, node_count: int) -> dict[str, list[float]]:
    vtks = branch_profiles.discover_vela_vtks(vtk_root)
    vtk = vtks.get(branch_profiles.bias_key(bias))
    if vtk is None:
        raise FileNotFoundError(f"missing VTK for {bias:g} V under {vtk_root}")
    return branch_profiles.load_vela_state(vtk, node_count)


def median(values: list[float]) -> float | None:
    finite = [value for value in values if math.isfinite(value)]
    return statistics.median(finite) if finite else None


def take(values: list[float], ids: list[int]) -> list[float]:
    return [values[idx] for idx in ids]


def diff_median(field: str,
                ids: list[int],
                candidate: dict[str, list[float]],
                reference: dict[str, list[float]]) -> float | None:
    return branch_profiles.diff_median(field, ids, candidate, reference)


def comparison_row(variant: Variant,
                   bias: float,
                   band: str,
                   ids: list[int],
                   candidate: dict[str, list[float]],
                   sentaurus: dict[str, list[float]],
                   baseline: dict[str, list[float]]) -> dict[str, Any]:
    row: dict[str, Any] = {
        "variant": variant.tag,
        "description": variant.description,
        "bias_V": bias,
        "band": band,
        "node_count": len(ids),
    }
    for field in ["psi", "psi_minus_phin", "phip_minus_psi"]:
        row[f"candidate_minus_sentaurus_{field}_median_V"] = diff_median(
            field, ids, candidate, sentaurus)
        row[f"candidate_minus_baseline_{field}_median_V"] = diff_median(
            field, ids, candidate, baseline)
    for carrier in ["electron", "hole"]:
        row[f"candidate_log10_{carrier}_over_sentaurus_median"] = (
            branch_profiles.median_log_ratio(
                take(candidate[carrier], ids), take(sentaurus[carrier], ids)))
        row[f"candidate_log10_{carrier}_over_baseline_median"] = (
            branch_profiles.median_log_ratio(
                take(candidate[carrier], ids), take(baseline[carrier], ids)))
    return row


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main() -> int:
    args = parse_args()
    base_config = args.base_config.resolve()
    base_dir = base_config.parent
    base = json.loads(base_config.read_text())
    biases = parse_float_list(args.biases)
    bands = branch_profiles.parse_bands(args.bands)
    variants = selected_variants(parse_tags(args.variants))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    cases_dir = args.out_dir / "cases"
    config_rows = []
    for variant in variants:
        case_dir = cases_dir / variant.tag
        config_path = prepare_config(base, base_dir, case_dir, variant, biases)
        config_rows.append({
            "variant": variant.tag,
            "description": variant.description,
            "config": str(config_path),
            "case_dir": str(case_dir.resolve()),
        })

    if args.prepare_only:
        (args.out_dir / "noimpact_variant_scan_summary.json").write_text(
            json.dumps({"prepared": config_rows}, indent=2) + "\n",
            encoding="utf-8",
        )
        return 0

    run_rows = []
    if not args.skip_run:
        for item in config_rows:
            case_dir = Path(item["case_dir"])
            completed = run_variant(args.runner.resolve(), Path(item["config"]), case_dir)
            (case_dir / "runner.log").write_text(completed.stdout, encoding="utf-8")
            run_rows.append({
                **item,
                "returncode": completed.returncode,
                "status": "ok" if completed.returncode == 0 else "failed",
            })
            if completed.returncode != 0:
                raise RuntimeError(
                    f"variant {item['variant']} failed; see {case_dir / 'runner.log'}")
    else:
        run_rows = [{**item, "returncode": None, "status": "skipped"} for item in config_rows]

    mesh = Path(base["mesh_file"])
    if not mesh.is_absolute():
        mesh = base_dir / mesh
    nodes = read_nodes(mesh)
    node_count = len(nodes)
    variant_by_tag = {variant.tag: variant for variant in variants}
    if "baseline" not in variant_by_tag:
        raise RuntimeError("baseline variant is required for comparison")

    comparison_rows: list[dict[str, Any]] = []
    for bias in biases:
        baseline_state = load_variant_state(cases_dir / "baseline", bias, node_count)
        sentaurus_state = branch_profiles.load_sentaurus_state(
            branch_profiles.sentaurus_dir(args.sentaurus_root, bias), node_count)
        for variant in variants:
            candidate_state = load_variant_state(cases_dir / variant.tag, bias, node_count)
            for band, x_min, x_max in bands:
                ids = branch_profiles.band_node_ids(nodes, x_min, x_max)
                comparison_rows.append(comparison_row(
                    variant,
                    bias,
                    band,
                    ids,
                    candidate_state,
                    sentaurus_state,
                    baseline_state,
                ))

    write_csv(
        args.out_dir / "noimpact_variant_branch_comparison.csv",
        comparison_rows,
        SUMMARY_FIELDS,
    )
    (args.out_dir / "noimpact_variant_scan_summary.json").write_text(
        json.dumps(clean_json({
            "base_config": str(base_config),
            "biases": biases,
            "bands": [
                {"name": name, "x_min_um": x_min, "x_max_um": x_max}
                for name, x_min, x_max in bands
            ],
            "variants": run_rows,
            "comparison_rows": len(comparison_rows),
        }), indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
