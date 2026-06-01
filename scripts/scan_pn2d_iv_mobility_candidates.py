#!/usr/bin/env python3
"""Phase 2A: scan pn2d IV mobility-model candidates with fixed guardrails.

This script regenerates (or reuses) a fresh pn2d reference workspace via
sentaurus_import, then evaluates a focused IV mobility candidate matrix by
overriding only the solver mobility block. The contact-boundary reconstruction
baseline from the reference deck is preserved so the closed 0.3 V terminal-sum
gate is not disturbed. It reports, per candidate:
- I(0.29)/I(0.30)
- IV window orders (0.2-0.3 V)
- terminal current sum at 0.3 V
- BV orders at 0.05 V (guardrail; BV mobility block is left untouched)
- strict Newton handoff status

Sentaurus pn2d IV physics declares DopingDep + HighFieldSaturation + Mobility,
so caughey_thomas and caughey_thomas_field variants are the faithful targets.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]

# Sentaurus-faithful Caughey-Thomas constants used by the promoted BV deck.
BV_CT_CONSTANTS = {
    "electron_mu_min_m2_V_s": 46.458,
    "hole_mu_min_m2_V_s": 39.961,
    "electron_nref_m3": 9.68e16,
    "hole_nref_m3": 2.23e17,
    "electron_alpha": 0.6052,
    "hole_alpha": 0.623,
}


@dataclass(frozen=True)
class Candidate:
    tag: str
    description: str
    mobility: dict[str, Any] | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO / "build" / "pn2d_iv_mobility_scan",
        help="Output directory for generated decks, curves, and reports.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=REPO / "reference_tcad" / "pn2d" / "pn2d_reference.json",
        help="pn2d reference config path.",
    )
    parser.add_argument(
        "--source-dir",
        type=Path,
        default=REPO / "reference_tcad" / "pn2d",
        help="Sentaurus source directory.",
    )
    parser.add_argument(
        "--runner",
        type=Path,
        default=None,
        help="Path to vela_example_runner executable.",
    )
    parser.add_argument(
        "--tdr-importer",
        type=Path,
        default=None,
        help="Path to sentaurus_import executable.",
    )
    parser.add_argument(
        "--keep-output",
        action="store_true",
        help="Keep existing output dir contents instead of deleting it first.",
    )
    parser.add_argument(
        "--skip-reference-generation",
        action="store_true",
        help="Reuse existing reference artifacts under output-dir/reference.",
    )
    return parser.parse_args()


def detect_runner(build_dir: Path) -> Path:
    exe = "vela_example_runner.exe" if sys.platform.startswith("win") else "vela_example_runner"
    for candidate in [build_dir / exe, build_dir / "src" / "tools" / exe]:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError("Could not find vela_example_runner in build output.")


def detect_importer(build_dir: Path) -> Path:
    exe = "sentaurus_import.exe" if sys.platform.startswith("win") else "sentaurus_import"
    candidate = build_dir / exe
    if candidate.is_file():
        return candidate
    raise FileNotFoundError("Could not find sentaurus_import executable in build output.")


def run(cmd: list[str], cwd: Path) -> None:
    print("[run]", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)


def read_curve(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def value_at_bias(rows: list[dict[str, str]], bias: float, column: str, tol: float = 1.0e-9) -> float:
    for row in rows:
        if abs(float(row["bias_V"]) - bias) <= tol:
            return float(row[column])
    raise ValueError(f"Missing bias {bias} in {column}")


def interp_at_bias(rows: list[dict[str, str]], bias: float, column: str, tol: float = 1.0e-12) -> float:
    points = sorted((float(row["bias_V"]), float(row[column])) for row in rows)
    if not points:
        raise ValueError(f"Curve has no rows for {column}")
    if bias < points[0][0] - tol or bias > points[-1][0] + tol:
        raise ValueError(f"Bias {bias} outside curve range for {column}")
    for v, y in points:
        if abs(v - bias) <= tol:
            return y
    for idx in range(1, len(points)):
        v0, y0 = points[idx - 1]
        v1, y1 = points[idx]
        if v0 <= bias <= v1:
            if abs(v1 - v0) <= tol:
                return y0
            alpha = (bias - v0) / (v1 - v0)
            return y0 + alpha * (y1 - y0)
    raise ValueError(f"Could not interpolate {column} at {bias}")


def strict_newton_handoff(rows: list[dict[str, str]]) -> bool:
    for row in rows:
        if row.get("solver_method") != "gummel_newton":
            return False
        if row.get("handoff_stage") != "newton":
            return False
        if int(row.get("newton_iterations", "0")) <= 0:
            return False
    return True


def apply_mobility_override(deck: dict[str, Any], mobility: dict[str, Any] | None) -> dict[str, Any]:
    merged = json.loads(json.dumps(deck))
    solver = merged.setdefault("solver", {})
    if mobility is None:
        solver.pop("mobility", None)
    else:
        solver["mobility"] = json.loads(json.dumps(mobility))
    return merged


def resolve_input_paths(deck: dict[str, Any], base_dir: Path) -> dict[str, Any]:
    resolved = json.loads(json.dumps(deck))
    path_keys = [
        "mesh_file",
        "node_doping_file",
        "node_mobility_file",
        "node_lifetime_file",
        "interface_charge_file",
    ]
    for key in path_keys:
        raw = resolved.get(key)
        if not isinstance(raw, str):
            continue
        candidate = Path(raw)
        if candidate.is_absolute():
            continue
        from_base = (base_dir / candidate).resolve()
        if from_base.exists():
            resolved[key] = str(from_base)
    return resolved


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n")


def evaluate_candidate(
    candidate: Candidate,
    candidate_dir: Path,
    base_iv_deck: dict[str, Any],
    base_bv_deck: dict[str, Any],
    base_deck_dir: Path,
    ref_iv: Path,
    ref_bv: Path,
    runner: Path,
    python_exe: str,
) -> dict[str, Any]:
    candidate_dir.mkdir(parents=True, exist_ok=True)

    # Only the IV mobility block changes; BV deck is left untouched as a guardrail.
    iv_deck = resolve_input_paths(
        apply_mobility_override(base_iv_deck, candidate.mobility),
        base_deck_dir,
    )
    bv_deck = resolve_input_paths(base_bv_deck, base_deck_dir)

    iv_deck["output_csv"] = f"pn2d_iv_{candidate.tag}.csv"
    bv_deck["output_csv"] = f"pn2d_bv_{candidate.tag}.csv"

    iv_cfg = candidate_dir / f"simulation_iv_{candidate.tag}.json"
    bv_cfg = candidate_dir / f"simulation_bv_{candidate.tag}.json"
    write_json(iv_cfg, iv_deck)
    write_json(bv_cfg, bv_deck)

    run([str(runner), "--config", str(iv_cfg)], cwd=candidate_dir)
    run([str(runner), "--config", str(bv_cfg)], cwd=candidate_dir)

    fine_cath = json.loads(json.dumps(iv_deck))
    fine_anode = json.loads(json.dumps(iv_deck))
    fine_cath["output_csv"] = f"pn2d_iv_fine_cathode_{candidate.tag}.csv"
    fine_anode["output_csv"] = f"pn2d_iv_fine_anode_{candidate.tag}.csv"
    fine_cath["sweep"]["start"] = 0.0
    fine_cath["sweep"]["stop"] = 0.3
    fine_cath["sweep"]["step"] = 0.01
    fine_anode["sweep"]["start"] = 0.0
    fine_anode["sweep"]["stop"] = 0.3
    fine_anode["sweep"]["step"] = 0.01
    fine_cath["sweep"]["current_contact"] = "Cathode"
    fine_anode["sweep"]["current_contact"] = "Anode"

    fine_cath_cfg = candidate_dir / f"simulation_iv_fine_cathode_{candidate.tag}.json"
    fine_anode_cfg = candidate_dir / f"simulation_iv_fine_anode_{candidate.tag}.json"
    write_json(fine_cath_cfg, fine_cath)
    write_json(fine_anode_cfg, fine_anode)
    run([str(runner), "--config", str(fine_cath_cfg)], cwd=candidate_dir)
    run([str(runner), "--config", str(fine_anode_cfg)], cwd=candidate_dir)

    iv_csv = candidate_dir / iv_deck["output_csv"]
    bv_csv = candidate_dir / bv_deck["output_csv"]
    fine_cath_csv = candidate_dir / fine_cath["output_csv"]
    fine_anode_csv = candidate_dir / fine_anode["output_csv"]

    iv_rows = read_curve(iv_csv)
    bv_rows = read_curve(bv_csv)
    fine_cath_rows = read_curve(fine_cath_csv)
    fine_anode_rows = read_curve(fine_anode_csv)

    candidate_ratio = abs(
        value_at_bias(fine_cath_rows, 0.29, "current_total_A_per_um")
        / value_at_bias(fine_cath_rows, 0.30, "current_total_A_per_um")
    )
    terminal_sum_0p3 = abs(
        value_at_bias(fine_cath_rows, 0.30, "current_total_A_per_um")
        + value_at_bias(fine_anode_rows, 0.30, "current_total_A_per_um")
    )

    iv_cmp_json = candidate_dir / f"pn2d_iv_window_{candidate.tag}.json"
    iv_cmp_md = candidate_dir / f"pn2d_iv_window_{candidate.tag}.md"
    run(
        [
            python_exe,
            str(REPO / "scripts" / "compare_reference_curves.py"),
            "--reference",
            str(ref_iv),
            "--candidate",
            str(iv_csv),
            "--output-json",
            str(iv_cmp_json),
            "--output-md",
            str(iv_cmp_md),
            "--kind",
            "iv",
            "--candidate-scale",
            "-1",
            "--candidate-column",
            "current_total_A_per_um",
            "--bias-min",
            "0.2",
            "--bias-max",
            "0.3",
            "--min-points",
            "2",
            "--require-trend-match",
        ],
        cwd=REPO,
    )
    bv_cmp_json = candidate_dir / f"pn2d_bv_0p05_{candidate.tag}.json"
    bv_cmp_md = candidate_dir / f"pn2d_bv_0p05_{candidate.tag}.md"
    run(
        [
            python_exe,
            str(REPO / "scripts" / "compare_reference_curves.py"),
            "--reference",
            str(ref_bv),
            "--candidate",
            str(bv_csv),
            "--output-json",
            str(bv_cmp_json),
            "--output-md",
            str(bv_cmp_md),
            "--kind",
            "iv",
            "--candidate-scale",
            "1",
            "--candidate-column",
            "current_total_A_per_um",
            "--bias-min",
            "0.05",
            "--bias-max",
            "0.05",
            "--min-points",
            "1",
        ],
        cwd=REPO,
    )

    iv_cmp = json.loads(iv_cmp_json.read_text())
    bv_cmp = json.loads(bv_cmp_json.read_text())

    return {
        "tag": candidate.tag,
        "description": candidate.description,
        "candidate_ratio_0p29_over_0p30": candidate_ratio,
        "iv_window_orders_0p2_to_0p3": iv_cmp["iv"].get("orders_of_magnitude"),
        "iv_window_trend_match": iv_cmp["iv"].get("trend_match"),
        "terminal_sum_abs_A_per_um_at_0p3": terminal_sum_0p3,
        "bv_orders_at_0p05": bv_cmp["iv"].get("orders_of_magnitude"),
        "strict_newton_handoff_iv": strict_newton_handoff(iv_rows),
        "strict_newton_handoff_bv": strict_newton_handoff(bv_rows),
        "strict_newton_handoff_fine_cathode": strict_newton_handoff(fine_cath_rows),
        "strict_newton_handoff_fine_anode": strict_newton_handoff(fine_anode_rows),
    }


def main() -> int:
    args = parse_args()
    build_dir = REPO / "build"
    runner = args.runner.resolve() if args.runner else detect_runner(build_dir)
    importer = args.tdr_importer.resolve() if args.tdr_importer else detect_importer(build_dir)

    out_dir = args.output_dir.resolve()
    if out_dir.exists() and not args.keep_output:
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    reference_dir = out_dir / "reference"
    if not args.skip_reference_generation:
        run(
            [
                sys.executable,
                str(REPO / "scripts" / "sentaurus_import.py"),
                "reference",
                "--config",
                str(args.config.resolve()),
                "--source-dir",
                str(args.source_dir.resolve()),
                "--output-dir",
                str(reference_dir),
                "--tdr-importer",
                str(importer),
                "--runner",
                str(runner),
            ],
            cwd=REPO,
        )

    ref_iv = reference_dir / "reference_curves" / "pn2d_iv_reference.csv"
    ref_bv = reference_dir / "reference_curves" / "pn2d_bv_reference.csv"
    if args.skip_reference_generation and (not ref_iv.exists() or not ref_bv.exists()):
        raise FileNotFoundError(
            "--skip-reference-generation requested but reference artifacts are missing."
        )
    base_iv_deck = json.loads((reference_dir / "vela" / "simulation_iv.json").read_text())
    base_bv_deck = json.loads((reference_dir / "vela" / "simulation_bv.json").read_text())
    base_deck_dir = reference_dir / "vela"

    reference_rows = read_curve(ref_iv)
    reference_ratio = abs(
        interp_at_bias(reference_rows, 0.29, "current_total")
        / interp_at_bias(reference_rows, 0.30, "current_total")
    )

    ct_field = dict(BV_CT_CONSTANTS)
    candidates = [
        Candidate("baseline", "Reference IV deck mobility (default model)", None),
        Candidate(
            "ct_silicon",
            "caughey_thomas with default silicon constants",
            {"model": "caughey_thomas"},
        ),
        Candidate(
            "ct_field_silicon",
            "caughey_thomas_field (DopingDep + HighFieldSaturation), silicon defaults",
            {"model": "caughey_thomas_field"},
        ),
        Candidate(
            "ct_bv_constants",
            "caughey_thomas with promoted BV Sentaurus constants",
            {"model": "caughey_thomas", **BV_CT_CONSTANTS},
        ),
        Candidate(
            "ct_field_bv_constants",
            "caughey_thomas_field with promoted BV Sentaurus constants",
            {"model": "caughey_thomas_field", **ct_field},
        ),
    ]

    rows: list[dict[str, Any]] = []
    for candidate in candidates:
        candidate_dir = out_dir / "candidates" / candidate.tag
        metrics = evaluate_candidate(
            candidate,
            candidate_dir,
            base_iv_deck,
            base_bv_deck,
            base_deck_dir,
            ref_iv,
            ref_bv,
            runner,
            sys.executable,
        )
        metrics["reference_ratio_0p29_over_0p30"] = reference_ratio
        metrics["ratio_abs_delta_to_reference"] = abs(
            metrics["candidate_ratio_0p29_over_0p30"] - reference_ratio
        )
        metrics["strict_newton_handoff_all"] = all(
            bool(metrics[key])
            for key in [
                "strict_newton_handoff_iv",
                "strict_newton_handoff_bv",
                "strict_newton_handoff_fine_cathode",
                "strict_newton_handoff_fine_anode",
            ]
        )
        rows.append(metrics)

    summary_csv = out_dir / "pn2d_iv_mobility_summary.csv"
    summary_json = out_dir / "pn2d_iv_mobility_summary.json"
    fieldnames = list(rows[0].keys()) if rows else []
    with summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    summary_json.write_text(json.dumps(rows, indent=2) + "\n")

    print("\n=== pn2d IV mobility candidate summary ===")
    print(f"reference ratio I(0.29)/I(0.30) = {reference_ratio:.6f}")
    for row in rows:
        print(
            f"{row['tag']}: ratio={row['candidate_ratio_0p29_over_0p30']:.6f}, "
            f"delta={row['ratio_abs_delta_to_reference']:.6f}, "
            f"iv_orders={row['iv_window_orders_0p2_to_0p3']}, "
            f"sum0p3={row['terminal_sum_abs_A_per_um_at_0p3']:.3e}, "
            f"bv_orders={row['bv_orders_at_0p05']}, "
            f"strict_newton={row['strict_newton_handoff_all']}"
        )
    print(f"\nSummary CSV: {summary_csv}")
    print(f"Summary JSON: {summary_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
