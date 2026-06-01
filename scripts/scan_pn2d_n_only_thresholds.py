#!/usr/bin/env python3
"""Round-2 pn2d n-contact-only threshold sweep without reference regeneration."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT_DIR = REPO / "build" / "pn2d_contact_relax_scan"
REF_DIR = OUT_DIR / "reference"
CAND_DIR = OUT_DIR / "candidates_round2"


def run(cmd: list[str], cwd: Path) -> None:
    print("[run]", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def at(rows: list[dict[str, str]], bias: float, col: str) -> float:
    for row in rows:
        if abs(float(row["bias_V"]) - bias) <= 1.0e-9:
            return float(row[col])
    raise RuntimeError(f"Missing bias {bias} in {col}")


def interp(rows: list[dict[str, str]], bias: float, col: str) -> float:
    pts = sorted((float(r["bias_V"]), float(r[col])) for r in rows)
    if not pts:
        raise RuntimeError(f"No rows for {col}")
    if bias < pts[0][0] or bias > pts[-1][0]:
        raise RuntimeError(f"Bias {bias} outside range for {col}")
    for x, y in pts:
        if abs(x - bias) <= 1.0e-12:
            return y
    for i in range(1, len(pts)):
        x0, y0 = pts[i - 1]
        x1, y1 = pts[i]
        if x0 <= bias <= x1:
            if abs(x1 - x0) <= 1.0e-12:
                return y0
            t = (bias - x0) / (x1 - x0)
            return y0 + t * (y1 - y0)
    raise RuntimeError(f"Could not interpolate {col} at {bias}")


def strict(rows: list[dict[str, str]]) -> bool:
    for row in rows:
        if row.get("solver_method") != "gummel_newton":
            return False
        if row.get("handoff_stage") != "newton":
            return False
        if int(row.get("newton_iterations", "0")) <= 0:
            return False
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run pn2d n-only threshold scans without regenerating references.")
    parser.add_argument(
        "--thresholds",
        default="0.08,0.10,0.12,0.15,0.20",
        help="Comma-separated threshold list in V (default: 0.08,0.10,0.12,0.15,0.20).",
    )
    parser.add_argument(
        "--strengths",
        default="",
        help="Optional comma-separated relaxation strengths in [0,1] to sweep instead of thresholds.",
    )
    parser.add_argument(
        "--strength-threshold",
        type=float,
        default=0.1,
        help="Bias threshold used when sweeping strengths (default: 0.1 V).",
    )
    parser.add_argument(
        "--summary-prefix",
        default="pn2d_contact_relax_round2_n_only_summary",
        help="Summary output basename under build/pn2d_contact_relax_scan.",
    )
    parser.add_argument(
        "--candidate-dirname",
        default="candidates_round2",
        help="Candidate subdirectory name under build/pn2d_contact_relax_scan.",
    )
    return parser.parse_args()


def parse_thresholds(raw: str) -> list[float]:
    values: list[float] = []
    for token in raw.split(","):
        text = token.strip()
        if not text:
            continue
        value = float(text)
        if value < 0.0:
            raise ValueError("thresholds must be non-negative")
        values.append(value)
    if not values:
        raise ValueError("threshold list is empty")
    return values


def parse_strengths(raw: str) -> list[float]:
    values: list[float] = []
    for token in raw.split(","):
        text = token.strip()
        if not text:
            continue
        value = float(text)
        if value < 0.0 or value > 1.0:
            raise ValueError("strengths must lie in [0, 1]")
        values.append(value)
    if not values:
        raise ValueError("strength list is empty")
    return values


def main() -> int:
    args = parse_args()
    ref_iv = REF_DIR / "reference_curves" / "pn2d_iv_reference.csv"
    ref_bv = REF_DIR / "reference_curves" / "pn2d_bv_reference.csv"
    if not ref_iv.exists() or not ref_bv.exists():
        raise FileNotFoundError("Reference artifacts missing under build/pn2d_contact_relax_scan/reference")

    base_iv = json.loads((REF_DIR / "vela" / "simulation_iv.json").read_text())
    base_bv = json.loads((REF_DIR / "vela" / "simulation_bv.json").read_text())
    runner = REPO / "build" / "vela_example_runner.exe"
    cmp_script = REPO / "scripts" / "compare_reference_curves.py"

    cand_dir = OUT_DIR / args.candidate_dirname
    cand_dir.mkdir(parents=True, exist_ok=True)

    ref_rows = read_csv(ref_iv)
    ref_ratio = abs(interp(ref_rows, 0.29, "current_total") / interp(ref_rows, 0.30, "current_total"))

    strengths = parse_strengths(args.strengths) if args.strengths.strip() else []
    thresholds = [] if strengths else parse_thresholds(args.thresholds)
    print(f"thresholds: {thresholds}", flush=True)
    print(f"strengths: {strengths}", flush=True)
    print(f"strength_threshold: {args.strength_threshold}", flush=True)
    print(f"candidate_dir: {cand_dir}", flush=True)
    print(f"summary_prefix: {args.summary_prefix}", flush=True)
    rows: list[dict[str, object]] = []

    baseline = {
        "tag": "baseline",
        "threshold": None,
        "solver_override": {},
    }
    if strengths:
        candidates = [baseline] + [
            {
                "tag": f"n_only_str{str(st).replace('.', 'p')}",
                "threshold": args.strength_threshold,
                "strength": st,
                "solver_override": {
                    "contact_boundary_reconstruction": "dominant_signed_contact_mean",
                    "contact_boundary_minority_electron_relaxation": True,
                    "contact_boundary_minority_electron_relaxation_contact_side": "n_contact_only",
                    "contact_boundary_minority_electron_relaxation_bias_threshold_V": args.strength_threshold,
                    "contact_boundary_minority_electron_relaxation_strength": st,
                },
            }
            for st in strengths
        ]
    else:
        candidates = [baseline] + [
            {
                "tag": f"n_only_th{str(th).replace('.', 'p')}",
                "threshold": th,
                "solver_override": {
                    "contact_boundary_reconstruction": "dominant_signed_contact_mean",
                    "contact_boundary_minority_electron_relaxation": True,
                    "contact_boundary_minority_electron_relaxation_contact_side": "n_contact_only",
                    "contact_boundary_minority_electron_relaxation_bias_threshold_V": th,
                },
            }
            for th in thresholds
        ]

    for cand in candidates:
        tag = str(cand["tag"])
        print(f"\n=== candidate: {tag} ===", flush=True)
        cdir = cand_dir / tag
        cdir.mkdir(parents=True, exist_ok=True)

        iv = json.loads(json.dumps(base_iv))
        bv = json.loads(json.dumps(base_bv))

        for deck in (iv, bv):
            deck["mesh_file"] = str((REF_DIR / "vela" / "mesh.json").resolve())
            deck["node_doping_file"] = str((REF_DIR / "vela" / "doping.csv").resolve())
            deck.setdefault("solver", {}).update(json.loads(json.dumps(cand["solver_override"])))

        iv["output_csv"] = f"pn2d_iv_{tag}.csv"
        bv["output_csv"] = f"pn2d_bv_{tag}.csv"

        iv_cfg = cdir / f"simulation_iv_{tag}.json"
        bv_cfg = cdir / f"simulation_bv_{tag}.json"
        iv_cfg.write_text(json.dumps(iv, indent=2) + "\n")
        bv_cfg.write_text(json.dumps(bv, indent=2) + "\n")

        run([str(runner), "--config", str(iv_cfg)], cdir)
        run([str(runner), "--config", str(bv_cfg)], cdir)

        fine_c = json.loads(json.dumps(iv))
        fine_a = json.loads(json.dumps(iv))
        fine_c["output_csv"] = f"pn2d_iv_fine_cathode_{tag}.csv"
        fine_a["output_csv"] = f"pn2d_iv_fine_anode_{tag}.csv"
        for deck, cc in ((fine_c, "Cathode"), (fine_a, "Anode")):
            deck["sweep"]["start"] = 0.0
            deck["sweep"]["stop"] = 0.3
            deck["sweep"]["step"] = 0.01
            deck["sweep"]["current_contact"] = cc

        fc_cfg = cdir / f"simulation_iv_fine_cathode_{tag}.json"
        fa_cfg = cdir / f"simulation_iv_fine_anode_{tag}.json"
        fc_cfg.write_text(json.dumps(fine_c, indent=2) + "\n")
        fa_cfg.write_text(json.dumps(fine_a, indent=2) + "\n")

        run([str(runner), "--config", str(fc_cfg)], cdir)
        run([str(runner), "--config", str(fa_cfg)], cdir)

        iv_csv = cdir / iv["output_csv"]
        bv_csv = cdir / bv["output_csv"]
        fc_csv = cdir / fine_c["output_csv"]
        fa_csv = cdir / fine_a["output_csv"]

        iv_cmp = cdir / f"pn2d_iv_window_{tag}.json"
        run(
            [
                sys.executable,
                str(cmp_script),
                "--reference",
                str(ref_iv),
                "--candidate",
                str(iv_csv),
                "--output-json",
                str(iv_cmp),
                "--output-md",
                str(cdir / f"pn2d_iv_window_{tag}.md"),
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
            REPO,
        )

        bv_cmp = cdir / f"pn2d_bv_0p05_{tag}.json"
        run(
            [
                sys.executable,
                str(cmp_script),
                "--reference",
                str(ref_bv),
                "--candidate",
                str(bv_csv),
                "--output-json",
                str(bv_cmp),
                "--output-md",
                str(cdir / f"pn2d_bv_0p05_{tag}.md"),
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
            REPO,
        )

        iv_cmp_data = json.loads(iv_cmp.read_text())
        bv_cmp_data = json.loads(bv_cmp.read_text())
        iv_rows = read_csv(iv_csv)
        bv_rows = read_csv(bv_csv)
        fc_rows = read_csv(fc_csv)
        fa_rows = read_csv(fa_csv)

        ratio = abs(
            at(fc_rows, 0.29, "current_total_A_per_um")
            / at(fc_rows, 0.30, "current_total_A_per_um")
        )
        terminal_sum = abs(
            at(fc_rows, 0.30, "current_total_A_per_um")
            + at(fa_rows, 0.30, "current_total_A_per_um")
        )

        row = {
            "tag": tag,
            "threshold": cand["threshold"],
            "candidate_ratio_0p29_over_0p30": ratio,
            "ratio_abs_delta_to_reference": abs(ratio - ref_ratio),
            "iv_window_orders_0p2_to_0p3": iv_cmp_data["iv"]["orders_of_magnitude"],
            "iv_window_trend_match": iv_cmp_data["iv"].get("trend_match"),
            "terminal_sum_abs_A_per_um_at_0p3": terminal_sum,
            "bv_orders_at_0p05": bv_cmp_data["iv"]["orders_of_magnitude"],
            "strict_newton_handoff_all": strict(iv_rows) and strict(bv_rows) and strict(fc_rows) and strict(fa_rows),
        }
        rows.append(row)
        print(
            f"{tag}: ratio={row['candidate_ratio_0p29_over_0p30']:.6f}, "
            f"delta={row['ratio_abs_delta_to_reference']:.6f}, "
            f"iv_orders={row['iv_window_orders_0p2_to_0p3']}, "
            f"sum0p3={row['terminal_sum_abs_A_per_um_at_0p3']:.3e}, "
            f"bv_orders={row['bv_orders_at_0p05']}, "
            f"strict_newton={row['strict_newton_handoff_all']}",
            flush=True,
        )

    out_csv = OUT_DIR / f"{args.summary_prefix}.csv"
    out_json = OUT_DIR / f"{args.summary_prefix}.json"
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    out_json.write_text(json.dumps(rows, indent=2) + "\n")

    print("\nSummary CSV:", out_csv, flush=True)
    print("Summary JSON:", out_json, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
