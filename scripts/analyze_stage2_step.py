#!/usr/bin/env python3
"""Stage 2 prerequisite: pick the pn2d IV deck vela_step that makes the official
gate reflect physics instead of linear-chord overshoot.

Uses the already-computed fine (step 0.01, 31-point) baseline curve and
sub-samples it at candidate deck steps. Each sub-sample point is a REAL solved
point, so this exactly predicts the gate the regression test would see if the
deck used that step. Gate is computed with the official compare_reference_curves
path.
"""
import csv
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import compare_reference_curves as crc  # noqa: E402

VD = REPO / "build" / "pn2d_probe" / "vela"
REF = REPO / "build" / "pn2d_probe" / "reference_curves" / "pn2d_iv_reference.csv"
SCALE = -1.0
BIAS_MIN, BIAS_MAX = 0.2, 0.3

# Steps to evaluate (must divide the 0.01 grid cleanly)
STEPS = [0.10, 0.05, 0.04, 0.025, 0.02, 0.01]
VARIANTS = [
    ("s0_baseline",     "tau=1e-7,  ni=1.0e10  (current deck)"),
    ("s2_ni145",        "tau=1e-7,  ni=1.45e10"),
    ("s3_tau1e6_ni145", "tau=1e-6,  ni=1.45e10"),
]


def fval(row, k):
    try:
        return float(row.get(k, "") or "nan")
    except (TypeError, ValueError):
        return float("nan")


def load_fine_rows(tag):
    fpath = VD / f"pn2d_iv_{tag}.csv"
    rows = [r for r in csv.DictReader(fpath.open()) if fval(r, "converged") == 1]
    return rows


def subsample(rows, step):
    """Keep rows whose bias is a multiple of `step` (within tolerance)."""
    out = []
    for r in rows:
        b = fval(r, "bias_V")
        n = round(b / step)
        if abs(b - n * step) <= 1.0e-6:
            out.append(r)
    return out


def gate(ref_rows, cand_rows):
    rep = crc.compare_series(
        "iv", ref_rows, cand_rows,
        candidate_scale=SCALE, bias_min=BIAS_MIN, bias_max=BIAS_MAX,
        candidate_column="current_total_A_per_um",
    )
    return rep.get("orders_of_magnitude")


def main():
    ref_rows = list(csv.DictReader(REF.open()))
    print("=" * 92)
    print("STAGE 2 PREREQUISITE — pn2d IV gate vs deck vela_step (sub-sampled real points)")
    print("=" * 92)
    print("Gate = official compare_reference_curves orders_of_magnitude, window [0.2,0.3].")
    print("Each step column is the gate the regression test WOULD see at that deck step.")
    print()
    hdr = f"{'variant':<16}{'desc':<34}" + "".join(f"{('s='+format(s,'.3g')):>9}" for s in STEPS)
    print(hdr)
    print("-" * 92)
    table = {}
    for tag, desc in VARIANTS:
        fine = load_fine_rows(tag)
        cells = []
        for s in STEPS:
            sub = subsample(fine, s)
            g = gate(ref_rows, sub)
            cells.append(g)
        table[tag] = cells
        cellstr = "".join(f"{g:>9.4f}" for g in cells)
        print(f"{tag:<16}{desc:<34}{cellstr}")

    print()
    print("Convergence check: gate should approach the true (fine) value as step -> 0.01.")
    print("The current deck step (0.10) is the leftmost column; note baseline non-monotonicity.")
    print()
    base = table["s0_baseline"]
    fine_val = base[-1]
    print(f"  baseline true (step 0.01) gate = {fine_val:.4f}")
    for s, g in zip(STEPS, base):
        art = g - fine_val
        flag = "  <-- artifact" if abs(art) > 0.05 else ""
        print(f"    step {s:<6} gate={g:.4f}  (artifact vs true {art:+.4f}){flag}")

    print()
    print("Recommended deck step: smallest step where baseline gate is within ~0.02 of")
    print("the true fine value AND ni=1.45e10 clearly ranks below baseline.")


if __name__ == "__main__":
    main()
