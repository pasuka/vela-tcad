#!/usr/bin/env python3
"""Full-range and segmented log10-magnitude IV RMS: Vela -20V BV sweep vs Sentaurus reference."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def load(path: Path, bcol: str, ccol: str, scale: float = 1.0) -> dict[float, float]:
    out: dict[float, float] = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            try:
                b = round(float(row[bcol]), 4)
                c = float(row[ccol]) * scale
            except (ValueError, KeyError):
                continue
            out[b] = c
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--candidate", type=Path, required=True)
    ap.add_argument("--floor", type=float, default=1.0e-25)
    ap.add_argument("--match-tol", type=float, default=0.06)
    args = ap.parse_args()

    ref = load(args.reference, "bias_V", "current_total")
    vel = load(args.candidate, "bias_V", "current_total_A_per_um")
    vb = sorted(vel.keys())

    def nearest(b: float) -> float:
        return min(vb, key=lambda x: abs(x - b))

    rows: list[tuple[float, float, float, float, float]] = []
    for b in sorted(ref.keys()):
        if b >= 0:
            continue
        nb = nearest(b)
        if abs(nb - b) > args.match_tol:
            continue
        rs, vs = ref[b], vel[nb]
        if abs(rs) < args.floor or abs(vs) < args.floor:
            continue
        e = math.log10(abs(vs)) - math.log10(abs(rs))
        rows.append((b, rs, vs, e, vs / rs))

    if not rows:
        print("no matched points")
        return

    def rms(lo: float, hi: float) -> tuple[float, int]:
        xs = [r[3] for r in rows if lo >= r[0] >= hi]
        if not xs:
            return float("nan"), 0
        return math.sqrt(sum(x * x for x in xs) / len(xs)), len(xs)

    allr = [r[3] for r in rows]
    print(f"matched points: {len(rows)}  bias range {rows[0][0]} .. {rows[-1][0]}")
    full = math.sqrt(sum(x * x for x in allr) / len(allr))
    print(f"FULL  (0..-20)   log10 RMS = {full:.4f} decades")
    for lo, hi, name in [
        (-0.1, -5.0, "0..-5"),
        (-5.0, -10.0, "-5..-10"),
        (-10.0, -13.2, "-10..-13.2"),
        (-13.2, -20.0, "-13.2..-20"),
    ]:
        r, n = rms(lo, hi)
        print(f"  seg {name:12s} RMS={r:.4f}  (n={n})")
    print("--- sample points ---")
    for target in [-1, -5, -10, -12, -13.2, -15, -18, -20]:
        cand = [r for r in rows if abs(r[0] - target) < 0.051]
        if cand:
            r = cand[0]
            print(
                f"  {r[0]:7.2f}V  Sent={r[1]:.3e}  Vela={r[2]:.3e}  "
                f"ratio={r[4]:.3f}  dlog10={r[3]:+.3f}"
            )


if __name__ == "__main__":
    main()
