#!/usr/bin/env python3
"""Compute ideality factors and IV slope comparison for pn2d anchor probe."""
import csv
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
VD = REPO / "build" / "pn2d_probe" / "vela"
RD = REPO / "build" / "pn2d_probe" / "reference_curves"
Vt = 0.02585  # kT/q at 300K


def f(row: dict, k: str) -> float:
    try:
        return float(row.get(k, "") or "nan")
    except (TypeError, ValueError):
        return float("nan")


def ideality_sequence(pts: list[tuple[float, float]], label: str) -> None:
    """Print ideality factors for adjacent pairs."""
    print(f"{label} ideality factor n = dV / (Vt * d·ln|I|):")
    prev = None
    for v, i in sorted(pts):
        if v < 0.05 or math.isnan(i) or i <= 0:
            prev = (v, i)
            continue
        if prev:
            vp, ip = prev
            if ip > 0 and i > 0:
                n = (v - vp) / (Vt * math.log(i / ip))
                print(f"  V={vp:.3f}->{v:.3f}: n={n:.4f}  I={i:.4e}")
        prev = (v, i)
    print()


def main() -> None:
    iv = [r for r in csv.DictReader(open(VD / "pn2d_iv_anchor_probe.csv"))
          if f(r, "converged") == 1]
    ref = list(csv.DictReader(open(RD / "pn2d_iv_reference.csv")))

    vela_pts = [(f(r, "bias_V"), abs(f(r, "current_total_A_per_um"))) for r in iv]
    ref_pts = [(f(r, "bias_V"), abs(f(r, "current_total"))) for r in ref]

    ideality_sequence(vela_pts, "Vela")
    # Sentaurus reference may have finer bias steps; sample every ~0.05 V
    ref_sampled: list[tuple[float, float]] = []
    last_v = -999.0
    for v, i in sorted(ref_pts):
        if v - last_v >= 0.04:
            ref_sampled.append((v, i))
            last_v = v
    ideality_sequence(ref_sampled, "Sentaurus")

    # Orders gap vs bias
    print("Orders gap (log10|I_sent/I_vela|) vs bias:")
    print(f"{'Bias':>6}  {'Vela I':>14}  {'Sent I':>14}  {'orders':>8}")
    from bisect import bisect_left
    ref_sorted = sorted(ref_pts)
    for v, vi in sorted(vela_pts):
        if v < 0.1 or vi <= 0:
            continue
        # Interpolate sentaurus
        xs = [x for x, _ in ref_sorted]
        idx = bisect_left(xs, v)
        if idx == 0:
            si = ref_sorted[0][1]
        elif idx >= len(xs):
            si = ref_sorted[-1][1]
        else:
            x0, y0 = ref_sorted[idx - 1]
            x1, y1 = ref_sorted[idx]
            si = y0 + (y1 - y0) * (v - x0) / (x1 - x0) if x1 != x0 else y0
        if si > 0:
            ord_ = math.log10(max(si, 1e-300) / max(vi, 1e-300))
            print(f"{v:>6.3f}  {vi:>14.4e}  {si:>14.4e}  {ord_:>8.4f}")


if __name__ == "__main__":
    main()
