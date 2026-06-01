#!/usr/bin/env python3
"""Iteration C analysis: isolate the I0 gap into ni and BGN contributions.

2x2 factorial (all no-SRH):
  c0: ni=1.0e10, slotboom
  c1: ni=1.45e10, slotboom
  c2: ni=1.0e10, no BGN
  c3: ni=1.45e10, no BGN
"""
import csv
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
VD = REPO / "build" / "pn2d_probe" / "vela"
RD = REPO / "build" / "pn2d_probe" / "reference_curves"

VARIANTS = [
    ("c0_ni10_slotboom",  "ni=1.0e10, slotboom"),
    ("c1_ni145_slotboom", "ni=1.45e10, slotboom"),
    ("c2_ni10_nobgn",     "ni=1.0e10, no BGN"),
    ("c3_ni145_nobgn",    "ni=1.45e10, no BGN"),
]
BIASES = [0.15, 0.20, 0.25, 0.27, 0.29, 0.30]


def f(row: dict, k: str) -> float:
    try:
        return float(row.get(k, "") or "nan")
    except (TypeError, ValueError):
        return float("nan")


def ref_interp(ref_rows: list, bias: float) -> float:
    pts = sorted((float(r["bias_V"]), float(r["current_total"])) for r in ref_rows)
    if bias <= pts[0][0]:
        return pts[0][1]
    if bias >= pts[-1][0]:
        return pts[-1][1]
    for i in range(1, len(pts)):
        x0, y0 = pts[i - 1]
        x1, y1 = pts[i]
        if x0 <= bias <= x1:
            t = (bias - x0) / (x1 - x0) if x1 != x0 else 0.0
            return y0 + t * (y1 - y0)
    return math.nan


def load(tag: str) -> dict:
    fpath = VD / f"pn2d_iv_{tag}.csv"
    if not fpath.exists():
        return {}
    rows = list(csv.DictReader(fpath.open()))
    return {round(f(r, "bias_V"), 3): abs(f(r, "current_total_A_per_um"))
            for r in rows if f(r, "converged") == 1}


def main() -> None:
    ref_rows = list(csv.DictReader(open(RD / "pn2d_iv_reference.csv")))
    data = {tag: load(tag) for tag, _ in VARIANTS}

    print("=" * 90)
    print("ITERATION C: I0 gap decomposition (all variants no-SRH) — current in A/um")
    print("=" * 90)
    hdr = f"{'Bias':>5}  {'Sentaurus':>12}  " + "  ".join(f"{lbl[:18]:>18}" for _, lbl in VARIANTS)
    print(hdr)
    for bias in BIASES:
        si = abs(ref_interp(ref_rows, bias))
        vals = []
        for tag, _ in VARIANTS:
            v = data[tag].get(round(bias, 3), math.nan)
            vals.append(v)
        vstr = "  ".join(f"{v:>18.4e}" if not math.isnan(v) else f"{'N/A':>18}" for v in vals)
        print(f"{bias:>5.2f}  {si:>12.4e}  {vstr}")

    print()
    print("Orders gap log10(I_sent/I_vela) at each bias:")
    hdr2 = f"{'Bias':>5}  " + "  ".join(f"{lbl[:18]:>18}" for _, lbl in VARIANTS)
    print(hdr2)
    for bias in BIASES:
        si = abs(ref_interp(ref_rows, bias))
        cells = []
        for tag, _ in VARIANTS:
            v = data[tag].get(round(bias, 3), math.nan)
            if math.isnan(v) or v <= 0:
                cells.append(f"{'N/A':>18}")
            else:
                cells.append(f"{math.log10(si / v):>18.4f}")
        print(f"{bias:>5.2f}  {'  '.join(cells)}")

    # Factorial decomposition at 0.30 V
    print()
    print("=" * 90)
    print("FACTORIAL DECOMPOSITION at 0.30 V (ratio of currents)")
    print("=" * 90)
    b = round(0.30, 3)
    c0 = data["c0_ni10_slotboom"].get(b, math.nan)
    c1 = data["c1_ni145_slotboom"].get(b, math.nan)
    c2 = data["c2_ni10_nobgn"].get(b, math.nan)
    c3 = data["c3_ni145_nobgn"].get(b, math.nan)
    si = abs(ref_interp(ref_rows, 0.30))
    print(f"  Sentaurus I(0.30)               = {si:.4e} A/um")
    print(f"  c0 (ni=1e10,  slotboom)         = {c0:.4e} A/um   gap={si/c0:.2f}x")
    print(f"  c1 (ni=1.45e10, slotboom)       = {c1:.4e} A/um   gap={si/c1:.2f}x")
    print(f"  c2 (ni=1e10,  no BGN)           = {c2:.4e} A/um   gap={si/c2:.2f}x")
    print(f"  c3 (ni=1.45e10, no BGN)         = {c3:.4e} A/um   gap={si/c3:.2f}x")
    print()
    print(f"  ni effect (c1/c0, slotboom)     = {c1/c0:.3f}x   (expected ~2.10x from (1.45/1.0)^2)")
    print(f"  ni effect (c3/c2, no BGN)       = {c3/c2:.3f}x")
    print(f"  BGN effect (c0/c2, ni=1e10)     = {c0/c2:.3f}x   (slotboom raises ni_eff in heavy doping)")
    print(f"  BGN effect (c1/c3, ni=1.45e10)  = {c1/c3:.3f}x")
    print()
    print(f"  Combined best (c1) residual gap = {si/c1:.2f}x  ({math.log10(si/c1):.3f} orders)")


if __name__ == "__main__":
    main()
