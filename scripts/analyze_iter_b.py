#!/usr/bin/env python3
"""Iteration B comparative IV analysis: baseline vs tau variants vs no-SRH."""
import csv
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
VD = REPO / "build" / "pn2d_probe" / "vela"
RD = REPO / "build" / "pn2d_probe" / "reference_curves"
Vt = 0.02585
ANCHORS = [0.27, 0.29, 0.30]

VARIANTS = [
    ("baseline",   "pn2d_iv_anchor_probe.csv",  "tau=1e-7 (default)"),
    ("b1_tau1e6",  "pn2d_iv_b1_tau1e6.csv",     "tau=1e-6 (10x)"),
    ("b2_tau1e5",  "pn2d_iv_b2_tau1e5.csv",     "tau=1e-5 (100x)"),
    ("b3_nosrh",   "pn2d_iv_b3_nosrh.csv",      "no SRH/Auger"),
]


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


def orders(ref: float, cand: float) -> float:
    return math.log10(max(abs(ref), 1e-300) / max(abs(cand), 1e-300))


def ideality(rows: list) -> list[tuple[float, float]]:
    """Return list of (V, n) for adjacent pairs."""
    pts = sorted((f(r, "bias_V"), abs(f(r, "current_total_A_per_um"))) for r in rows
                 if f(r, "converged") == 1)
    result = []
    for i in range(1, len(pts)):
        v0, i0 = pts[i - 1]
        v1, i1 = pts[i]
        if v0 >= 0.05 and i0 > 0 and i1 > 0:
            n = (v1 - v0) / (Vt * math.log(i1 / i0))
            result.append((v1, n))
    return result


def main() -> None:
    ref_rows = list(csv.DictReader(open(RD / "pn2d_iv_reference.csv")))

    # --- Orders gap at anchor biases ---
    print("=" * 80)
    print("ITERATION B: orders gap at anchor biases (log10|I_sent/I_vela|)")
    print("=" * 80)
    header = f"{'Variant':<22}  {'0.27 V':>8}  {'0.29 V':>8}  {'0.30 V':>8}  {'slope Δ':>9}"
    print(header)
    print("-" * 65)

    for tag, fname, label in VARIANTS:
        fpath = VD / fname
        if not fpath.exists():
            print(f"  {label:<22}: FILE NOT FOUND")
            continue
        rows = list(csv.DictReader(fpath.open()))
        ords = {}
        for bias in ANCHORS:
            row = next((r for r in rows if abs(f(r, "bias_V") - bias) < 0.006), None)
            if row:
                vi = abs(f(row, "current_total_A_per_um"))
                ri = abs(ref_interp(ref_rows, bias))
                ords[bias] = orders(ri, vi)
            else:
                ords[bias] = float("nan")
        row29 = next((r for r in rows if abs(f(r, "bias_V") - 0.29) < 0.006), None)
        row30 = next((r for r in rows if abs(f(r, "bias_V") - 0.30) < 0.006), None)
        if row29 and row30:
            i29 = abs(f(row29, "current_total_A_per_um"))
            i30 = abs(f(row30, "current_total_A_per_um"))
            vr = i29 / i30 if i30 else math.nan
            r29 = abs(ref_interp(ref_rows, 0.29))
            r30 = abs(ref_interp(ref_rows, 0.30))
            sr = r29 / r30 if r30 else math.nan
            slope_delta = abs(vr - sr)
        else:
            slope_delta = math.nan
        o27 = f"{ords.get(0.27, math.nan):>8.4f}" if not math.isnan(ords.get(0.27, math.nan)) else "     N/A"
        o29 = f"{ords.get(0.29, math.nan):>8.4f}" if not math.isnan(ords.get(0.29, math.nan)) else "     N/A"
        o30 = f"{ords.get(0.30, math.nan):>8.4f}" if not math.isnan(ords.get(0.30, math.nan)) else "     N/A"
        sd  = f"{slope_delta:>9.6f}" if not math.isnan(slope_delta) else "      N/A"
        print(f"  {label:<22}  {o27}  {o29}  {o30}  {sd}")

    # --- Ideality factor comparison ---
    print()
    print("=" * 80)
    print("IDEALITY FACTOR n = dV/(Vt·d·ln|I|) at key bias points")
    print("=" * 80)
    # Sentaurus reference
    ref_pts = sorted((f(r, "bias_V"), abs(f(r, "current_total"))) for r in ref_rows)
    print("\nSentaurus (reference):")
    prev = None
    for v, i in ref_pts:
        if v < 0.05 or i <= 0:
            prev = (v, i)
            continue
        if prev:
            vp, ip = prev
            if v - vp >= 0.04 and ip > 0:
                n = (v - vp) / (Vt * math.log(i / ip))
                if 0.1 < v <= 0.35:
                    print(f"  V={vp:.3f}->{v:.3f}: n={n:.4f}  I={i:.4e}")
        prev = (v, i)

    for tag, fname, label in VARIANTS:
        fpath = VD / fname
        if not fpath.exists():
            continue
        rows = [r for r in csv.DictReader(fpath.open()) if f(r, "converged") == 1]
        nfactors = ideality(rows)
        print(f"\n{label}:")
        for v, n in nfactors:
            if 0.1 < v <= 0.31:
                print(f"  ->V={v:.3f}: n={n:.4f}")

    # --- Full IV table ---
    print()
    print("=" * 80)
    print("IV TABLE (A/um) — all variants vs Sentaurus")
    print("=" * 80)
    header2 = f"{'Bias':>6}  {'Sentaurus':>14}  " + "  ".join(
        f"{t[2][:14]:>14}" for t in VARIANTS
    )
    print(header2)
    biases = [0.10, 0.15, 0.17, 0.20, 0.25, 0.27, 0.29, 0.30]
    all_rows = {}
    for tag, fname, label in VARIANTS:
        fpath = VD / fname
        if fpath.exists():
            all_rows[tag] = list(csv.DictReader(fpath.open()))
    for bias in biases:
        si = abs(ref_interp(ref_rows, bias))
        row_vals = []
        for tag, fname, label in VARIANTS:
            rows = all_rows.get(tag, [])
            row = next((r for r in rows if abs(f(r, "bias_V") - bias) < 0.006), None)
            vi = abs(f(row, "current_total_A_per_um")) if row else math.nan
            row_vals.append(vi)
        vals_str = "  ".join(
            f"{v:>14.4e}" if not math.isnan(v) else f"{'N/A':>14}" for v in row_vals
        )
        print(f"{bias:>6.2f}  {si:>14.4e}  {vals_str}")


if __name__ == "__main__":
    main()
