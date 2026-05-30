#!/usr/bin/env python3
"""Stage 0 decision table: compute the OFFICIAL IV gate metric for the four levers.

Gate metric replicates scripts/compare_reference_curves.orders_of_magnitude:
  orders = max over aligned bias points in [bias_min, bias_max] of
           | log10( |cand * scale| / |ref| ) |
with scale = -1.0, cand = current_total_A_per_um, ref = Sentaurus current_total,
bias window [0.2, 0.3]. The regression test asserts this <= 0.50.

Two flavors are reported:
  - gate_2pt : official deck uses vela_step 0.1, so only biases 0.20 and 0.30
               fall in the window. This is the exact number the test gates on.
  - window_full : max over ALL 0.20-0.30 points (stricter, informational).

Also reports local slope ratio I(0.29)/I(0.30) (Sentaurus ~0.632) and handoff.
"""
import csv
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import compare_reference_curves as crc  # noqa: E402

HERE = REPO / "build" / "pn2d_probe" / "vela"
REF = REPO / "build" / "pn2d_probe" / "reference_curves" / "pn2d_iv_reference.csv"
SCALE = -1.0
BIAS_MIN, BIAS_MAX = 0.2, 0.3

VARIANTS = [
    ("s0_baseline",     "tau=1e-7,  ni=1.0e10  (current deck)"),
    ("s1_tau1e6",       "tau=1e-6,  ni=1.0e10"),
    ("s2_ni145",        "tau=1e-7,  ni=1.45e10"),
    ("s3_tau1e6_ni145", "tau=1e-6,  ni=1.45e10"),
]


def fval(row, k):
    try:
        return float(row.get(k, "") or "nan")
    except (TypeError, ValueError):
        return float("nan")


def load_ref_rows():
    return list(csv.DictReader(REF.open()))


def load_ref_pairs():
    rows = load_ref_rows()
    return sorted((float(r["bias_V"]), float(r["current_total"])) for r in rows)


def ref_interp(pts, bias):
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


def load_variant(tag):
    fpath = HERE / f"pn2d_iv_{tag}.csv"
    rows = list(csv.DictReader(fpath.open()))
    data = {}
    handoff_ok = True
    for r in rows:
        if fval(r, "converged") != 1:
            continue
        b = round(fval(r, "bias_V"), 4)
        data[b] = {"i": fval(r, "current_total_A_per_um")}
        if r.get("handoff_stage", "") != "newton" or fval(r, "newton_iterations") <= 0:
            handoff_ok = False
    return rows, data, handoff_ok


def official_gate(ref_rows, cand_rows):
    """Exact replica of the regression gate via compare_reference_curves."""
    rep = crc.compare_series(
        "iv", ref_rows, cand_rows,
        candidate_scale=SCALE,
        bias_min=BIAS_MIN, bias_max=BIAS_MAX,
        candidate_column="current_total_A_per_um",
    )
    return rep.get("orders_of_magnitude"), rep.get("points_compared"), rep.get("trend_match")


def main():
    ref_rows = load_ref_rows()
    pts = load_ref_pairs()
    si29 = abs(ref_interp(pts, 0.29))
    si30 = abs(ref_interp(pts, 0.30))
    sent_slope = si29 / si30 if si30 else math.nan

    print("=" * 104)
    print("STAGE 0 DECISION TABLE — OFFICIAL IV gate metric (test asserts gate <= 0.50)")
    print("=" * 104)
    print(f"Sentaurus I(0.30)={si30:.4e} A/um   slope I(0.29)/I(0.30)={sent_slope:.4f}   "
          f"gate aligns at reference bias points in [{BIAS_MIN},{BIAS_MAX}]")
    print()
    hdr = (f"{'variant':<16}{'desc':<38}{'GATE':>8}{'pts':>4}{'trend':>6}"
           f"{'I(0.30)':>12}{'gap_x':>8}{'slopeR':>8}{'slopeD':>8}{'ho':>4}")
    print(hdr)
    print("-" * 104)
    rows_out = []
    for tag, desc in VARIANTS:
        cand_rows, data, handoff_ok = load_variant(tag)
        gate, npts, trend = official_gate(ref_rows, cand_rows)
        i30 = data.get(0.30, {}).get("i", math.nan)
        i29 = data.get(0.29, {}).get("i", math.nan)
        gapx = si30 / abs(i30) if i30 and math.isfinite(i30) else math.nan
        slopeR = abs(i29) / abs(i30) if i30 and math.isfinite(i29) and math.isfinite(i30) else math.nan
        slopeD = abs(slopeR - sent_slope)
        ho = "OK" if handoff_ok else "FAIL"
        tm = "yes" if trend else "no"
        print(f"{tag:<16}{desc:<38}{gate:>8.4f}{npts:>4}{tm:>6}"
              f"{abs(i30):>12.4e}{gapx:>8.2f}{slopeR:>8.4f}{slopeD:>8.4f}{ho:>4}")
        rows_out.append((tag, gate, gapx, slopeR, slopeD, handoff_ok, trend))

    print()
    print("Gate decision (regression test: IV orders_of_magnitude <= 0.50):")
    base_gate = next(r[1] for r in rows_out if r[0] == "s0_baseline")
    base_slope = next(r[3] for r in rows_out if r[0] == "s0_baseline")
    for tag, gate, gapx, slopeR, slopeD, ho, trend in rows_out:
        pass_gate = "PASS" if gate <= 0.50 else "FAIL(>0.50)"
        dgate = gate - base_gate
        dslope = slopeR - base_slope
        print(f"  {tag:<16} GATE={gate:.4f} ({pass_gate:>11})  "
              f"Δgate_vs_base={dgate:+.4f}  Δslope_vs_base={dslope:+.4f}  "
              f"trend={'ok' if trend else 'NO'}  handoff={'OK' if ho else 'FAIL'}")

    # Coarse step-0.1 gate replica (matches the regression deck exactly)
    print()
    print("=" * 104)
    print("AUTHORITATIVE gate at official deck step (vela_step=0.1) — matches the regression test")
    print("=" * 104)
    gate_variants = [
        ("s0g_baseline",     "tau=1e-7,  ni=1.0e10  (current deck)"),
        ("s1g_tau1e6",       "tau=1e-6,  ni=1.0e10"),
        ("s2g_ni145",        "tau=1e-7,  ni=1.45e10"),
        ("s3g_tau1e6_ni145", "tau=1e-6,  ni=1.45e10"),
    ]
    print(f"{'variant':<18}{'desc':<38}{'GATE':>8}{'pts':>4}{'decision':>14}")
    print("-" * 104)
    for tag, desc in gate_variants:
        fpath = HERE / f"pn2d_iv_{tag}.csv"
        if not fpath.exists():
            print(f"{tag:<18}{desc:<38}{'(missing)':>8}")
            continue
        cand_rows = list(csv.DictReader(fpath.open()))
        gate, npts, trend = official_gate(ref_rows, cand_rows)
        decision = "PASS" if gate <= 0.50 else "FAIL(>0.50)"
        print(f"{tag:<18}{desc:<38}{gate:>8.4f}{npts:>4}{decision:>14}")


if __name__ == "__main__":
    main()
