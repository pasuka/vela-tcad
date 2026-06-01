#!/usr/bin/env python3
"""Quick anode vs cathode contact-edge comparison for pn2d anchor biases."""
import csv
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
VD = REPO / "build" / "pn2d_probe" / "vela"
RD = REPO / "build" / "pn2d_probe" / "reference_curves"
ANCHORS = [0.27, 0.29, 0.30]


def f(row: dict, k: str, d: float = 0.0) -> float:
    try:
        return float(row.get(k, "") or d)
    except (TypeError, ValueError):
        return d


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
    eps = 1e-300
    return abs(math.log10(max(abs(cand), eps) / max(abs(ref), eps)))


def analyze(ce_file: Path, label: str, ref_rows: list) -> None:
    if not ce_file.exists():
        print(f"{label}: file not found: {ce_file}")
        return
    ce = list(csv.DictReader(ce_file.open()))
    print(f"\n{'='*66}")
    print(f"  {label} contact edges  ({ce_file.name})")
    print(f"{'='*66}")
    print(f"{'Bias':>6}  {'N':>2}  {'I_tot (A/m)':>14}  {'e_drift':>14}  {'e_diff':>14}  {'drift_frac':>10}  {'branch':>11}")
    for bias in ANCHORS:
        rows = [r for r in ce if abs(f(r, "bias_V") - bias) < 0.006]
        if not rows:
            print(f"{bias:>6.3f}  -- no data")
            continue
        ct = sum(f(r, "current_total") for r in rows)
        ed = sum(f(r, "current_electron_drift") for r in rows)
        edi = sum(f(r, "current_electron_diffusion") for r in rows)
        hd = sum(f(r, "current_hole_drift") for r in rows)
        hdi = sum(f(r, "current_hole_diffusion") for r in rows)
        frac = abs(ed) / (abs(ed) + abs(edi) + 1e-300)
        branch = rows[0].get("electron_branch", "?")
        print(f"{bias:>6.3f}  {len(rows):>2}  {ct:>14.4e}  {ed:>14.4e}  {edi:>14.4e}  {frac:>10.4f}  {branch:>11}")
        print(f"{'':>6}      {'h_drift':>14}: {hd:>14.4e}  h_diff: {hdi:>14.4e}")
    # top edges at 0.30 V
    rows_30 = [r for r in ce if abs(f(r, "bias_V") - 0.30) < 0.006]
    if rows_30:
        rows_30_sorted = sorted(rows_30, key=lambda r: abs(f(r, "current_total")), reverse=True)
        print(f"\n  Top edges at 0.30 V:")
        for r in rows_30_sorted[:5]:
            eid = r.get("edge_id")
            n0, n1 = r.get("node0"), r.get("node1")
            sign = r.get("outward_sign")
            ct = f(r, "current_total")
            ed = f(r, "current_electron_drift")
            edi = f(r, "current_electron_diffusion")
            bu = f(r, "bernoulli_u")
            print(f"    edge {eid:>6} [{n0},{n1}] sign={sign:>3}  u={bu:>8.3f}  I={ct:+.3e}  ed={ed:+.3e}  edi={edi:+.3e}")


def main() -> None:
    ref_rows: list = []
    ref_iv = RD / "pn2d_iv_reference.csv"
    if ref_iv.exists():
        ref_rows = list(csv.DictReader(ref_iv.open()))

    # IV summary
    iv_csv = VD / "pn2d_iv_anchor_probe.csv"
    if iv_csv.exists():
        iv = list(csv.DictReader(iv_csv.open()))
        print("IV anchor summary:")
        print(f"{'Bias':>6}  {'Vela I (A/um)':>16}  {'Sent I (A/um)':>16}  {'orders':>8}  {'handoff':>10}")
        for bias in ANCHORS:
            row = next((r for r in iv if abs(f(r, "bias_V") - bias) < 0.006), None)
            if row is None:
                print(f"{bias:>6.3f}  no data")
                continue
            vi = f(row, "current_total_A_per_um")
            ri = ref_interp(ref_rows, bias) if ref_rows else math.nan
            ord_ = orders(abs(ri), abs(vi)) if not math.isnan(ri) else math.nan
            hs = row.get("handoff_stage", "?")
            print(f"{bias:>6.3f}  {vi:>16.6e}  {abs(ri):>16.6e}  {ord_:>8.4f}  {hs:>10}")
        row29 = next((r for r in iv if abs(f(r, "bias_V") - 0.29) < 0.006), None)
        row30 = next((r for r in iv if abs(f(r, "bias_V") - 0.30) < 0.006), None)
        if row29 and row30:
            i29 = abs(f(row29, "current_total_A_per_um"))
            i30 = abs(f(row30, "current_total_A_per_um"))
            r29 = abs(ref_interp(ref_rows, 0.29)) if ref_rows else math.nan
            r30 = abs(ref_interp(ref_rows, 0.30)) if ref_rows else math.nan
            vr = i29 / i30 if i30 else math.nan
            sr = r29 / r30 if r30 else math.nan
            print(f"\nI(0.29)/I(0.30)  Vela={vr:.6f}  Sent={sr:.6f}  delta={abs(vr-sr):.6f}")

    analyze(VD / "pn2d_iv_anchor_probe_contact_edges.csv", "CATHODE", ref_rows)
    analyze(VD / "pn2d_iv_anchor_probe_anode_edges.csv", "ANODE", ref_rows)


if __name__ == "__main__":
    main()
