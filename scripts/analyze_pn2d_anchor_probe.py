#!/usr/bin/env python3
"""Analyze pn2d anchor-bias contact-edge diagnostic bundle.

Reads pn2d_iv_anchor_probe.csv and pn2d_iv_anchor_probe_contact_edges.csv
from the given vela directory and prints a structured summary of:
  - Total currents at anchor biases (0.27, 0.29, 0.30 V)
  - Sentaurus reference at the same biases (from reference_curves/)
  - Per-contact per-edge drift/diffusion balance at each anchor
  - Dominant electron-branch edges at the anode contact

Usage:
    python scripts/analyze_pn2d_anchor_probe.py [--vela-dir build/pn2d_probe/vela]
"""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

ANCHOR_BIASES = [0.27, 0.29, 0.30]


def _f(row: dict, key: str, default: float = 0.0) -> float:
    v = row.get(key, "")
    if v is None or v == "":
        return default
    try:
        return float(v)
    except (ValueError, TypeError):
        return default


def read_csv(path: Path) -> list[dict]:
    return list(csv.DictReader(path.open(encoding="utf-8")))


def find_row_at_bias(rows: list[dict], bias: float, key: str = "bias_V") -> dict | None:
    best, best_err = None, 1e300
    for row in rows:
        try:
            err = abs(float(row[key]) - bias)
        except (KeyError, ValueError):
            continue
        if err < best_err:
            best, best_err = row, err
    if best is None or best_err > max(1e-9, abs(bias) * 1e-6):
        return None
    return best


def sentaurus_current_at(ref_rows: list[dict], bias: float) -> float:
    """Linearly interpolate Sentaurus reference current (A/um) at given bias."""
    pts: list[tuple[float, float]] = []
    for row in ref_rows:
        try:
            pts.append((float(row["bias_V"]), float(row["current_total"])))
        except (KeyError, ValueError):
            continue
    pts.sort()
    if not pts:
        return math.nan
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--vela-dir",
        type=Path,
        default=REPO / "build" / "pn2d_probe" / "vela",
        help="Directory containing the probe CSV outputs",
    )
    parser.add_argument(
        "--ref-dir",
        type=Path,
        default=REPO / "build" / "pn2d_probe" / "reference_curves",
        help="Directory containing Sentaurus reference curves",
    )
    args = parser.parse_args()

    vela_dir: Path = args.vela_dir
    ref_dir: Path = args.ref_dir

    iv_csv = vela_dir / "pn2d_iv_anchor_probe.csv"
    ce_csv = vela_dir / "pn2d_iv_anchor_probe_contact_edges.csv"

    if not iv_csv.exists():
        print(f"ERROR: IV CSV not found: {iv_csv}")
        return
    if not ce_csv.exists():
        print(f"ERROR: Contact-edge CSV not found: {ce_csv}")
        return

    iv_rows = read_csv(iv_csv)
    ce_rows = read_csv(ce_csv)

    # Load Sentaurus reference
    ref_iv_csv = ref_dir / "pn2d_iv_reference.csv"
    ref_rows = read_csv(ref_iv_csv) if ref_iv_csv.exists() else []

    # ── IV summary ───────────────────────────────────────────────────────────
    print("=" * 72)
    print("pn2d anchor-bias IV summary")
    print("=" * 72)
    print(f"{'Bias (V)':>10}  {'Vela I (A/um)':>18}  {'Sent I (A/um)':>18}  {'orders':>8}  {'handoff':>12}  {'newton':>7}")
    print("-" * 80)
    for bias in ANCHOR_BIASES:
        row = find_row_at_bias(iv_rows, bias)
        if row is None:
            print(f"{bias:>10.3f}  {'N/A':>18}")
            continue
        v_i = _f(row, "current_total_A_per_um")
        hs = row.get("handoff_stage", "?")
        ni = row.get("newton_iterations", "?")
        s_i = sentaurus_current_at(ref_rows, bias) if ref_rows else math.nan
        ord_ = orders(s_i, v_i) if ref_rows and not math.isnan(s_i) else math.nan
        sent_str = f"{s_i:.6e}" if not math.isnan(s_i) else "N/A"
        ord_str = f"{ord_:.4f}" if not math.isnan(ord_) else "N/A"
        print(f"{bias:>10.3f}  {v_i:>18.6e}  {sent_str:>18}  {ord_str:>8}  {hs:>12}  {ni:>7}")

    # ── Contact-edge analysis per anchor ─────────────────────────────────────
    print()
    print("=" * 72)
    print("Contact-edge diagnostics at anchor biases")
    print("=" * 72)

    # Group CE rows by (bias, contact)
    by_bias_contact: dict[tuple[float, str], list[dict]] = defaultdict(list)
    for row in ce_rows:
        try:
            b = round(float(row["bias_V"]), 6)
        except (KeyError, ValueError):
            continue
        contact = row.get("current_contact", "?")
        by_bias_contact[(b, contact)].append(row)

    for bias in ANCHOR_BIASES:
        # Find nearest bias key in ce data
        best_key = None
        best_err = 1e300
        for (b, c) in by_bias_contact:
            err = abs(b - bias)
            if err < best_err:
                best_err = err
                best_key_bias = b
        if best_err > 0.005:
            print(f"\nNo contact-edge data near V={bias:.3f} V")
            continue

        print(f"\n── V = {best_key_bias:.4f} V ──")
        contacts_at_bias = {c for (b, c) in by_bias_contact if b == best_key_bias}
        for contact in sorted(contacts_at_bias):
            edges = by_bias_contact[(best_key_bias, contact)]
            n_edges = len(edges)
            total_sum = sum(_f(e, "current_total") for e in edges)
            e_drift_sum = sum(_f(e, "current_electron_drift") for e in edges)
            e_diff_sum = sum(_f(e, "current_electron_diffusion") for e in edges)
            h_drift_sum = sum(_f(e, "current_hole_drift") for e in edges)
            h_diff_sum = sum(_f(e, "current_hole_diffusion") for e in edges)
            qf_edges = [e for e in edges if e.get("electron_branch") == "quasi_fermi"]
            dens_edges = [e for e in edges if e.get("electron_branch") == "density"]

            print(f"  Contact: {contact}  ({n_edges} edges)")
            print(f"    Total current sum:          {total_sum:+.6e} A/m")
            print(f"    Electron drift sum:         {e_drift_sum:+.6e} A/m")
            print(f"    Electron diffusion sum:     {e_diff_sum:+.6e} A/m")
            print(f"    Hole drift sum:             {h_drift_sum:+.6e} A/m")
            print(f"    Hole diffusion sum:         {h_diff_sum:+.6e} A/m")
            if n_edges > 0:
                frac_drift = abs(e_drift_sum) / (abs(e_drift_sum) + abs(e_diff_sum) + 1e-300)
                print(f"    Electron drift fraction:    {frac_drift:.4f}")
                print(f"    QF-branch edges: {len(qf_edges)}/{n_edges}  density-branch: {len(dens_edges)}/{n_edges}")

            # Top 5 edges by absolute total current
            sorted_edges = sorted(edges, key=lambda e: abs(_f(e, "current_total")), reverse=True)
            print(f"    Top 5 edges by |current_total|:")
            for e in sorted_edges[:5]:
                eid = e.get("edge_id", "?")
                n0, n1 = e.get("node0", "?"), e.get("node1", "?")
                i_tot = _f(e, "current_total")
                i_ed = _f(e, "current_electron_drift")
                i_edi = _f(e, "current_electron_diffusion")
                eb = e.get("electron_branch", "?")
                sign = e.get("outward_sign", "?")
                print(f"      edge {eid:>6} [{n0},{n1}] sign={sign:>3}  I_tot={i_tot:+.4e}"
                      f"  e_drift={i_ed:+.4e}  e_diff={i_edi:+.4e}  e_branch={eb}")

    # ── Ratio I(0.29)/I(0.30) ────────────────────────────────────────────────
    print()
    print("=" * 72)
    row_029 = find_row_at_bias(iv_rows, 0.29)
    row_030 = find_row_at_bias(iv_rows, 0.30)
    if row_029 and row_030:
        i029 = abs(_f(row_029, "current_total_A_per_um"))
        i030 = abs(_f(row_030, "current_total_A_per_um"))
        ratio = i029 / i030 if i030 != 0 else math.nan
        s_i029 = sentaurus_current_at(ref_rows, 0.29) if ref_rows else math.nan
        s_i030 = sentaurus_current_at(ref_rows, 0.30) if ref_rows else math.nan
        s_ratio = abs(s_i029) / abs(s_i030) if (ref_rows and s_i030 != 0 and not math.isnan(s_i030)) else math.nan
        print(f"I(0.29)/I(0.30) slope ratio:")
        print(f"  Vela:      {ratio:.6f}")
        if not math.isnan(s_ratio):
            print(f"  Sentaurus: {s_ratio:.6f}")
            print(f"  Delta:     {abs(ratio - s_ratio):.6f}  (target: < 0.06)")
    print("=" * 72)


if __name__ == "__main__":
    main()
