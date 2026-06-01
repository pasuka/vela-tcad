"""Probe pn2d IV per-bias contact-side drift/diffusion decomposition.

Reads existing default/iv_bgn_none/iv_recomb_none Vela IV CSVs and the
Sentaurus reference CSV; emits a per-bias table of:

  bias_V, Vela I_total, ref I_total, ratio,
  I_e_total, I_e_drift, I_e_diff,
  I_h_total, I_h_drift, I_h_diff,
  (drift+diff cancellation magnitude on electron side),
  (apparent doubling of |drift| relative to |total|).

The goal is to localize the IV high-forward-bias slope deviation by
following how the contact-side carrier-current composition evolves with
bias for each model toggle.

Run with the vela-tcad python (no extra deps; uses csv + bisect).
"""

from __future__ import annotations

import bisect
import csv
from pathlib import Path
from typing import Dict, List, Tuple

BUILD_DIR = Path("build/pn2d_tdr_tie_probe")
VELA_DIR = BUILD_DIR / "vela"
REF_CSV = BUILD_DIR / "reference_curves" / "pn2d_iv_reference.csv"


def read_vela_iv(path: Path) -> List[Dict[str, float]]:
    rows: List[Dict[str, float]] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for raw in reader:
            try:
                bias = float(raw["bias_V"])
            except (KeyError, ValueError):
                continue
            if raw.get("handoff_stage", "") != "newton":
                continue
            rows.append(
                {
                    "bias_V": bias,
                    "I_total": float(raw["current_total_A_per_um"]),
                    "I_e_total": float(raw["current_electron_A_per_um"]),
                    "I_e_drift": float(raw["current_electron_drift_A_per_um"]),
                    "I_e_diff": float(raw["current_electron_diffusion_A_per_um"]),
                    "I_h_total": float(raw["current_hole_A_per_um"]),
                    "I_h_drift": float(raw["current_hole_drift_A_per_um"]),
                    "I_h_diff": float(raw["current_hole_diffusion_A_per_um"]),
                }
            )
    rows.sort(key=lambda r: r["bias_V"])
    return rows


def read_reference_iv(path: Path) -> Tuple[List[float], List[float]]:
    biases: List[float] = []
    currents: List[float] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for raw in reader:
            try:
                b = float(raw["bias_V"])
                i = float(raw["current_total"])
            except (KeyError, ValueError):
                continue
            biases.append(b)
            currents.append(i)
    paired = sorted(zip(biases, currents))
    biases = [b for b, _ in paired]
    currents = [i for _, i in paired]
    return biases, currents


def interp_ref(biases: List[float], currents: List[float], target_bias: float) -> float:
    idx = bisect.bisect_left(biases, target_bias)
    if idx <= 0:
        return currents[0]
    if idx >= len(biases):
        return currents[-1]
    b0, b1 = biases[idx - 1], biases[idx]
    i0, i1 = currents[idx - 1], currents[idx]
    if b1 == b0:
        return i0
    t = (target_bias - b0) / (b1 - b0)
    return i0 + t * (i1 - i0)


def fmt(x: float) -> str:
    if x == 0.0:
        return "0"
    return f"{x: .3e}"


def main() -> None:
    ref_b, ref_i = read_reference_iv(REF_CSV)
    cases = [
        ("default", VELA_DIR / "pn2d_iv_default.csv"),
        ("iv_bgn_none", VELA_DIR / "pn2d_iv_iv_bgn_none.csv"),
        ("iv_recomb_none", VELA_DIR / "pn2d_iv_iv_recomb_none.csv"),
    ]

    for name, path in cases:
        if not path.exists():
            print(f"# SKIP {name}: missing {path}")
            continue
        rows = read_vela_iv(path)
        print(f"\n== {name} ({path.name}) ==")
        print(
            f"{'bias':>6} | {'I_total':>11} {'ref_I':>11} {'ratio':>7} | "
            f"{'I_e':>11} {'I_e_drift':>11} {'I_e_diff':>11} {'|cancel/total|':>14} | "
            f"{'I_h':>11} {'I_h_drift':>11} {'I_h_diff':>11}"
        )
        for r in rows:
            b = r["bias_V"]
            if b <= 1e-9:
                continue
            ref = interp_ref(ref_b, ref_i, b)
            ratio = (r["I_total"] / ref) if ref != 0 else float("nan")
            # Cancellation magnitude on electron side: how much larger is the
            # raw drift than the net electron current (signal-to-noise proxy)?
            cancel = (
                abs(r["I_e_drift"]) / abs(r["I_e_total"])
                if r["I_e_total"] != 0
                else float("nan")
            )
            print(
                f"{b:6.3f} | {fmt(r['I_total'])} {fmt(ref)} {ratio:7.3f} | "
                f"{fmt(r['I_e_total'])} {fmt(r['I_e_drift'])} {fmt(r['I_e_diff'])} "
                f"{cancel:14.2f} | "
                f"{fmt(r['I_h_total'])} {fmt(r['I_h_drift'])} {fmt(r['I_h_diff'])}"
            )


if __name__ == "__main__":
    main()
