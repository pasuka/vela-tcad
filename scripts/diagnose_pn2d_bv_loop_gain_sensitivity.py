#!/usr/bin/env python3
"""Step B: loop-gain / sensitivity (dG/dstate) decomposition of the PN2D BV
converged-state current deficit, instead of yet another converged-value compare.

Background (Steps 1-6, repo memory): at -13.2 V the Vela terminal current is
~0.73x Sentaurus 2018. Steps 4-5 localized it to a carrier-density deficit
(electron ~0.80x, hole ~0.71x) that is purely an absolute quasi-Fermi-to-
potential *level* offset (d(psi-phin) ~ -6 mV, d(phip-psi) ~ -9 mV) with psi,
ni/BGN, mobility, the SG flux form, and every per-edge gradient field matching
Sentaurus to <2%. Step 6 established 0.73x is a TRUE r~=0 fixed point of Vela's
discretized equations -- a sub-percent state difference amplified by a strong
nonlinearity -- and asked (item B) for a loop-gain / sensitivity comparison
rather than another converged-value comparison.

This script separates the two candidate amplifiers of that small state offset:

  (1) BOLTZMANN forward sensitivity of avalanche generation to the quasi-Fermi
      LEVEL.  Because G = alpha(F) * flux and flux ~ carrier density
      ~ exp((psi - phin)/Vt), the open-loop forward gain is
          d ln G / d(psi - phin) = 1/Vt  (~38.7 /V at 300 K).
      So a level offset dL converts to a generation ratio exp(dL/Vt).  We verify
      the measured electron/hole generation (density) ratios equal
      exp(d level / Vt) and report the per-mV gain (~3.9%/mV).

  (2) AVALANCHE-MULTIPLICATION loop gain.  The ionization integral
          I_ion = integral(alpha dx)  along the junction-normal streamline sets
      the multiplication M = 1/(1 - I_ion), whose sensitivity is dM/dI = M^2.
      We integrate alpha(F) along central mesh rows for the Vela and Sentaurus
      states under BOTH quasi-Fermi-gradient driving (what Vela's deck uses) and
      electric-field driving (the low-density fallback Sentaurus uses), split
      into the high-field active band (E >= threshold) and the low-field tail,
      and compare M_V vs M_S and M_efield vs M_qf.

If (2) shows M_V ~ M_S (active fields already match, tail negligible) then the
multiplication loop gain is NOT the amplifier, and the 0.73x is a direct
Boltzmann-sensitive generation deficit fed by the ~7 mV quasi-Fermi level offset
-- i.e. the lever is absolute (psi - phi) alignment in the depletion region, not
any avalanche/source-geometry knob.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_flux_forms as fluxforms

sgdiag = fluxforms.sgdiag

# Van Overstraeten silicon defaults (Sentaurus 2018); Vela matches these.
ELECTRON_ALPHA = (7.03e7, 7.03e7, 1.231e8, 1.231e8)
HOLE_ALPHA = (1.582e8, 6.71e7, 2.036e8, 1.693e8)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--noimpact-vtk-root", type=Path, default=None)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--junction-x-um", type=float, default=1.140625)
    parser.add_argument("--field-threshold", type=float, default=1.0e7)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--material-ni-m3", type=float, default=1.4638914958767616e16)
    parser.add_argument("--bandgap-narrowing", choices=["none", "old_slotboom"], default="old_slotboom")
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def stat_block(values: list[float]) -> dict[str, float] | None:
    clean = [v for v in values if v is not None and math.isfinite(v)]
    if not clean:
        return None
    clean.sort()
    return {
        "count": len(clean),
        "median": statistics.median(clean),
        "mean": statistics.fmean(clean),
        "p10": clean[int(0.10 * (len(clean) - 1))],
        "p90": clean[int(0.90 * (len(clean) - 1))],
        "min": clean[0],
        "max": clean[-1],
    }


# ---------------------------------------------------------------------------
# Metric 1: Boltzmann forward sensitivity d ln G / d level = 1/Vt
# ---------------------------------------------------------------------------
def boltzmann_sensitivity(
    vela: dict[str, list[float]],
    sent: dict[str, list[float]],
    support: list[int],
    vt: float,
) -> dict[str, Any]:
    rows: list[dict[str, float]] = []
    for nid in support:
        level_n_v = vela["psi"][nid] - vela["phin"][nid]
        level_n_s = sent["psi"][nid] - sent["phin"][nid]
        level_p_v = vela["phip"][nid] - vela["psi"][nid]
        level_p_s = sent["phip"][nid] - sent["psi"][nid]
        d_level_n = level_n_v - level_n_s
        d_level_p = level_p_v - level_p_s
        n_v, n_s = vela["n"][nid], sent["n"][nid]
        p_v, p_s = vela["p"][nid], sent["p"][nid]
        if n_s <= 0.0 or p_s <= 0.0:
            continue
        rows.append({
            "d_level_n_mV": d_level_n * 1.0e3,
            "d_level_p_mV": d_level_p * 1.0e3,
            "n_ratio_measured": n_v / n_s,
            "p_ratio_measured": p_v / p_s,
            "n_ratio_boltzmann": math.exp(d_level_n / vt),
            "p_ratio_boltzmann": math.exp(d_level_p / vt),
        })

    def col(key: str) -> list[float]:
        return [r[key] for r in rows]

    # residual between the measured density ratio and the pure Boltzmann
    # prediction exp(d level / Vt): ~0 confirms G ratio = exp(d level / Vt).
    n_closure = [r["n_ratio_measured"] / r["n_ratio_boltzmann"] for r in rows]
    p_closure = [r["p_ratio_measured"] / r["p_ratio_boltzmann"] for r in rows]
    return {
        "node_count": len(rows),
        "forward_gain_per_volt": 1.0 / vt,
        "forward_gain_percent_per_mV": (1.0 / vt) * 1.0e-3 * 100.0,
        "d_level_n_mV": stat_block(col("d_level_n_mV")),
        "d_level_p_mV": stat_block(col("d_level_p_mV")),
        "electron_generation_ratio_measured": stat_block(col("n_ratio_measured")),
        "hole_generation_ratio_measured": stat_block(col("p_ratio_measured")),
        "electron_generation_ratio_boltzmann_pred": stat_block(col("n_ratio_boltzmann")),
        "hole_generation_ratio_boltzmann_pred": stat_block(col("p_ratio_boltzmann")),
        "electron_measured_over_boltzmann": stat_block(n_closure),
        "hole_measured_over_boltzmann": stat_block(p_closure),
    }


# ---------------------------------------------------------------------------
# Metric 2: ionization integral I_ion = integral(alpha dx) and M = 1/(1-I_ion)
# along central junction-normal mesh rows, qf vs efield driving.
# ---------------------------------------------------------------------------
def discover_rows(nodes: dict[int, dict[str, float]]) -> list[list[int]]:
    by_y: dict[float, list[int]] = defaultdict(list)
    for nid, node in nodes.items():
        by_y[round(node["y_um"], 6)].append(nid)
    max_count = max(len(ids) for ids in by_y.values())
    rows: list[list[int]] = []
    for _, ids in sorted(by_y.items()):
        if len(ids) >= 0.8 * max_count:
            rows.append(sorted(ids, key=lambda n: nodes[n]["x_um"]))
    return rows


def integrate_row(
    row: list[int],
    nodes: dict[int, dict[str, float]],
    state: dict[str, list[float]],
    field_threshold: float,
) -> dict[str, float]:
    acc = {
        "n_qf": 0.0, "n_ef": 0.0, "p_qf": 0.0, "p_ef": 0.0,
        "n_qf_active": 0.0, "n_qf_tail": 0.0,
        "n_ef_active": 0.0, "n_ef_tail": 0.0,
        "length_m": 0.0,
    }
    for a, b in zip(row[:-1], row[1:]):
        dx = (nodes[b]["x_um"] - nodes[a]["x_um"]) * 1.0e-6
        dy = (nodes[b]["y_um"] - nodes[a]["y_um"]) * 1.0e-6
        h = math.hypot(dx, dy)
        if h <= 1.0e-30:
            continue
        e_field = sgdiag.field_between(state["psi"], a, b, h)
        e_qf_n = sgdiag.field_between(state["phin"], a, b, h)
        e_qf_p = sgdiag.field_between(state["phip"], a, b, h)
        an_qf = sgdiag.van_overstraeten_alpha(e_qf_n, *ELECTRON_ALPHA)
        an_ef = sgdiag.van_overstraeten_alpha(e_field, *ELECTRON_ALPHA)
        ap_qf = sgdiag.van_overstraeten_alpha(e_qf_p, *HOLE_ALPHA)
        ap_ef = sgdiag.van_overstraeten_alpha(e_field, *HOLE_ALPHA)
        acc["n_qf"] += an_qf * h
        acc["n_ef"] += an_ef * h
        acc["p_qf"] += ap_qf * h
        acc["p_ef"] += ap_ef * h
        acc["length_m"] += h
        if e_field >= field_threshold:
            acc["n_qf_active"] += an_qf * h
            acc["n_ef_active"] += an_ef * h
        else:
            acc["n_qf_tail"] += an_qf * h
            acc["n_ef_tail"] += an_ef * h
    return acc


def multiplication(i_ion: float) -> float:
    return 1.0 / (1.0 - min(i_ion, 0.999999))


def ionization_integral(
    nodes: dict[int, dict[str, float]],
    rows: list[list[int]],
    state: dict[str, list[float]],
    field_threshold: float,
) -> dict[str, Any]:
    per_row = [integrate_row(row, nodes, state, field_threshold) for row in rows]

    def med(key: str) -> float:
        return statistics.median([r[key] for r in per_row])

    i_n_qf = med("n_qf")
    i_n_ef = med("n_ef")
    i_p_qf = med("p_qf")
    i_p_ef = med("p_ef")
    return {
        "row_count": len(rows),
        "I_ion_electron_qf": i_n_qf,
        "I_ion_electron_efield": i_n_ef,
        "I_ion_hole_qf": i_p_qf,
        "I_ion_hole_efield": i_p_ef,
        "I_ion_electron_qf_active": med("n_qf_active"),
        "I_ion_electron_qf_tail": med("n_qf_tail"),
        "I_ion_electron_efield_active": med("n_ef_active"),
        "I_ion_electron_efield_tail": med("n_ef_tail"),
        "M_electron_qf": multiplication(i_n_qf),
        "M_electron_efield": multiplication(i_n_ef),
        "M_hole_qf": multiplication(i_p_qf),
        "M_hole_efield": multiplication(i_p_ef),
        "dM_dI_electron_qf_sensitivity": multiplication(i_n_qf) ** 2,
        "tail_fraction_of_efield_integral_electron":
            (med("n_ef_tail") / i_n_ef) if i_n_ef > 0.0 else None,
    }


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    nodes, _, _ = sgdiag.read_mesh(args.mesh)
    node_count = len(nodes)
    vt = sgdiag.K_B_OVER_Q * args.temperature_k

    vela = fluxforms.load_vela_state(
        fluxforms.discover_vela_vtks(args.vela_vtk_root)[fluxforms.bias_key(args.bias)],
        node_count)
    sent = fluxforms.load_sentaurus_state(args.sentaurus_dir, node_count)
    noimpact = None
    if args.noimpact_vtk_root is not None:
        noimpact = fluxforms.load_vela_state(
            fluxforms.discover_vela_vtks(args.noimpact_vtk_root)[fluxforms.bias_key(args.bias)],
            node_count)

    support = [
        int(row["node_id"]) for row in read_csv_rows(args.support_csv)
        if row.get("support_class") not in (None, "", "inactive")
    ]

    rows = discover_rows(nodes)

    summary: dict[str, Any] = {
        "bias_V": args.bias,
        "vt_V": vt,
        "field_threshold_V_per_m": args.field_threshold,
        "support_node_count": len(support),
        "metric1_boltzmann_forward_sensitivity":
            boltzmann_sensitivity(vela, sent, support, vt),
        "metric2_ionization_integral": {
            "vela": ionization_integral(nodes, rows, vela, args.field_threshold),
            "sentaurus": ionization_integral(nodes, rows, sent, args.field_threshold),
        },
    }
    if noimpact is not None:
        summary["metric2_ionization_integral"]["vela_no_impact"] = \
            ionization_integral(nodes, rows, noimpact, args.field_threshold)

    # Cross-engine multiplication ratio vs the measured terminal-current ratio.
    iv = summary["metric2_ionization_integral"]["vela"]
    isd = summary["metric2_ionization_integral"]["sentaurus"]
    summary["multiplication_cross_check"] = {
        "M_electron_qf_vela_over_sentaurus":
            iv["M_electron_qf"] / isd["M_electron_qf"] if isd["M_electron_qf"] else None,
        "M_electron_efield_over_qf_vela":
            iv["M_electron_efield"] / iv["M_electron_qf"] if iv["M_electron_qf"] else None,
        "note": "If both ~1.0, avalanche multiplication is NOT the amplifier; "
                "the deficit is direct Boltzmann-sensitive generation from the "
                "quasi-Fermi level offset (metric 1).",
    }

    out_json = args.out_dir / "loop_gain_sensitivity_summary.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"\nWrote {out_json}")


if __name__ == "__main__":
    main()
