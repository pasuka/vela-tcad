#!/usr/bin/env python3
"""Probe: does Sentaurus' SRH generation contain a field/doping-dependent term
that Vela's plain midgap Shockley-Read-Hall model lacks?

Context (mesh-refinement verification result): the residual 0.73x BV deficit is
MESH-CONVERGED (uniform 2x edge refinement moved the -13.2 V current only 0.11%),
so it is NOT a discretization artifact but a genuine, mesh-independent MODEL
difference confined to the depletion region (x ~ 0.78-1.22 um).  Vela's
``RecombinationModel`` is plain midgap SRH with CONSTANT lifetimes and NO field
enhancement:

    R_plain = (n*p - ni_eff^2) / ( taup*(n+ni_eff) + taun*(p+ni_eff) )

Sentaurus 2018 PN defaults frequently add a doping-dependent lifetime
(Scharfetter tau(N)) and/or Hurkx trap-assisted field enhancement
R_srh -> R_plain * (1 + Gamma(E)), both of which BOOST depletion-shoulder
generation -- exactly the missing pre-multiplication seed Step B localized.

This script is fully SELF-CONTAINED on Sentaurus' own 1943-node mesh (no
cross-mesh interpolation).  It reconstructs Vela's *plain* SRH rate from
Sentaurus' OWN exported state (n, p, doping -> ni_eff via old_slotboom BGN, the
same model Vela uses) and compares it to Sentaurus' OWN exported
``srhRecombination``.  Because Vela's exported SRH IS exactly R_plain by
construction, the ratio

    g = R_srh_Sentaurus / R_plain_Vela(Sentaurus state)

is the model-level enhancement factor Vela is missing.  Binning g by |E|,
|net doping|, and x separates the candidates:

  * g ~ 1 across all |E| and doping  -> Sentaurus also uses plain midgap SRH;
    SRH is fully exonerated even in the shoulders (move to the continuity
    e/h-partition candidate).
  * g rising with |E| (flat in doping) -> Hurkx field enhancement: THE missing
    seed; its magnitude in the shoulders quantifies the deficit.
  * g rising with |net doping| (flat in |E|) -> doping-dependent lifetime tau(N).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_flux_forms as fluxforms

sgdiag = fluxforms.sgdiag

CM3_TO_M3 = 1.0e6           # cm^-3 -> m^-3   and   cm^-3 s^-1 -> m^-3 s^-1
VPCM_TO_VPM = 1.0e2         # V cm^-1 -> V m^-1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-dir", type=Path, required=True,
                        help="Sentaurus state dir (has nodes.csv, doping.csv, fields/).")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--material-ni-cm3", type=float, default=1.4638914958767616e10)
    parser.add_argument("--taun", type=float, default=1.0e-5)
    parser.add_argument("--taup", type=float, default=3.0e-6)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--bandgap-narrowing", default="old_slotboom")
    parser.add_argument("--bias", type=float, default=-13.2)
    return parser.parse_args()


def read_scalar_column(path: Path, node_count: int, scale: float) -> list[float]:
    values: dict[int, float] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        value_col = next((n for n in reader.fieldnames or [] if n != "node_id"), None)
        if value_col is None:
            raise RuntimeError(f"{path} has no value column")
        for row in reader:
            values[int(row["node_id"])] = float(row[value_col]) * scale
    return [values[i] for i in range(node_count)]


def read_nodes(path: Path) -> tuple[list[float], list[float]]:
    xs: dict[int, float] = {}
    ys: dict[int, float] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            nid = int(row["id"])
            xs[nid] = float(row["x_um"])
            ys[nid] = float(row["y_um"])
    n = max(xs) + 1
    return [xs[i] for i in range(n)], [ys[i] for i in range(n)]


def read_net_doping_m3(path: Path, node_count: int) -> list[float]:
    donors = [0.0] * node_count
    acceptors = [0.0] * node_count
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            nid = int(row["node_id"])
            donors[nid] = float(row.get("donors_cm3", 0.0) or 0.0) * CM3_TO_M3
            acceptors[nid] = float(row.get("acceptors_cm3", 0.0) or 0.0) * CM3_TO_M3
    return [donors[i] - acceptors[i] for i in range(node_count)]


def plain_srh_rate(n: float, p: float, ni: float, taun: float, taup: float) -> float:
    den = taup * (n + ni) + taun * (p + ni)
    if abs(den) < 1.0e-100:
        return 0.0
    return (n * p - ni * ni) / den


def stat_block(values: list[float]) -> dict[str, float] | None:
    clean = sorted(v for v in values if v is not None and math.isfinite(v))
    if not clean:
        return None
    return {
        "count": len(clean),
        "min": clean[0],
        "p10": clean[max(0, int(0.10 * (len(clean) - 1)))],
        "median": statistics.median(clean),
        "p90": clean[min(len(clean) - 1, int(0.90 * (len(clean) - 1)))],
        "max": clean[-1],
    }


def junction_x_um(xs: list[float], net: list[float]) -> float:
    """Estimate junction position as the |net doping| minimum along the device."""
    best_x = None
    best_abs = None
    for x, nd in zip(xs, net):
        if best_abs is None or abs(nd) < best_abs:
            best_abs = abs(nd)
            best_x = x
    return best_x if best_x is not None else 0.0


def main() -> None:
    args = parse_args()
    sd = args.sentaurus_dir
    fields = sd / "fields"
    nodes_csv = sd / "nodes.csv"
    doping_csv = sd / "doping.csv"

    xs, ys = read_nodes(nodes_csv)
    node_count = len(xs)

    vt = sgdiag.K_B_OVER_Q * args.temperature_k
    ni_m3 = args.material_ni_cm3 * CM3_TO_M3

    net = read_net_doping_m3(doping_csv, node_count)
    ni_eff = fluxforms.effective_ni_from_doping(
        doping_csv, node_count, ni_m3, vt, args.bandgap_narrowing)

    n = read_scalar_column(fields / "eDensity_region0.csv", node_count, CM3_TO_M3)
    p = read_scalar_column(fields / "hDensity_region0.csv", node_count, CM3_TO_M3)
    efield = read_scalar_column(fields / "ElectricField_region0.csv", node_count, VPCM_TO_VPM)
    srh_sent = read_scalar_column(fields / "srhRecombination_region0.csv", node_count, CM3_TO_M3)

    junction = junction_x_um(xs, net)

    # Per-node enhancement factor g = R_srh_Sentaurus / R_plain_Vela.
    rows: list[dict[str, Any]] = []
    for i in range(node_count):
        r_plain = plain_srh_rate(n[i], p[i], ni_eff[i], args.taun, args.taup)
        absE = abs(efield[i])
        g = None
        if abs(r_plain) > 0.0 and math.isfinite(r_plain) and math.isfinite(srh_sent[i]):
            g = srh_sent[i] / r_plain
        rows.append({
            "node": i,
            "x_um": xs[i],
            "y_um": ys[i],
            "dist_to_junction_um": xs[i] - junction,
            "net_doping_m3": net[i],
            "abs_net_doping_m3": abs(net[i]),
            "ni_eff_m3": ni_eff[i],
            "n_m3": n[i],
            "p_m3": p[i],
            "absE_Vpm": absE,
            "R_plain_vela": r_plain,
            "R_srh_sentaurus": srh_sent[i],
            "g_enhancement": g,
        })

    # Only nodes where SRH is net GENERATION (reverse-bias depletion, r_plain<0)
    # and the magnitude is non-negligible carry the seed.  Filter for the
    # depletion-relevant set: r_plain < 0 and |r_plain| above a small floor.
    gen = [r for r in rows
           if r["g_enhancement"] is not None
           and r["R_plain_vela"] < 0.0
           and abs(r["R_plain_vela"]) > 1.0e15]

    def bin_by(rows_in: list[dict[str, Any]], key: str, edges: list[float]) -> list[dict[str, Any]]:
        out = []
        labels = []
        for k in range(len(edges) + 1):
            lo = edges[k - 1] if k > 0 else None
            hi = edges[k] if k < len(edges) else None
            labels.append((lo, hi))
        for lo, hi in labels:
            sel = []
            for r in rows_in:
                v = r[key]
                if (lo is None or v >= lo) and (hi is None or v < hi):
                    sel.append(r["g_enhancement"])
            block = stat_block(sel)
            out.append({"lo": lo, "hi": hi, "g_stats": block})
        return out

    field_edges = [1.0e5, 1.0e6, 3.0e6, 1.0e7, 2.0e7, 3.0e7]
    doping_edges = [1.0e21, 1.0e22, 1.0e23, 5.0e23]

    # Depletion-window restriction (the localized seed region).
    depletion = [r for r in gen if abs(r["dist_to_junction_um"]) <= 0.25]

    summary = {
        "bias_V": args.bias,
        "sentaurus_dir": str(sd),
        "node_count": node_count,
        "junction_x_um": junction,
        "vt_V": vt,
        "ni_material_m3": ni_m3,
        "taun_s": args.taun,
        "taup_s": args.taup,
        "model_note": (
            "g = R_srh_Sentaurus / R_plain_Vela on Sentaurus state. "
            "g~1 everywhere => identical plain SRH (SRH exonerated). "
            "g rising with |E| => Hurkx field enhancement (missing seed). "
            "g rising with |net doping| (flat in E) => doping-dependent tau(N)."
        ),
        "generation_node_count": len(gen),
        "depletion_window_node_count": len(depletion),
        "g_overall_generation": stat_block([r["g_enhancement"] for r in gen]),
        "g_depletion_window": stat_block([r["g_enhancement"] for r in depletion]),
        "g_by_absE_generation": bin_by(gen, "absE_Vpm", field_edges),
        "g_by_absdoping_generation": bin_by(gen, "abs_net_doping_m3", doping_edges),
        "g_by_absE_depletion": bin_by(depletion, "absE_Vpm", field_edges),
        "g_by_absdoping_depletion": bin_by(depletion, "abs_net_doping_m3", doping_edges),
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.out_dir / "srh_field_enhancement_summary.json"
    with summary_path.open("w") as handle:
        json.dump(summary, handle, indent=2)

    # Per-node CSV along the central junction-normal cut (smallest |y|).
    ymin = min(abs(y) for y in ys)
    cut = sorted((r for r in rows if abs(r["y_um"]) <= ymin + 1.0e-9),
                 key=lambda r: r["x_um"])
    cut_path = args.out_dir / "srh_field_enhancement_cut.csv"
    with cut_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in cut:
            writer.writerow(r)

    print(json.dumps({
        "summary": str(summary_path),
        "cut_csv": str(cut_path),
        "junction_x_um": junction,
        "generation_nodes": len(gen),
        "g_overall": summary["g_overall_generation"],
        "g_depletion": summary["g_depletion_window"],
    }, indent=2))


if __name__ == "__main__":
    main()
