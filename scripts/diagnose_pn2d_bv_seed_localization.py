#!/usr/bin/env python3
"""Step B-followup: localize the ~7 mV depletion-region quasi-Fermi LEVEL offset
(the primary pre-multiplication carrier "seed") that Step B showed is
Boltzmann-amplified into the 0.73x BV current deficit.

Step B established: avalanche multiplication matches Sentaurus (M_V/M_S = 1.007),
per-edge alpha, SG flux, fields, mobility, ni/BGN all match; the 0.73x is the
matched multiplication faithfully carrying a ~0.79x primary carrier seed,
expressed as a d(psi-phin) ~ -6 mV / d(phip-psi) ~ -9 mV depletion-region level
offset with forward gain 1/Vt = 38.7 /V.

This script separates the two candidate origins of that level offset using the
exported nodal source fields that BOTH engines write:

  * SRH net rate (Vela ``SRHRecombination`` vs Sentaurus ``srhRecombination``).
    In deep reverse depletion np << ni^2 so the SRH GENERATION rate saturates at
    ~ -ni/(tau_n+tau_p), i.e. it is *density-insensitive*.  If R_srh,V ~ R_srh,S
    the primary thermal seed magnitude matches and the lever is transport/BC, not
    generation strength; if R_srh,V ~ 0.79x the seed itself is short (tau/ni/SRH
    formulation).
  * Avalanche generation (Vela ``AvalancheGeneration`` vs Sentaurus
    ``ImpactIonization``) -- the *density-sensitive* channel, expected ~0.76x.

It also traces the quasi-Fermi level offset along the junction-normal cut from
the anode contact to the cathode contact.  A level offset that is ~0 at both
contacts and bulges inside the depletion region is a VOLUME (generation/transport)
effect; an offset already present at a contact is a BOUNDARY/contact-reconstruction
effect.
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

CM3_TO_M3 = 1.0e6  # cm^-3 s^-1 -> m^-3 s^-1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--junction-x-um", type=float, default=1.140625)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_scalar_csv(path: Path, scale: float, node_count: int) -> list[float]:
    values: dict[int, float] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        value_col = next((n for n in reader.fieldnames or [] if n != "node_id"), None)
        if value_col is None:
            raise RuntimeError(f"{path} has no value column")
        for row in reader:
            values[int(row["node_id"])] = float(row[value_col]) * scale
    return [values.get(i, 0.0) for i in range(node_count)]


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


def signed_ratio(a: float, b: float) -> float | None:
    if b == 0.0 or not math.isfinite(a) or not math.isfinite(b):
        return None
    return a / b


def central_row(nodes: dict[int, dict[str, float]]) -> list[int]:
    by_y: dict[float, list[int]] = defaultdict(list)
    for nid, node in nodes.items():
        by_y[round(node["y_um"], 6)].append(nid)
    max_count = max(len(ids) for ids in by_y.values())
    full = {y: ids for y, ids in by_y.items() if len(ids) >= 0.8 * max_count}
    median_y = statistics.median(sorted(full.keys()))
    best_y = min(full.keys(), key=lambda y: abs(y - median_y))
    return sorted(full[best_y], key=lambda n: nodes[n]["x_um"])


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    nodes, _, contact_by_node = sgdiag.read_mesh(args.mesh)
    node_count = len(nodes)
    vt = sgdiag.K_B_OVER_Q * args.temperature_k

    vela = fluxforms.load_vela_state(
        fluxforms.discover_vela_vtks(args.vela_vtk_root)[fluxforms.bias_key(args.bias)],
        node_count)
    sent = fluxforms.load_sentaurus_state(args.sentaurus_dir, node_count)

    vela_scalars = sgdiag.parse_vtk_scalars(
        fluxforms.discover_vela_vtks(args.vela_vtk_root)[fluxforms.bias_key(args.bias)])
    vela_srh = vela_scalars["SRHRecombination"]
    vela_aval = vela_scalars["AvalancheGeneration"]
    sent_srh = read_scalar_csv(
        args.sentaurus_dir / "fields" / "srhRecombination_region0.csv", CM3_TO_M3, node_count)
    sent_aval = read_scalar_csv(
        args.sentaurus_dir / "fields" / "ImpactIonization_region0.csv", CM3_TO_M3, node_count)

    support = [
        int(row["node_id"]) for row in read_csv_rows(args.support_csv)
        if row.get("support_class") not in (None, "", "inactive")
    ]

    # ---- per-node ratios on the active/depletion support set ----------------
    srh_ratio, aval_ratio = [], []
    n_ratio, p_ratio = [], []
    level_n_delta_mV, level_p_delta_mV = [], []
    srh_v_vals, srh_s_vals = [], []
    for nid in support:
        srh_ratio.append(signed_ratio(vela_srh[nid], sent_srh[nid]))
        aval_ratio.append(signed_ratio(vela_aval[nid], sent_aval[nid]))
        if sent["n"][nid] > 0:
            n_ratio.append(vela["n"][nid] / sent["n"][nid])
        if sent["p"][nid] > 0:
            p_ratio.append(vela["p"][nid] / sent["p"][nid])
        level_n_delta_mV.append(((vela["psi"][nid] - vela["phin"][nid])
                                 - (sent["psi"][nid] - sent["phin"][nid])) * 1.0e3)
        level_p_delta_mV.append(((vela["phip"][nid] - vela["psi"][nid])
                                 - (sent["phip"][nid] - sent["psi"][nid])) * 1.0e3)
        srh_v_vals.append(vela_srh[nid])
        srh_s_vals.append(sent_srh[nid])

    # ---- 1D junction-normal cut from anode to cathode -----------------------
    row = central_row(nodes)
    x_contact = {name: [] for name in {n for s in contact_by_node.values() for n in s}}
    for nid, names in contact_by_node.items():
        for name in names:
            x_contact[name].append(nodes[nid]["x_um"])
    contact_x = {name: (min(xs), max(xs)) for name, xs in x_contact.items()}

    profile_fields = [
        "node_id", "x_um", "dist_to_junction_um", "is_contact",
        "level_psi_minus_phin_delta_mV", "level_phip_minus_psi_delta_mV",
        "vela_srh", "sent_srh", "srh_ratio",
        "vela_aval", "sent_aval", "aval_ratio",
        "n_ratio", "p_ratio",
    ]
    profile_rows: list[dict[str, Any]] = []
    for nid in row:
        names = contact_by_node.get(nid, set())
        profile_rows.append({
            "node_id": nid,
            "x_um": nodes[nid]["x_um"],
            "dist_to_junction_um": nodes[nid]["x_um"] - args.junction_x_um,
            "is_contact": ";".join(sorted(names)) if names else "",
            "level_psi_minus_phin_delta_mV":
                ((vela["psi"][nid] - vela["phin"][nid])
                 - (sent["psi"][nid] - sent["phin"][nid])) * 1.0e3,
            "level_phip_minus_psi_delta_mV":
                ((vela["phip"][nid] - vela["psi"][nid])
                 - (sent["phip"][nid] - sent["psi"][nid])) * 1.0e3,
            "vela_srh": vela_srh[nid],
            "sent_srh": sent_srh[nid],
            "srh_ratio": signed_ratio(vela_srh[nid], sent_srh[nid]),
            "vela_aval": vela_aval[nid],
            "sent_aval": sent_aval[nid],
            "aval_ratio": signed_ratio(vela_aval[nid], sent_aval[nid]),
            "n_ratio": (vela["n"][nid] / sent["n"][nid]) if sent["n"][nid] > 0 else None,
            "p_ratio": (vela["p"][nid] / sent["p"][nid]) if sent["p"][nid] > 0 else None,
        })

    profile_path = args.out_dir / "seed_localization_cut.csv"
    with profile_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=profile_fields)
        writer.writeheader()
        writer.writerows(profile_rows)

    # level offset at the contact-adjacent cut nodes vs the depletion interior
    contact_nodes = [r for r in profile_rows if r["is_contact"]]
    interior = [r for r in profile_rows
                if not r["is_contact"] and abs(r["dist_to_junction_um"]) <= 0.25]

    summary: dict[str, Any] = {
        "bias_V": args.bias,
        "vt_V": vt,
        "support_node_count": len(support),
        "contacts_x_um": contact_x,
        "srh_sign_convention": {
            "vela_srh_median_on_support": statistics.median(srh_v_vals),
            "sentaurus_srh_median_on_support": statistics.median(srh_s_vals),
            "note": "Negative => net generation (np<ni^2) in depletion.",
        },
        "channel_ratios_on_support": {
            "srh_rate_vela_over_sentaurus": stat_block([r for r in srh_ratio if r is not None]),
            "avalanche_gen_vela_over_sentaurus": stat_block([r for r in aval_ratio if r is not None]),
            "electron_density_vela_over_sentaurus": stat_block(n_ratio),
            "hole_density_vela_over_sentaurus": stat_block(p_ratio),
            "level_psi_minus_phin_delta_mV": stat_block(level_n_delta_mV),
            "level_phip_minus_psi_delta_mV": stat_block(level_p_delta_mV),
        },
        "level_offset_origin": {
            "at_contacts_psi_minus_phin_mV":
                stat_block([r["level_psi_minus_phin_delta_mV"] for r in contact_nodes]),
            "at_contacts_phip_minus_psi_mV":
                stat_block([r["level_phip_minus_psi_delta_mV"] for r in contact_nodes]),
            "depletion_interior_psi_minus_phin_mV":
                stat_block([r["level_psi_minus_phin_delta_mV"] for r in interior]),
            "depletion_interior_phip_minus_psi_mV":
                stat_block([r["level_phip_minus_psi_delta_mV"] for r in interior]),
            "note": "If contact offsets ~0 and interior offsets ~ -6/-9 mV the "
                    "seed offset is a VOLUME (generation/transport) effect, not a "
                    "boundary/contact-reconstruction gauge offset.",
        },
    }

    out_json = args.out_dir / "seed_localization_summary.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"\nWrote {out_json}")
    print(f"Wrote {profile_path}")


if __name__ == "__main__":
    main()
