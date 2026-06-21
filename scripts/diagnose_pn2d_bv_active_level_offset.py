#!/usr/bin/env python3
"""Step 5: localize the active-region carrier-density deficit as an absolute
quasi-Fermi-to-potential level offset, and split it into the electrostatic part
and the ni/BGN part.

Step 4 showed the active high-field avalanche flux deficit (~0.73x) is driven by
carrier density (electron ~0.80x, hole ~0.71x), while every gradient term
(dpsi/h, dphin/h, dphip/h) and mobility match Sentaurus to <2%. Because
``n = ni * exp((psi - phin) / Vt)`` and ``p = ni * exp((phip - psi) / Vt)``, a
multiplicative density deficit at matched gradients must come from an absolute
*level* offset in ``(psi - phin)`` / ``(phip - psi)`` and/or an ``ni`` (BGN)
difference. This script measures, per active node:

* the electrostatic level delta ``d(psi - phin)`` and ``d(phip - psi)`` in mV,
* the density-implied level delta ``Vt * ln(n_V / n_S)`` and ``Vt * ln(p_V / p_S)``,
* their difference, which is exactly ``Vt * ln(ni_V / ni_S)`` (the BGN/ni part),
* the raw ``dpsi``, ``dphin``, ``dphip`` to expose any common-mode gauge shift
  (which cancels in densities) versus a physical split.

It then reports the distribution versus distance to the junction plane to decide
whether the offset is a uniform constant (global reference/gauge or ni constant)
or a junction-localized physical offset.
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--junction-x-um", type=float, default=1.140625)
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


def vt_ln_ratio(a: float, b: float, vt: float) -> float | None:
    if a <= 0.0 or b <= 0.0:
        return None
    return vt * math.log(a / b)


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
    ni_model = fluxforms.effective_ni_from_doping(
        args.doping_csv, node_count, args.material_ni_m3, vt, args.bandgap_narrowing)

    support = [
        int(row["node_id"]) for row in read_csv_rows(args.support_csv)
        if row.get("support_class") not in (None, "", "inactive")
    ]

    fields = [
        "node_id", "x_um", "y_um", "dist_to_junction_um",
        "dpsi_mV", "dphin_mV", "dphip_mV",
        "level_psi_minus_phin_delta_mV", "level_phip_minus_psi_delta_mV",
        "density_level_electron_mV", "density_level_hole_mV",
        "ni_term_electron_mV", "ni_term_hole_mV",
        "electron_density_ratio_V_over_S", "hole_density_ratio_V_over_S",
    ]

    rows_out: list[dict[str, Any]] = []
    for nid in support:
        node = nodes[nid]
        psi_v, phin_v, phip_v = vela["psi"][nid], vela["phin"][nid], vela["phip"][nid]
        psi_s, phin_s, phip_s = sent["psi"][nid], sent["phin"][nid], sent["phip"][nid]
        n_v, p_v = vela["n"][nid], vela["p"][nid]
        n_s, p_s = sent["n"][nid], sent["p"][nid]

        level_n_v = psi_v - phin_v
        level_n_s = psi_s - phin_s
        level_p_v = phip_v - psi_v
        level_p_s = phip_s - psi_s

        dens_level_n = vt_ln_ratio(n_v, n_s, vt)
        dens_level_p = vt_ln_ratio(p_v, p_s, vt)
        level_n_delta = level_n_v - level_n_s
        level_p_delta = level_p_v - level_p_s

        rows_out.append({
            "node_id": nid,
            "x_um": node["x_um"],
            "y_um": node["y_um"],
            "dist_to_junction_um": node["x_um"] - args.junction_x_um,
            "dpsi_mV": (psi_v - psi_s) * 1.0e3,
            "dphin_mV": (phin_v - phin_s) * 1.0e3,
            "dphip_mV": (phip_v - phip_s) * 1.0e3,
            "level_psi_minus_phin_delta_mV": level_n_delta * 1.0e3,
            "level_phip_minus_psi_delta_mV": level_p_delta * 1.0e3,
            "density_level_electron_mV": None if dens_level_n is None else dens_level_n * 1.0e3,
            "density_level_hole_mV": None if dens_level_p is None else dens_level_p * 1.0e3,
            "ni_term_electron_mV": None if dens_level_n is None else (dens_level_n - level_n_delta) * 1.0e3,
            "ni_term_hole_mV": None if dens_level_p is None else (dens_level_p - level_p_delta) * 1.0e3,
            "electron_density_ratio_V_over_S": None if n_s <= 0 else n_v / n_s,
            "hole_density_ratio_V_over_S": None if p_s <= 0 else p_v / p_s,
        })

    detail_path = args.out_dir / "level_offset_nodes.csv"
    with detail_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows_out)

    def col(name: str) -> list[float]:
        return [r[name] for r in rows_out if r[name] is not None]

    summary: dict[str, Any] = {
        "bias_V": args.bias,
        "junction_x_um": args.junction_x_um,
        "active_node_count": len(rows_out),
        "vt_V": vt,
        "gauge_check_common_mode": {
            "dpsi_mV": stat_block(col("dpsi_mV")),
            "dphin_mV": stat_block(col("dphin_mV")),
            "dphip_mV": stat_block(col("dphip_mV")),
            "note": "If dpsi==dphin==dphip the offset is a common-mode gauge shift "
                    "that cancels in densities; divergence is physical.",
        },
        "electrostatic_level_offset": {
            "level_psi_minus_phin_delta_mV": stat_block(col("level_psi_minus_phin_delta_mV")),
            "level_phip_minus_psi_delta_mV": stat_block(col("level_phip_minus_psi_delta_mV")),
        },
        "density_implied_level_offset": {
            "density_level_electron_mV": stat_block(col("density_level_electron_mV")),
            "density_level_hole_mV": stat_block(col("density_level_hole_mV")),
        },
        "ni_bgn_term": {
            "ni_term_electron_mV": stat_block(col("ni_term_electron_mV")),
            "ni_term_hole_mV": stat_block(col("ni_term_hole_mV")),
        },
        "density_ratio_V_over_S": {
            "electron": stat_block(col("electron_density_ratio_V_over_S")),
            "hole": stat_block(col("hole_density_ratio_V_over_S")),
        },
    }

    # Junction-distance bins to test uniform-vs-localized.
    bins = [(-1.0, -0.25), (-0.25, -0.1), (-0.1, -0.03), (-0.03, 0.03),
            (0.03, 0.1), (0.1, 0.25), (0.25, 1.0)]
    bin_rows: list[dict[str, Any]] = []
    for lo, hi in bins:
        subset = [r for r in rows_out if lo <= r["dist_to_junction_um"] < hi]
        if not subset:
            continue
        def med(name: str) -> float | None:
            vals = [r[name] for r in subset if r[name] is not None]
            return statistics.median(vals) if vals else None
        bin_rows.append({
            "dist_lo_um": lo,
            "dist_hi_um": hi,
            "node_count": len(subset),
            "median_level_psi_minus_phin_delta_mV": med("level_psi_minus_phin_delta_mV"),
            "median_level_phip_minus_psi_delta_mV": med("level_phip_minus_psi_delta_mV"),
            "median_ni_term_electron_mV": med("ni_term_electron_mV"),
            "median_electron_density_ratio": med("electron_density_ratio_V_over_S"),
            "median_hole_density_ratio": med("hole_density_ratio_V_over_S"),
        })
    summary["junction_distance_bins"] = bin_rows

    summary_path = args.out_dir / "level_offset_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))

    print(json.dumps(summary, indent=2))
    print(f"\nWrote {summary_path}")
    print(f"Wrote {detail_path}")


if __name__ == "__main__":
    main()
