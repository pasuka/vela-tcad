#!/usr/bin/env python3
"""Quantitative residual table linking pn2d density mismatch to quasi-Fermi drift.

For a single reverse bias this script proves, node by node, that the Vela vs.
Sentaurus carrier-density discrepancy is a pure exponential image of the
minority quasi-Fermi (QF) potential difference, not an electrostatic-potential
or doping error.

Physics used (Boltzmann statistics, identical n_i assumed for both codes):

    n = n_i * exp( (psi - phi_n) / Vt )
    p = n_i * exp( (phi_p - psi) / Vt )

Hence the density ratio (Vela / Sentaurus) at a node is predicted purely from
the state variables:

    ratio_e_pred = exp( ( (psi_v - phin_v) - (psi_s - phin_s) ) / Vt )
    ratio_p_pred = exp( ( (phip_v - psi_v) - (phip_s - psi_s) ) / Vt )

The residual is the measured density ratio divided by this prediction; a value
near 1.0 means the density difference is *fully* explained by the QF/psi state
difference (i.e. the discretization sets a slightly different minority QF, and
the density just follows exponentially).

Input: the ``coarse_node_field_compare_aligned.csv`` produced by the coarse7x3
comparison report. Pure standard library; no numpy required.

Units: potential / QF [V], density [cm^-3], Vt [V].
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import defaultdict

# Boltzmann constant over elementary charge [V/K].
K_B_OVER_Q = 8.617333262145e-5

# Quantity labels as they appear in the aligned CSV.
Q_POTENTIAL = "potential"
Q_EDENS = "electron_density"
Q_HDENS = "hole_density"
Q_EQF = "electron_qf"
Q_HQF = "hole_qf"
WANTED = {Q_POTENTIAL, Q_EDENS, Q_HDENS, Q_EQF, Q_HQF}

DEFAULT_CSV = os.path.join(
    "build-release",
    "reference_tcad",
    "pn2d_sentaurus2018_coarse7x3",
    "reports",
    "coarse_previous_full20_vector_current_20260630",
    "coarse_node_field_compare_aligned.csv",
)


def _fnum(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return float("nan")


def load_bias_nodes(path, bias):
    """Return {node_id: {"x","y", quantity: {"s": sen, "v": vela}}} for one bias."""
    nodes = defaultdict(dict)
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = [c.strip() for c in next(reader)]
        idx = {name: i for i, name in enumerate(header)}
        for row in reader:
            if len(row) < len(header):
                continue
            if abs(_fnum(row[idx["bias_V"]]) - bias) > 1e-9:
                continue
            quantity = row[idx["quantity"]].strip()
            if quantity not in WANTED:
                continue
            nid = int(_fnum(row[idx["node_id"]]))
            entry = nodes[nid]
            entry.setdefault("x", _fnum(row[idx["x_um"]]))
            entry.setdefault("y", _fnum(row[idx["y_um"]]))
            entry[quantity] = {
                "s": _fnum(row[idx["sentaurus_value"]]),
                "v": _fnum(row[idx["vela_value_scaled_to_sentaurus_units"]]),
            }
    return nodes


def safe_ratio(vela, sen):
    if sen == 0 or math.isnan(sen) or math.isnan(vela):
        return float("nan")
    return vela / sen


def build_rows(nodes, vt):
    rows = []
    for nid in sorted(nodes):
        n = nodes[nid]
        psi = n.get(Q_POTENTIAL, {})
        edens = n.get(Q_EDENS, {})
        hdens = n.get(Q_HDENS, {})
        eqf = n.get(Q_EQF, {})
        hqf = n.get(Q_HQF, {})

        d_psi = psi.get("v", float("nan")) - psi.get("s", float("nan"))
        d_phin = eqf.get("v", float("nan")) - eqf.get("s", float("nan"))
        d_phip = hqf.get("v", float("nan")) - hqf.get("s", float("nan"))

        # exponent arguments for the Boltzmann prediction
        arg_e = (d_psi - d_phin) / vt
        arg_p = (d_phip - d_psi) / vt
        ratio_e_pred = math.exp(arg_e) if abs(arg_e) < 700 else float("inf")
        ratio_p_pred = math.exp(arg_p) if abs(arg_p) < 700 else float("inf")

        ratio_e_meas = safe_ratio(edens.get("v", float("nan")), edens.get("s", float("nan")))
        ratio_p_meas = safe_ratio(hdens.get("v", float("nan")), hdens.get("s", float("nan")))

        resid_e = safe_ratio(ratio_e_meas, ratio_e_pred)
        resid_p = safe_ratio(ratio_p_meas, ratio_p_pred)

        rows.append(
            {
                "node": nid,
                "x": n.get("x", float("nan")),
                "y": n.get("y", float("nan")),
                "d_psi_mV": d_psi * 1e3,
                "d_phin_mV": d_phin * 1e3,
                "d_phip_mV": d_phip * 1e3,
                "ratio_e_meas": ratio_e_meas,
                "ratio_e_pred": ratio_e_pred,
                "resid_e": resid_e,
                "ratio_p_meas": ratio_p_meas,
                "ratio_p_pred": ratio_p_pred,
                "resid_p": resid_p,
            }
        )
    return rows


def fmt(value, spec="{: .4f}"):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "   nan  "
    return spec.format(value)


def print_table(rows, bias, vt):
    print(f"# pn2d density <- quasi-Fermi residual table  (bias = {bias:g} V, Vt = {vt*1e3:.3f} mV)")
    print("# ratio = Vela/Sentaurus ; pred = exp(d(psi-phi)/Vt) ; resid = meas/pred (1.0 = fully explained)")
    print()
    header = (
        f"{'node':>4} {'x':>5} {'y':>5} "
        f"{'dpsi':>7} {'dphin':>7} {'dphip':>7} | "
        f"{'e_meas':>8} {'e_pred':>8} {'e_res':>8} | "
        f"{'p_meas':>8} {'p_pred':>8} {'p_res':>8}"
    )
    print(header)
    print("-" * len(header))
    for r in rows:
        print(
            f"{r['node']:>4} {r['x']:>5.2f} {r['y']:>5.2f} "
            f"{fmt(r['d_psi_mV'], '{: .3f}'):>7} {fmt(r['d_phin_mV'], '{: .3f}'):>7} {fmt(r['d_phip_mV'], '{: .3f}'):>7} | "
            f"{fmt(r['ratio_e_meas']):>8} {fmt(r['ratio_e_pred']):>8} {fmt(r['resid_e']):>8} | "
            f"{fmt(r['ratio_p_meas']):>8} {fmt(r['ratio_p_pred']):>8} {fmt(r['resid_p']):>8}"
        )
    print()
    _print_summary("electron", [r["resid_e"] for r in rows])
    _print_summary("hole", [r["resid_p"] for r in rows])


def _print_summary(label, resids):
    finite = [abs(v - 1.0) for v in resids if not math.isnan(v)]
    if not finite:
        print(f"# {label:>8}: no finite residuals")
        return
    print(
        f"# {label:>8} residual |meas/pred - 1|:"
        f"  max={max(finite):.3e}  mean={sum(finite)/len(finite):.3e}  n={len(finite)}"
    )


def write_csv(rows, path):
    fields = [
        "node", "x", "y", "d_psi_mV", "d_phin_mV", "d_phip_mV",
        "ratio_e_meas", "ratio_e_pred", "resid_e",
        "ratio_p_meas", "ratio_p_pred", "resid_p",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", default=DEFAULT_CSV, help="Aligned comparison CSV.")
    parser.add_argument("--bias", type=float, default=-5.0, help="Bias [V] to tabulate.")
    parser.add_argument("--temperature-k", type=float, default=300.0, help="Temperature [K] for Vt.")
    parser.add_argument("--out-csv", default=None, help="Optional path to write the table as CSV.")
    args = parser.parse_args()

    vt = K_B_OVER_Q * args.temperature_k
    nodes = load_bias_nodes(args.csv, args.bias)
    if not nodes:
        raise SystemExit(f"No rows found for bias {args.bias} V in {args.csv}")

    rows = build_rows(nodes, vt)
    print_table(rows, args.bias, vt)

    if args.out_csv:
        write_csv(rows, args.out_csv)
        print(f"# wrote {args.out_csv}")


if __name__ == "__main__":
    main()
