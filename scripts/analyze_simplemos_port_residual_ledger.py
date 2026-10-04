"""Read-only re-analysis of the frozen M46/M54 ledgers (no new solves).

Box-method hypothesis evaluated throughout (Dirichlet contacts, doping-well terminal algorithm):

    I_default(c) - I_direct(c) = sum_{i in W_c, i not a contact node} r_i        (per carrier)
    sum_c I_default,total(c)   = -sum_{i in no contact well} (r_i^e + r_i^h)

where W_c is the connected same-doping-sign node set of contact c and r_i would be the
discrete continuity residual left by the Sentaurus Newton solve.  The frozen ledgers do
not export r_i independently.  They directly provide the substrate ``default - direct``
electron-current observable; interpreting it as a p-well residual sum remains a
falsifiable box-method hypothesis for M60.

Inputs (frozen, committed):
  reference_tcad/simplemos_sentaurus2022/terminal_common_mode_attribution/m54_terminal_component_ledger.csv
  reference_tcad/simplemos_sentaurus2022/terminal_common_mode_attribution/m54_state_decomposition_ledger.csv
  reference_tcad/simplemos_sentaurus2022/full_matrix_requalification/m46_pointwise_delta.csv
"""
from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "reference_tcad" / "simplemos_sentaurus2022"
LEDGER = ROOT / "terminal_common_mode_attribution" / "m54_terminal_component_ledger.csv"
STATE = ROOT / "terminal_common_mode_attribution" / "m54_state_decomposition_ledger.csv"
DELTA = ROOT / "full_matrix_requalification" / "m46_pointwise_delta.csv"

BURST_FRACTION = 0.02   # |default-direct substrate eCurrent| / Id_S burst threshold
CLEAN_FRACTION = 0.005  # neighbours below this are used to estimate the smooth baseline


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return sxy / math.sqrt(sxx * syy)


def lagrange_excluded(xs, ys, k):
    """Value at xs[k] of the cubic through the four nearest neighbours (k excluded), minus ys[k]."""
    idx = [j for j in (k - 2, k - 1, k + 1, k + 2) if 0 <= j < len(xs)]
    if len(idx) < 4:
        return None
    acc = 0.0
    for j in idx:
        w = 1.0
        for m in idx:
            if m != j:
                w *= (xs[k] - xs[m]) / (xs[j] - xs[m])
        acc += w * ys[j]
    return ys[k] - acc


def load():
    comp = defaultdict(dict)
    meta = {}
    with LEDGER.open(newline="") as f:
        for r in csv.DictReader(f):
            comp[r["state"]][(r["contact"], r["component"])] = (
                float(r["default_A_per_um"]), float(r["direct_A_per_um"]))
            meta[r["state"]] = (r["case"], r["device"], float(r["drain_voltage_V"]),
                                round(float(r["gate_voltage_V"]), 3))
    kcl = {}
    with STATE.open(newline="") as f:
        for r in csv.DictReader(f):
            kcl[r["state"]] = (float(r["total_default_four_terminal_sum_A_per_um"]),
                               float(r["total_direct_four_terminal_sum_A_per_um"]))
    cur = {}
    with DELTA.open(newline="") as f:
        for r in csv.DictReader(f):
            cur[(r["case"], round(float(r["gate_voltage_V"]), 3))] = (
                float(r["sentaurus_current_A_per_um"]), float(r["m46_vela_current_A_per_um"]))

    rows = []
    for st, (case, dev, vd, vg) in meta.items():
        sub_e_def, sub_e_dir = comp[st][("substrate", "electron")]
        drn_def, _ = comp[st][("drain", "total")]
        s_id, v_id = cur[(case, vg)]
        if abs(s_id - drn_def) > 1e-12 * abs(s_id):
            raise RuntimeError(f"{st}: M46 Sentaurus drain {s_id} != M54 default drain {drn_def}")
        rows.append(dict(state=st, case=case, dev=dev, vd=vd, vg=vg, s_id=s_id, v_id=v_id,
                         dex=math.log10(abs(v_id) / abs(s_id)), dI=v_id - s_id,
                         sub_e_def=sub_e_def, sub_e_dir=sub_e_dir,
                         observable_gap_e=sub_e_def - sub_e_dir,
                         kcl_def=kcl[st][0], kcl_dir=kcl[st][1]))
    rows.sort(key=lambda r: (r["case"], r["vg"]))
    bycase = defaultdict(list)
    for r in rows:
        bycase[r["case"]].append(r)
    return rows, bycase


def main() -> None:
    rows, bycase = load()

    print("=== A. low-Vg states: substrate default-minus-Direct electron-current observable ===")
    print(f"{'case':<12}{'Vg':>5}{'Id_S':>11}{'Id_V':>11}{'dex':>8}{'subE_def':>11}{'subE_dir':>11}"
          f"{'obs_gap':>11}{'gap/Id_S':>9}")
    for r in rows:
        if r["vg"] <= 0.25:
            print(f"{r['case']:<12}{r['vg']:>5.2f}{r['s_id']:>11.2e}{r['v_id']:>11.2e}{r['dex']:>8.4f}"
                  f"{r['sub_e_def']:>11.2e}{r['sub_e_dir']:>11.2e}{r['observable_gap_e']:>11.2e}"
                  f"{r['observable_gap_e'] / r['s_id']:>9.3f}")

    print(f"\n=== B. flagged points (|observable_gap|/Id_S >= {BURST_FRACTION}): excess drain error association ===")
    print(f"{'case':<12}{'Vg':>5}{'gap/Id':>8}{'dex':>8}{'base':>8}{'excess':>8}{'excess_lin':>11}{'ratio':>7}")
    xs, ys = [], []
    for case, rs in sorted(bycase.items()):
        for r in rs:
            frac = r["observable_gap_e"] / r["s_id"]
            if abs(frac) < BURST_FRACTION:
                continue
            clean = [q for q in rs if q is not r and abs(q["vg"] - r["vg"]) <= 0.1001
                     and abs(q["observable_gap_e"] / q["s_id"]) < CLEAN_FRACTION]
            if not clean:
                continue
            base = sum(q["dex"] for q in clean) / len(clean)
            excess_lin = 10 ** (r["dex"] - base) - 1.0
            xs.append(frac)
            ys.append(excess_lin)
            print(f"{case:<12}{r['vg']:>5.2f}{frac:>8.3f}{r['dex']:>8.4f}{base:>8.4f}{r['dex'] - base:>8.4f}"
                  f"{excess_lin:>11.4f}{excess_lin / frac:>7.3f}")
    slope = sum(a * b for a, b in zip(xs, ys)) / sum(a * a for a in xs)
    print(f"least-squares excess_lin = k * (observable_gap/Id): k = {slope:.3f} "
          f"(n={len(xs)}, Pearson = {pearson(xs, ys):.3f}); report drain scopes separately")

    print("\n=== C. default-method KCL residual (sum of four default totals): floating-well leakage test ===")
    print(f"{'case':<12}{'med|KCL_def|':>13}{'max|KCL_def|':>13}{'med|KCL_dir|':>13}")
    for case, rs in sorted(bycase.items()):
        kd = sorted(abs(r["kcl_def"]) for r in rs)
        kr = sorted(abs(r["kcl_dir"]) for r in rs)
        print(f"{case:<12}{kd[len(kd) // 2]:>13.2e}{kd[-1]:>13.2e}{kr[len(kr) // 2]:>13.2e}")

    print("\n=== D. engine-internal roughness at Vg<=0.30: log10(Id) minus 4-neighbour cubic (dex) ===")
    print(f"{'case':<12}{'Vg':>5}{'rough_S':>9}{'rough_V':>9}{'gap/Id_S':>9}")
    for case, rs in sorted(bycase.items()):
        vgs = [r["vg"] for r in rs]
        ls = [math.log10(abs(r["s_id"])) for r in rs]
        lv = [math.log10(abs(r["v_id"])) for r in rs]
        for k, r in enumerate(rs):
            if r["vg"] > 0.30:
                continue
            ds, dv = lagrange_excluded(vgs, ls, k), lagrange_excluded(vgs, lv, k)
            if ds is None:
                continue
            flag = " <-- burst" if abs(r["observable_gap_e"] / r["s_id"]) >= BURST_FRACTION else ""
            print(f"{case:<12}{r['vg']:>5.2f}{ds:>9.4f}{dv:>9.4f}{r['observable_gap_e'] / r['s_id']:>9.3f}{flag}")

    print("\n=== E. per-curve maximum: all points vs burst-free points ===")
    print(f"{'case':<12}{'max_all':>8}{'Vg':>5}{'burst':>6}{'max_clean':>10}{'Vg':>5}"
          f"{'Id_S@0.05':>11}{'gap@0.05':>11}{'gap/Id':>7}{'Id_S@0':>10}")
    clean_max = {}
    for case, rs in sorted(bycase.items()):
        m = max(rs, key=lambda r: abs(r["dex"]))
        cl = [r for r in rs if abs(r["observable_gap_e"] / r["s_id"]) < BURST_FRACTION]
        mc = max(cl, key=lambda r: abs(r["dex"]))
        clean_max[case] = mc["dex"]
        r05 = next(r for r in rs if abs(r["vg"] - 0.05) < 1e-9)
        r00 = next(r for r in rs if abs(r["vg"]) < 1e-9)
        burst = "yes" if abs(m["observable_gap_e"] / m["s_id"]) >= BURST_FRACTION else "no"
        print(f"{case:<12}{m['dex']:>8.4f}{m['vg']:>5.2f}{burst:>6}{mc['dex']:>10.4f}{mc['vg']:>5.2f}"
              f"{r05['s_id']:>11.2e}{r05['observable_gap_e']:>11.2e}{r05['observable_gap_e'] / r05['s_id']:>7.3f}{r00['s_id']:>10.2e}")

    print("\n=== F. NWell-pair growth of the burst-free maximum ===")
    for lo, hi in (("n17", "n21"), ("n18", "n22"), ("n19", "n23"), ("n20", "n24")):
        for vd in ("vd_0p05", "vd_1"):
            a, b = clean_max[f"{lo}_{vd}"], clean_max[f"{hi}_{vd}"]
            print(f"{lo}/{hi} {vd:<8} {a:.4f} -> {b:.4f}   growth {b - a:+.4f}")

    print("\n=== G. smooth error levels per curve (median dex in Vg windows) ===")
    print(f"{'case':<12}{'0.25-0.50':>10}{'0.75-0.85':>10}{'2.50':>7}")
    for case, rs in sorted(bycase.items()):
        w1 = sorted(r["dex"] for r in rs if 0.25 <= r["vg"] <= 0.5)
        w2 = sorted(r["dex"] for r in rs if 0.75 <= r["vg"] <= 0.85)
        e = next(r["dex"] for r in rs if abs(r["vg"] - 2.5) < 1e-9)
        print(f"{case:<12}{w1[len(w1) // 2]:>10.4f}{w2[len(w2) // 2]:>10.4f}{e:>7.4f}")


if __name__ == "__main__":
    main()
