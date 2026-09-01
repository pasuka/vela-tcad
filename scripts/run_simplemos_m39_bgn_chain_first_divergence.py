#!/usr/bin/env python3
"""SimpleMOS M39: BGN-chain first-divergence audit (read-only).

Motivation (M27 + M29 factorial):
- M27: the HFS on/off *response* agrees between Sentaurus and Vela to 3.93%,
  and the absolute deep-off gap survives HFS-off nearly unchanged (0.483%).
- M29 originally labelled two cells as ``BGN off``.  M39 must first verify the
  model selected by sdevice instead of trusting that label.

The audit is fail-closed and starts one stage before the numeric chain:

    C-1 model selection -> C0 N(doping) -> C1 deltaEg(N)
                        -> C2/C3 nie(deltaEg) -> C4/C5 density/band-edge usage

M39 audits each link on Sentaurus' OWN frozen n23 exports (no cross-mesh
interpolation, no absolute-phin global comparisons; per-node identities only):

  C-1 parse the sdevice log and deck to identify the actual BGN model
  C0  doping.csv vs Donor/AcceptorConcentration fields (input identity);
      missing observables are NOT treated as zero
  C1  Vela's SlotboomBandgapNarrowing::deltaEg applied to the exported doping
      vs the exported BandgapNarrowing field, for candidate N conventions
      {donors+acceptors, |donors-acceptors|, max(donors,acceptors)}
  C2  exported BandGap field vs Eg0 and Eg0-deltaEg hypotheses
  C3  Boltzmann density-reconstruction residuals on the Sentaurus state:
        r_e = ln n - ln nie_vela - (psi-phin)/Vt
        r_h = ln p - ln nie_vela - (phip-psi)/Vt
      with nie_vela = ni_vela*exp(deltaEg_sent/(2 Vt)); joint regression on
      u = (psi-phin)/Vt (thermal-voltage convention) and g = deltaEg/(2 Vt)
      (band-split / BGN-usage in the density chain)
  C4  values at the frozen M26/M32 target nodes and edges

States: M27 full, M27 no_hfs, M29 bgn_on_srh_off, and the M29 cell that
was labelled bgn_off_srh_on.  M29 exports omit donor/acceptor fields, so an
explicit M27 same-mesh input-doping oracle is used only for downstream
diagnostics; C0 remains unscored for those states.

No Sentaurus rerun, no Vela solve, no default-model change.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BASE = REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
OUT_DEFAULT = BASE / "m39_bgn_chain_first_divergence"
PORTABLE_DEFAULT = (REPO / "reference_tcad/simplemos_sentaurus2022"
                    / "bgn_chain_first_divergence")
MESH_JSON = BASE / "m8a_model_ablation/vela/mesh.json"
M29_BUNDLE = BASE / "m29_bgn_srh_factorial/sentaurus_raw/sentaurus_bundle"
M27_BUNDLE = BASE / "m27_cross_solver_response_audit/sentaurus_raw/sentaurus_bundle"

KB = 1.380649e-23
Q = 1.602176634e-19

# Frozen Vela Si material contract (materials_sentaurus2022.json).
NI_VELA_CM3 = 14638914958.767616
EG0_EV = 1.12

# Frozen Vela old_slotboom config (bandgapNarrowingConfig in
# src/physics/BandgapNarrowing.cpp).
BGN_COEF_EV = 9.0e-3
BGN_NREF_CM3 = 1.0e17
BGN_SMOOTHING = 0.5
BGN_OFFSET_EV = 0.0

STATES = {
    "m27_full": {
        "export": BASE / "m27_cross_solver_response_audit/sentaurus_exports/full",
        "deck": M27_BUNDLE / "n23_vd_0p05_vg_0p05_full/n23_vd_0p05_vg_0p05_full_des.cmd",
        "log": M27_BUNDLE / "n23_vd_0p05_vg_0p05_full/n23_vd_0p05_vg_0p05_full.log_des.log",
        "declared_model": "OldSlotboom",
        "expected_model": "OldSlotboom",
        "doping_policy": "native",
    },
    "m27_no_hfs": {
        "export": BASE / "m27_cross_solver_response_audit/sentaurus_exports/no_hfs",
        "deck": M27_BUNDLE / "n23_vd_0p05_vg_0p05_no_hfs/n23_vd_0p05_vg_0p05_no_hfs_des.cmd",
        "log": M27_BUNDLE / "n23_vd_0p05_vg_0p05_no_hfs/n23_vd_0p05_vg_0p05_no_hfs.log_des.log",
        "declared_model": "OldSlotboom",
        "expected_model": "OldSlotboom",
        "doping_policy": "native",
    },
    "m29_bgn_on_srh_off": {
        "export": BASE / "m29_bgn_srh_factorial/sentaurus_exports/bgn_on_srh_off",
        "deck": M29_BUNDLE / "n23_m29_bgn_on_srh_off/n23_m29_bgn_on_srh_off_des.cmd",
        "log": M29_BUNDLE / "n23_m29_bgn_on_srh_off/n23_m29_bgn_on_srh_off.log_des.log",
        "declared_model": "OldSlotboom",
        "expected_model": "OldSlotboom",
        "doping_policy": "qualified_m27_oracle",
    },
    "m29_bgn_off_srh_on": {
        "export": BASE / "m29_bgn_srh_factorial/sentaurus_exports/bgn_off_srh_on",
        "deck": M29_BUNDLE / "n23_m29_bgn_off_srh_on/n23_m29_bgn_off_srh_on_des.cmd",
        "log": M29_BUNDLE / "n23_m29_bgn_off_srh_on/n23_m29_bgn_off_srh_on.log_des.log",
        "declared_model": "none",
        "expected_model": "none",
        "doping_policy": "qualified_m27_oracle",
    },
}
DOPING_ORACLE_EXPORT = STATES["m27_full"]["export"]

# Frozen diagnostic targets (Vela mesh node ids, current n23 topology).
TARGET_NODES = {
    "drain_first_layer": [1089, 1090, 1092, 1099, 1100],
    "m26_upstream_source": [1051, 1061, 1048, 1064, 1070],
    "m32_interior_feedback": [837, 846, 836, 91, 831],
}
TARGET_EDGES = [(2416, 1100, 1101), (2415, 1099, 1102), (2392, 1091, 1092)]


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def field_path(export: Path, name: str, region: int) -> Path:
    return export / "fields" / f"{name}_region{region}.csv"


def load_scalar_field(export: Path, name: str, region: int) -> dict[int, float]:
    path = field_path(export, name, region)
    if not path.exists():
        return {}
    out: dict[int, float] = {}
    for row in read_csv_rows(path):
        out[int(row["node_id"])] = float(row["component0"])
    return out


def sha256_file(path: Path) -> str | None:
    if not path.exists():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_bgn_model(log_path: Path | None) -> dict[str, str | bool | None]:
    """Return the effective sdevice BGN model from the solver log.

    Omitting EffectiveIntrinsicDensity from a deck is not a valid indication
    that BGN is disabled: sdevice 2022.03-SP2 selects Bennett/Wilson by
    default in the M29 runs.  The solver log is therefore authoritative.
    """
    if log_path is None or not log_path.exists():
        return {"status": "unobservable", "model": None, "log_line": None}
    text = log_path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"Bandgap narrowing model:\s*([^\r\n]+)", text)
    if not match:
        return {"status": "unobservable", "model": None, "log_line": None}
    line = match.group(1).strip()
    if "OldSlotboom" in line:
        model = "OldSlotboom"
    elif "Bennett/Wilson" in line:
        model = "BennettWilson"
    elif "without bandgap narrowing" in line.lower():
        model = "none"
    else:
        model = "unknown"
    return {"status": "observed", "model": model, "log_line": line}


def deck_bgn_declaration(deck_path: Path | None) -> dict[str, str | bool | None]:
    if deck_path is None or not deck_path.exists():
        return {"status": "unobservable", "explicit_model": None,
                "has_effective_intrinsic_density": None}
    text = deck_path.read_text(encoding="utf-8", errors="replace")
    has_eid = "EffectiveIntrinsicDensity" in text
    model = "OldSlotboom" if re.search(
        r"EffectiveIntrinsicDensity\s*\(\s*OldSlotboom\s*\)", text) else None
    return {"status": "observed", "explicit_model": model,
            "has_effective_intrinsic_density": has_eid}


def load_doping(export: Path) -> dict[int, tuple[float, float]]:
    return {
        int(row["node_id"]): (float(row["donors_cm3"]),
                              float(row["acceptors_cm3"]))
        for row in read_csv_rows(export / "doping.csv")
    }


def coordinates(export: Path) -> dict[int, tuple[float, float]]:
    return {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in read_csv_rows(export / "nodes.csv")
    }


def coordinate_identity(lhs: dict[int, tuple[float, float]],
                        rhs: dict[int, tuple[float, float]],
                        tol_um: float = 1e-12) -> dict[str, float | int | bool]:
    shared = sorted(set(lhs) & set(rhs))
    max_diff = 0.0
    for node_id in shared:
        max_diff = max(max_diff,
                       abs(lhs[node_id][0] - rhs[node_id][0]),
                       abs(lhs[node_id][1] - rhs[node_id][1]))
    return {
        "same_node_ids": set(lhs) == set(rhs),
        "shared_node_count": len(shared),
        "max_coord_diff_um": max_diff,
        "qualified": set(lhs) == set(rhs) and max_diff <= tol_um,
    }


def vela_delta_eg(n_cm3: float) -> float:
    if n_cm3 <= 0.0 or BGN_COEF_EV == 0.0:
        return 0.0
    x = math.log(n_cm3 / BGN_NREF_CM3)
    delta = BGN_OFFSET_EV + BGN_COEF_EV * (x + math.sqrt(x * x + BGN_SMOOTHING))
    return max(delta, 0.0)


def percentile(sorted_values: list[float], q: float) -> float:
    if not sorted_values:
        return math.nan
    pos = (len(sorted_values) - 1) * q
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return sorted_values[lo]
    return sorted_values[lo] + (sorted_values[hi] - sorted_values[lo]) * (pos - lo)


def stats(values: list[float]) -> dict[str, float]:
    if not values:
        return {"count": 0}
    ordered = sorted(values)
    mean = sum(values) / len(values)
    var = sum((v - mean) ** 2 for v in values) / max(len(values) - 1, 1)
    return {
        "count": len(values),
        "mean": mean,
        "std": math.sqrt(var),
        "min": ordered[0],
        "p05": percentile(ordered, 0.05),
        "median": percentile(ordered, 0.50),
        "p95": percentile(ordered, 0.95),
        "max": ordered[-1],
        "max_abs": max(abs(ordered[0]), abs(ordered[-1])),
    }


def solve_normal(ata: list[list[float]], atb: list[float]) -> list[float]:
    n = len(atb)
    a = [row[:] + [atb[i]] for i, row in enumerate(ata)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) < 1e-300:
            raise ValueError("singular normal system")
        a[col], a[pivot] = a[pivot], a[col]
        inv = 1.0 / a[col][col]
        for j in range(col, n + 1):
            a[col][j] *= inv
        for row in range(n):
            if row != col and a[row][col] != 0.0:
                factor = a[row][col]
                for j in range(col, n + 1):
                    a[row][j] -= factor * a[col][j]
    return [a[i][n] for i in range(n)]


def linear_fit(y: list[float], x: list[float]) -> dict[str, float]:
    """Least-squares y ~= a + b*x."""
    n = len(y)
    if n < 4:
        return {"count": n}
    mx = sum(x) / n
    my = sum(y) / n
    sxx = sum((v - mx) ** 2 for v in x)
    if sxx <= 0.0:
        return {"count": n, "offset_a": my, "slope_b": 0.0}
    b = sum((x[k] - mx) * (y[k] - my) for k in range(n)) / sxx
    a = my - b * mx
    resid = [y[k] - a - b * x[k] for k in range(n)]
    return {
        "count": n,
        "offset_a": a,
        "slope_b": b,
        "rms_after_fit": math.sqrt(sum(v * v for v in resid) / n),
        "max_abs_resid": max(abs(v) for v in resid),
    }


def joint_fit(residual: list[float], u: list[float], g: list[float]) -> dict[str, float]:
    """Least-squares residual ~= a + b*u + c*g; returns coefficients and fit rms."""
    n = len(residual)
    if n < 8:
        return {"count": n}
    cols = [[1.0] * n, u, g]
    ata = [[sum(cols[i][k] * cols[j][k] for k in range(n)) for j in range(3)] for i in range(3)]
    atb = [sum(cols[i][k] * residual[k] for k in range(n)) for i in range(3)]
    try:
        coef = solve_normal(ata, atb)
    except ValueError:
        return {"count": n, "error": float("nan")}
    fitted = [coef[0] + coef[1] * u[k] + coef[2] * g[k] for k in range(n)]
    resid = [residual[k] - fitted[k] for k in range(n)]
    rms_before = math.sqrt(sum(v * v for v in residual) / n)
    rms_after = math.sqrt(sum(v * v for v in resid) / n)
    return {
        "count": n,
        "offset_a": coef[0],
        "slope_u": coef[1],
        "slope_g": coef[2],
        "rms_before": rms_before,
        "rms_after_fit": rms_after,
    }


def analyze_state(name: str, spec: dict, mesh_nodes: list[dict], out_dir: Path,
                  oracle_nodes: dict[int, tuple[float, float]],
                  oracle_doping: dict[int, tuple[float, float]]) -> dict:
    export = spec["export"]
    meta = json.loads((export / "metadata.json").read_text(encoding="utf-8-sig"))
    silicon_regions = [r["index"] for r in meta["regions"]
                       if r.get("material") == "Silicon"]
    if not silicon_regions:
        raise RuntimeError(f"{name}: no Silicon region in metadata")
    region = silicon_regions[0]

    nodes = coordinates(export)
    native_doping = load_doping(export)

    # Mesh-id identity check: export node coordinates vs Vela mesh node ids.
    max_coord_diff = 0.0
    shared = 0
    for node in mesh_nodes:
        nid = int(node["id"])
        if nid in nodes:
            shared += 1
            dx = abs(node["x"] - nodes[nid][0])
            dy = abs(node["y"] - nodes[nid][1])
            max_coord_diff = max(max_coord_diff, dx, dy)
    id_identity = shared == len(mesh_nodes) and max_coord_diff < 1e-9

    donors_f = load_scalar_field(export, "DonorConcentration", region)
    acceptors_f = load_scalar_field(export, "AcceptorConcentration", region)
    bgn = load_scalar_field(export, "BandgapNarrowing", region)
    bandgap = load_scalar_field(export, "BandGap", region)
    psi = load_scalar_field(export, "ElectrostaticPotential", region)
    phin = load_scalar_field(export, "eQuasiFermiPotential", region)
    phip = load_scalar_field(export, "hQuasiFermiPotential", region)
    ndens = load_scalar_field(export, "eDensity", region)
    pdens = load_scalar_field(export, "hDensity", region)
    temp = load_scalar_field(export, "LatticeTemperature", region)

    si_nodes = sorted(psi.keys())
    temps = [temp[i] for i in si_nodes if i in temp]
    t_lattice = sum(temps) / len(temps) if temps else 300.0
    vt = KB * t_lattice / Q

    # --- C-1: effective model selection ----------------------------------
    log_path = spec.get("log")
    deck_path = spec.get("deck")
    effective_model = parse_bgn_model(log_path)
    deck_model = deck_bgn_declaration(deck_path)
    observed_model = effective_model["model"]
    expected_model = spec["expected_model"]
    model_selection = {
        "declared_model": spec["declared_model"],
        "expected_model": expected_model,
        "effective_model": observed_model,
        "matches_expected": observed_model == expected_model,
        "solver_log": str(log_path) if log_path else None,
        "solver_log_sha256": sha256_file(log_path) if log_path else None,
        "solver_log_line": effective_model["log_line"],
        "deck": str(deck_path) if deck_path else None,
        "deck_sha256": sha256_file(deck_path) if deck_path else None,
        "deck_declaration": deck_model,
    }

    # --- C0: doping input identity and observability ---------------------
    fields_complete = (bool(si_nodes) and all(i in donors_f and i in acceptors_f
                                              for i in si_nodes))
    csv_complete = bool(si_nodes) and all(i in native_doping for i in si_nodes)
    csv_nonzero = any(d != 0.0 or a != 0.0 for d, a in native_doping.values())
    c0_observable = fields_complete and csv_complete
    c0_max_rel = None
    if c0_observable:
        c0_max_rel = 0.0
        for i in si_nodes:
            d_csv, a_csv = native_doping[i]
            for lhs, rhs in ((donors_f[i], d_csv), (acceptors_f[i], a_csv)):
                denom = max(abs(lhs), abs(rhs), 1.0)
                c0_max_rel = max(c0_max_rel, abs(lhs - rhs) / denom)

    oracle_identity = coordinate_identity(nodes, oracle_nodes)
    if spec["doping_policy"] == "native":
        if not c0_observable or not csv_nonzero:
            raise RuntimeError(f"{name}: native doping observables are incomplete")
        doping = native_doping
        doping_source = "native_export"
        downstream_doping_qualified = True
    elif spec["doping_policy"] == "qualified_m27_oracle":
        if not oracle_identity["qualified"]:
            raise RuntimeError(f"{name}: M27 doping oracle is not mesh-identical")
        doping = oracle_doping
        doping_source = "m27_full_same_mesh_input_oracle"
        downstream_doping_qualified = True
    else:
        raise RuntimeError(f"{name}: unknown doping policy {spec['doping_policy']}")

    c0 = {
        "status": "scored" if c0_observable else "not_scored_missing_fields",
        "native_doping_csv_complete": csv_complete,
        "native_doping_csv_has_nonzero_value": csv_nonzero,
        "donor_acceptor_fields_complete": fields_complete,
        "max_relative_difference": c0_max_rel,
        "downstream_doping_source": doping_source,
        "downstream_doping_qualified": downstream_doping_qualified,
        "oracle_mesh_identity": oracle_identity,
    }

    # --- C1: deltaEg operator -------------------------------------------
    conventions = {
        "total": lambda d, a: d + a,
        "abs_net": lambda d, a: abs(d - a),
        "max_species": lambda d, a: max(d, a),
    }
    c1_rows = []
    c1_diffs = {k: [] for k in conventions}
    has_bgn_field = len(bgn) > 0
    for i in si_nodes:
        if i not in doping:
            continue
        d_cm3, a_cm3 = doping[i]
        sent = bgn.get(i, 0.0)
        row = {
            "node_id": i,
            "x_um": nodes[i][0],
            "y_um": nodes[i][1],
            "donors_cm3": d_cm3,
            "acceptors_cm3": a_cm3,
            "deltaEg_sentaurus_eV": sent,
        }
        for key, fn in conventions.items():
            model = vela_delta_eg(fn(d_cm3, a_cm3))
            row[f"deltaEg_vela_{key}_eV"] = model
            row[f"diff_{key}_eV"] = model - sent
            c1_diffs[key].append(model - sent)
        c1_rows.append(row)

    c1_stats = {k: stats(v) for k, v in c1_diffs.items()}
    c1_best = min(c1_stats, key=lambda k: c1_stats[k].get("max_abs", math.inf)) \
        if has_bgn_field else None
    if not has_bgn_field:
        c1_status = "not_scored_missing_bandgap_narrowing_field"
    elif observed_model == "OldSlotboom":
        c1_status = "scored_same_model"
    else:
        c1_status = "diagnostic_cross_model_not_parity"

    # --- C2: BandGap convention ------------------------------------------
    c2 = {}
    if bandgap:
        raw = [bandgap[i] for i in si_nodes if i in bandgap]
        plus_bgn = [bandgap[i] + bgn.get(i, 0.0) for i in si_nodes if i in bandgap]
        c2 = {
            "bandgap_raw": stats(raw),
            "bandgap_plus_deltaEg": stats(plus_bgn),
            "eg0_contract_eV": EG0_EV,
        }

    # --- C3: Boltzmann reconstruction residuals ---------------------------
    ln_ni = math.log(NI_VELA_CM3)
    c3_rows = []
    r_e, r_h, u_e, u_h, g_arr = [], [], [], [], []
    for i in si_nodes:
        if any(i not in f for f in (phin, phip, ndens, pdens)):
            continue
        n_i, p_i = ndens[i], pdens[i]
        if n_i <= 0.0 or p_i <= 0.0:
            continue
        d_eg = bgn.get(i, 0.0)
        gval = d_eg / (2.0 * vt)
        ue = (psi[i] - phin[i]) / vt
        uh = (phip[i] - psi[i]) / vt
        re = math.log(n_i) - (ln_ni + gval) - ue
        rh = math.log(p_i) - (ln_ni + gval) - uh
        r_e.append(re)
        r_h.append(rh)
        u_e.append(ue)
        u_h.append(uh)
        g_arr.append(gval)
        c3_rows.append({
            "node_id": i,
            "x_um": nodes[i][0],
            "y_um": nodes[i][1],
            "deltaEg_eV": d_eg,
            "u_e": ue,
            "u_h": uh,
            "r_e": re,
            "r_h": rh,
            "r_sum_ref_free": re + rh,
            "r_diff_split": re - rh,
        })

    c3 = {
        "thermal_voltage_V": vt,
        "lattice_temperature_K": t_lattice,
        "r_e": stats(r_e),
        "r_h": stats(r_h),
        "r_sum_ref_free": stats([a + b for a, b in zip(r_e, r_h)]),
        "r_diff_split": stats([a - b for a, b in zip(r_e, r_h)]),
        "fit_r_e": joint_fit(r_e, u_e, g_arr),
        "fit_r_h": joint_fit(r_h, u_h, g_arr),
    }

    # --- C5: band-edge split of deltaEg (conduction share) -----------------
    # sdevice exports Ec/Ev per node.  Writing Ec + psi = k_c - a*deltaEg and
    # Ev + psi = k_v + b*deltaEg, the fitted 'a' is the conduction-band share
    # of the doping-dependent narrowing actually used by Sentaurus transport.
    # Vela's VariableNi SG drift implies a symmetric a = b = 0.5.
    ec = load_scalar_field(export, "ConductionBandEnergy", region)
    ev = load_scalar_field(export, "ValenceBandEnergy", region)
    c5 = {}
    if ec and ev:
        ids = [i for i in si_nodes if i in ec and i in ev and i in psi]
        d_eg_arr = [bgn.get(i, 0.0) for i in ids]
        ec_psi = [ec[i] + psi[i] for i in ids]
        ev_psi = [ev[i] + psi[i] for i in ids]
        gap = [ec[i] - ev[i] for i in ids]
        fit_c = linear_fit(ec_psi, d_eg_arr)
        fit_v = linear_fit(ev_psi, d_eg_arr)
        fit_g = linear_fit(gap, d_eg_arr)
        c5 = {
            "conduction_share_a": -fit_c.get("slope_b", math.nan),
            "valence_share_b": fit_v.get("slope_b", math.nan),
            "ec_fit": fit_c,
            "ev_fit": fit_v,
            "gap_slope_vs_deltaEg": fit_g.get("slope_b", math.nan),
            "gap_offset_eV": fit_g.get("offset_a", math.nan),
        }

    # --- C4: frozen targets ------------------------------------------------
    c3_by_node = {row["node_id"]: row for row in c3_rows}
    c1_by_node = {row["node_id"]: row for row in c1_rows}
    c4_nodes = []
    for group, ids in TARGET_NODES.items():
        for nid in ids:
            entry = {"group": group, "node_id": nid}
            if nid in c1_by_node:
                src = c1_by_node[nid]
                entry.update({
                    "x_um": src["x_um"],
                    "y_um": src["y_um"],
                    "deltaEg_sentaurus_eV": src["deltaEg_sentaurus_eV"],
                    "deltaEg_vela_total_eV": src["deltaEg_vela_total_eV"],
                    "diff_total_eV": src["diff_total_eV"],
                })
            if nid in c3_by_node:
                src = c3_by_node[nid]
                entry.update({"r_e": src["r_e"], "r_h": src["r_h"]})
            c4_nodes.append(entry)

    c4_edges = []
    for edge_id, n0, n1 in TARGET_EDGES:
        entry = {"edge_id": edge_id, "n0": n0, "n1": n1}
        if n0 in c1_by_node and n1 in c1_by_node:
            s0 = c1_by_node[n0]["deltaEg_sentaurus_eV"]
            s1 = c1_by_node[n1]["deltaEg_sentaurus_eV"]
            v0 = c1_by_node[n0]["deltaEg_vela_total_eV"]
            v1 = c1_by_node[n1]["deltaEg_vela_total_eV"]
            entry.update({
                "d_deltaEg_sentaurus_eV": s1 - s0,
                "d_deltaEg_vela_total_eV": v1 - v0,
                "d_ln_nie_sentaurus": (s1 - s0) / (2.0 * vt),
                "d_ln_nie_vela": (v1 - v0) / (2.0 * vt),
            })
        c4_edges.append(entry)

    # --- write per-state artifacts ----------------------------------------
    state_dir = out_dir / name
    state_dir.mkdir(parents=True, exist_ok=True)
    if c1_rows:
        with (state_dir / "c1_delta_eg_nodes.csv").open("w", newline="",
                                                        encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(c1_rows[0]),
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(c1_rows)
    if c3_rows:
        with (state_dir / "c3_reconstruction_nodes.csv").open(
                "w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(c3_rows[0]),
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(c3_rows)

    return {
        "export_dir": str(export),
        "silicon_region": region,
        "silicon_node_count": len(si_nodes),
        "mesh_id_identity": id_identity,
        "mesh_id_max_coord_diff_um": max_coord_diff,
        "c_minus_1_model_selection": model_selection,
        "has_bandgap_narrowing_field": has_bgn_field,
        "c0_doping": c0,
        "c1_status": c1_status,
        "c1_delta_eg_diff_stats_eV": c1_stats,
        "c1_best_convention": c1_best,
        "c2_bandgap": c2,
        "c3_reconstruction": c3,
        "c5_band_edge_split": c5,
        "c4_target_nodes": c4_nodes,
        "c4_target_edges": c4_edges,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUT_DEFAULT)
    parser.add_argument("--portable-dir", type=Path, default=PORTABLE_DEFAULT)
    args = parser.parse_args()

    mesh_nodes = json.loads(MESH_JSON.read_text(encoding="utf-8-sig"))["nodes"]
    oracle_nodes = coordinates(DOPING_ORACLE_EXPORT)
    oracle_doping = load_doping(DOPING_ORACLE_EXPORT)

    summary: dict = {
        "schema": "vela.simplemos.m39_bgn_chain_first_divergence.v2",
        "constants": {
            "kb_J_per_K": KB,
            "q_C": Q,
            "ni_vela_cm3": NI_VELA_CM3,
            "bgn_coefficient_eV": BGN_COEF_EV,
            "bgn_reference_cm3": BGN_NREF_CM3,
            "bgn_smoothing": BGN_SMOOTHING,
        },
        "doping_oracle": {
            "export_dir": str(DOPING_ORACLE_EXPORT),
            "nodes_sha256": sha256_file(DOPING_ORACLE_EXPORT / "nodes.csv"),
            "doping_sha256": sha256_file(DOPING_ORACLE_EXPORT / "doping.csv"),
        },
        "states": {},
    }

    for name, spec in STATES.items():
        export = spec["export"]
        if not export.exists():
            summary["states"][name] = {"skipped": f"missing export {export}"}
            print(f"[M39] SKIP {name}: {export} missing")
            continue
        result = analyze_state(name, spec, mesh_nodes, args.output_dir,
                               oracle_nodes, oracle_doping)
        summary["states"][name] = result

        c1 = result["c1_delta_eg_diff_stats_eV"]
        c3 = result["c3_reconstruction"]
        c0 = result["c0_doping"]
        cm1 = result["c_minus_1_model_selection"]
        c0_value = c0["max_relative_difference"]
        c0_text = f"{c0_value:.3e}" if c0_value is not None else c0["status"]
        print(f"\n[M39] state {name} (Si nodes {result['silicon_node_count']}, "
              f"mesh-id identity {result['mesh_id_identity']}, "
              f"C0 {c0_text})")
        print(f"    C-1 model: declared {cm1['declared_model']}, "
              f"effective {cm1['effective_model']}, "
              f"matches expected {cm1['matches_expected']}")
        print(f"    C0 source: {c0['downstream_doping_source']} "
              f"({c0['status']})")
        if result["has_bandgap_narrowing_field"]:
            for key in ("total", "abs_net", "max_species"):
                st = c1[key]
                print(f"    C1 deltaEg[{key:12s}] max|diff| {st['max_abs']:.6e} eV, "
                      f"median {st['median']:+.6e} eV")
            print(f"    C1 best convention: {result['c1_best_convention']} "
                  f"({result['c1_status']})")
        else:
            print("    C1 skipped: no BandgapNarrowing field (BGN off state)")
        fe, fh = c3["fit_r_e"], c3["fit_r_h"]
        print(f"    C3 Vt {c3['thermal_voltage_V']:.9f} V; "
              f"r_e mean {c3['r_e']['mean']:+.6e} std {c3['r_e']['std']:.3e}; "
              f"r_h mean {c3['r_h']['mean']:+.6e} std {c3['r_h']['std']:.3e}")
        if "slope_u" in fe:
            print(f"    C3 fit r_e: a {fe['offset_a']:+.3e}, b(u) {fe['slope_u']:+.3e}, "
                  f"c(g) {fe['slope_g']:+.3e}, rms {fe['rms_before']:.3e}->"
                  f"{fe['rms_after_fit']:.3e}")
            print(f"    C3 fit r_h: a {fh['offset_a']:+.3e}, b(u) {fh['slope_u']:+.3e}, "
                  f"c(g) {fh['slope_g']:+.3e}, rms {fh['rms_before']:.3e}->"
                  f"{fh['rms_after_fit']:.3e}")
        c5 = result.get("c5_band_edge_split") or {}
        if c5:
            print(f"    C5 conduction share a {c5['conduction_share_a']:+.6f}, "
                  f"valence share b {c5['valence_share_b']:+.6f}, "
                  f"gap slope {c5['gap_slope_vs_deltaEg']:+.6f} "
                  f"(ec fit rms {c5['ec_fit'].get('rms_after_fit', float('nan')):.3e})")

    # Acceptance gates are deliberately asymmetric.  M27 is the qualified
    # OldSlotboom numeric chain.  The M29 "off" cell must expose the model
    # selection mismatch rather than being allowed to masquerade as no-BGN.
    m27_names = ("m27_full", "m27_no_hfs")
    m27_results = [summary["states"].get(name, {}) for name in m27_names]
    m29_off = summary["states"].get("m29_bgn_off_srh_on", {})
    gates = {
        "m27_models_are_oldslotboom": all(
            r.get("c_minus_1_model_selection", {}).get("effective_model") == "OldSlotboom"
            for r in m27_results),
        "m27_c0_exact": all(
            r.get("c0_doping", {}).get("status") == "scored" and
            r.get("c0_doping", {}).get("max_relative_difference", math.inf) <= 1e-15
            for r in m27_results),
        "m27_c1_total_exact": all(
            r.get("c1_delta_eg_diff_stats_eV", {}).get("total", {}).get(
                "max_abs", math.inf) <= 1e-12 for r in m27_results),
        "m27_c3_joint_fit_closes": all(
            max(r.get("c3_reconstruction", {}).get("fit_r_e", {}).get(
                    "rms_after_fit", math.inf),
                r.get("c3_reconstruction", {}).get("fit_r_h", {}).get(
                    "rms_after_fit", math.inf)) <= 1e-10
            for r in m27_results),
        "m27_c5_symmetric_split": all(
            abs(r.get("c5_band_edge_split", {}).get(
                    "conduction_share_a", math.inf) - 0.5) <= 1e-10 and
            abs(r.get("c5_band_edge_split", {}).get(
                    "valence_share_b", math.inf) - 0.5) <= 1e-10
            for r in m27_results),
        "m29_off_label_mismatch_detected": (
            m29_off.get("c_minus_1_model_selection", {}).get("declared_model") == "none" and
            m29_off.get("c_minus_1_model_selection", {}).get("effective_model") ==
            "BennettWilson" and
            not m29_off.get("c_minus_1_model_selection", {}).get(
                "matches_expected", True)),
        "m29_missing_doping_not_silently_scored": (
            m29_off.get("c0_doping", {}).get("status") ==
            "not_scored_missing_fields" and
            m29_off.get("c0_doping", {}).get("downstream_doping_source") ==
            "m27_full_same_mesh_input_oracle"),
    }
    summary["acceptance_gates"] = gates
    summary["all_acceptance_gates_pass"] = all(gates.values())
    summary["first_divergence"] = {
        "stage": "C-1_model_selection_contract",
        "state": "m29_bgn_off_srh_on",
        "declared": "none",
        "observed": m29_off.get("c_minus_1_model_selection", {}).get(
            "effective_model"),
        "finding": (
            "Omitting EffectiveIntrinsicDensity(OldSlotboom) selected the "
            "Sentaurus default Bennett/Wilson BGN model; the M29 cell is not "
            "a no-BGN control."
        ),
        "impact": (
            "M29 terminal currents remain measured data, but the prior "
            "BGN-on versus BGN-off causal label and effect estimate are invalid."
        ),
        "m27_oldslotboom_chain": (
            "C0 input doping, C1 total-impurity OldSlotboom deltaEg, C3 frozen "
            "density identity, and C5 symmetric band split all pass."
        ),
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    ledger_path = args.output_dir / "m39_model_selection_ledger.csv"
    with ledger_path.open("w", newline="", encoding="utf-8") as stream:
        fields = ["state", "declared_model", "expected_model", "effective_model",
                  "matches_expected", "deck_explicit_model", "c0_status",
                  "doping_source", "c1_status"]
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for name, result in summary["states"].items():
            if "skipped" in result:
                continue
            cm1 = result["c_minus_1_model_selection"]
            c0 = result["c0_doping"]
            writer.writerow({
                "state": name,
                "declared_model": cm1["declared_model"],
                "expected_model": cm1["expected_model"],
                "effective_model": cm1["effective_model"],
                "matches_expected": cm1["matches_expected"],
                "deck_explicit_model": cm1["deck_declaration"]["explicit_model"],
                "c0_status": c0["status"],
                "doping_source": c0["downstream_doping_source"],
                "c1_status": result["c1_status"],
            })

    out_json = args.output_dir / "m39_summary.json"
    out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8",
                        newline="\n")
    m27_full = summary["states"]["m27_full"]
    m27_no_hfs = summary["states"]["m27_no_hfs"]
    m29_on = summary["states"]["m29_bgn_on_srh_off"]
    m29_off = summary["states"]["m29_bgn_off_srh_on"]
    portable = {
        "schema": "vela.simplemos.sdevice.m39_bgn_chain_first_divergence_report.v1",
        "status": "complete" if summary["all_acceptance_gates_pass"] else "failed",
        "scope": "read_only_frozen_state_audit_no_solver_rerun",
        "first_divergence": summary["first_divergence"],
        "m27_qualified_oldslotboom_chain": {
            "states": ["m27_full", "m27_no_hfs"],
            "silicon_nodes_per_state": m27_full["silicon_node_count"],
            "mesh_id_identity": (m27_full["mesh_id_identity"] and
                                 m27_no_hfs["mesh_id_identity"]),
            "c0_max_relative_difference": max(
                m27_full["c0_doping"]["max_relative_difference"],
                m27_no_hfs["c0_doping"]["max_relative_difference"]),
            "c1_total_impurity_max_abs_difference_eV": max(
                m27_full["c1_delta_eg_diff_stats_eV"]["total"]["max_abs"],
                m27_no_hfs["c1_delta_eg_diff_stats_eV"]["total"]["max_abs"]),
            "c1_abs_net_max_abs_difference_eV": max(
                m27_full["c1_delta_eg_diff_stats_eV"]["abs_net"]["max_abs"],
                m27_no_hfs["c1_delta_eg_diff_stats_eV"]["abs_net"]["max_abs"]),
            "c1_best_convention": m27_full["c1_best_convention"],
            "c3_max_joint_fit_rms": max(
                m27_full["c3_reconstruction"]["fit_r_e"]["rms_after_fit"],
                m27_full["c3_reconstruction"]["fit_r_h"]["rms_after_fit"],
                m27_no_hfs["c3_reconstruction"]["fit_r_e"]["rms_after_fit"],
                m27_no_hfs["c3_reconstruction"]["fit_r_h"]["rms_after_fit"]),
            "c5_conduction_share": m27_full["c5_band_edge_split"][
                "conduction_share_a"],
            "c5_valence_share": m27_full["c5_band_edge_split"][
                "valence_share_b"],
            "c5_gap_slope": m27_full["c5_band_edge_split"][
                "gap_slope_vs_deltaEg"],
        },
        "m29_contract_audit": {
            "bgn_on_effective_model": m29_on[
                "c_minus_1_model_selection"]["effective_model"],
            "bgn_on_solver_log_line": m29_on[
                "c_minus_1_model_selection"]["solver_log_line"],
            "bgn_on_solver_log_sha256": m29_on[
                "c_minus_1_model_selection"]["solver_log_sha256"],
            "bgn_off_label_effective_model": m29_off[
                "c_minus_1_model_selection"]["effective_model"],
            "bgn_off_solver_log_line": m29_off[
                "c_minus_1_model_selection"]["solver_log_line"],
            "bgn_off_solver_log_sha256": m29_off[
                "c_minus_1_model_selection"]["solver_log_sha256"],
            "bgn_off_deck_sha256": m29_off[
                "c_minus_1_model_selection"]["deck_sha256"],
            "bgn_off_label_matches_expected": m29_off[
                "c_minus_1_model_selection"]["matches_expected"],
            "bgn_off_c1_status": m29_off["c1_status"],
            "native_c0_status": m29_off["c0_doping"]["status"],
            "downstream_doping_source": m29_off["c0_doping"][
                "downstream_doping_source"],
            "oracle_mesh_identity": m29_off["c0_doping"][
                "oracle_mesh_identity"],
            "bennett_wilson_vs_vela_oldslotboom_total_max_abs_difference_eV":
                m29_off["c1_delta_eg_diff_stats_eV"]["total"]["max_abs"],
        },
        "input_provenance": summary["doping_oracle"],
        "acceptance": {
            **gates,
            "all_checks_pass": summary["all_acceptance_gates_pass"],
        },
        "claim_limits": [
            "M29 currents are retained as measured results, but its BGN off causal label is invalid.",
            "The M29 downstream C1 result uses a mesh-identical M27 input-doping oracle because M29 omitted donor/acceptor exports.",
            "The Bennett/Wilson versus OldSlotboom C1 difference is diagnostic and is not a same-model parity score.",
            "M39 does not change a Vela default and does not rerun Sentaurus or Vela.",
        ],
    }
    args.portable_dir.mkdir(parents=True, exist_ok=True)
    portable_path = args.portable_dir / "m39_bgn_chain_first_divergence_report.json"
    portable_path.write_text(json.dumps(portable, indent=2) + "\n",
                             encoding="utf-8", newline="\n")
    portable_ledger = args.portable_dir / "m39_model_selection_ledger.csv"
    portable_ledger.write_bytes(ledger_path.read_bytes())
    print(f"\n[M39] summary written to {out_json}")
    print(f"[M39] portable report written to {portable_path}")
    print(f"[M39] first divergence: {summary['first_divergence']['stage']} "
          f"({summary['first_divergence']['declared']} -> "
          f"{summary['first_divergence']['observed']})")
    print(f"[M39] acceptance gates: {sum(gates.values())}/{len(gates)}")
    return 0 if summary["all_acceptance_gates_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
