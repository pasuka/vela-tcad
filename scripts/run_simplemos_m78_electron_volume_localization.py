"""Read-only M78 attribution of the existing M65/M74 electron-volume A/B.

No device runner, subprocess, remote access, or state writer is used. Sparse
linear solves are fixed-state Poisson observers, never self-consistent solves.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

import run_simplemos_m73_material_partitioned_poisson_ledger as m73
import run_simplemos_m74_electron_poisson_charge_volume as m74

REPO, ROOT = m74.REPO, m74.ROOT
SCRIPT = Path(__file__).resolve()
SLUG = "m78_electron_volume_localization"
CONTRACT = ROOT / f"simplemos_{SLUG}_contract_v1.json"
FREEZE = ROOT / f"simplemos_{SLUG}_contract_freeze_v1.json"
EVIDENCE = ROOT / f"simplemos_{SLUG}_evidence.json"
OUT = ROOT / "electron_volume_localization"
LOCAL = REPO / "build-release" / SLUG
DOC = REPO / "docs/validation/simplemos_m78_electron_volume_localization_2026-09-05.md"
STAGES = ("equilibrium", "drain", "gate")
FIELDS = ("psi", "phin", "phip", "log10_n", "log10_p")
read_json, write_json = m74.read_json, m74.write_json
read_csv, write_csv = m74.read_csv, m74.write_csv
sha256, portable = m74.sha256, m74.portable


def distribution(values):
    a = np.asarray(values, dtype=float)
    return {"median": float(np.median(a)), "p95": float(np.percentile(a, 95)),
            "maximum": float(np.max(a))}


def norms(values):
    a = np.asarray(values, dtype=float)
    if not len(a):
        raise ValueError("empty spatial support")
    return {"mean": float(np.mean(a)), "l1": float(np.sum(np.abs(a))),
            "rms": float(np.sqrt(np.mean(a*a))), "p95_abs": float(np.percentile(np.abs(a), 95)),
            "max_abs": float(np.max(np.abs(a)))}


def check_hash(path, expected):
    if not path.is_file() or sha256(path) != expected:
        raise ValueError(f"identity gate failed; do not recompute: {path}")


def historical_inputs():
    """Check historical data hashes, excluding intentionally rebuilt executables."""
    inputs = set()
    for ev in (m74.EVIDENCE, m73.EVIDENCE, m73.M72_EVIDENCE, m73.M69_EVIDENCE,
               m74.M65_EVIDENCE):
        evidence = read_json(ev)
        if evidence["status"] != "frozen":
            raise ValueError(f"unqualified evidence: {ev}")
        inputs.add(ev)
        artifacts = evidence.get("artifacts", evidence.get("hashes", evidence.get("result_hashes")))
        if artifacts is None:
            raise ValueError(f"unknown historical evidence schema: {ev}")
        for rel, expected in artifacts.items():
            path = REPO / rel
            check_hash(path, expected)
            inputs.add(path)
    # All 48 baseline states were frozen by M69, including equilibrium.
    f69 = ROOT / "simplemos_m69_stage_residual_localization_contract_freeze_v2.json"
    f73 = m73.FREEZE
    baseline = read_json(f69)["upstream_hashes"]
    inputs.update((f69, f73, m74.CONTRACT, m74.FREEZE, m74.M65_VELA, m74.M65_SENT))
    manifest = read_json(m74.PORTABLE_MANIFEST)
    check_hash(m74.CONTRACT, manifest["contract_sha256"])
    check_hash(m74.RUN_MANIFEST, manifest["source_sha256"])
    inputs.add(m74.RUN_MANIFEST)
    if len(manifest["workflows"]) != 16:
        raise ValueError("expected 16 M74 workflows")
    rows = []
    for workflow in manifest["workflows"]:
        if [p["phase"] for p in workflow["phases"]] != list(m74.PHASES):
            raise ValueError("candidate stage path mismatch")
        for stage, phase in zip(STAGES, workflow["phases"], strict=True):
            for kind in ("state", "curve", "config"):
                path = REPO / phase[kind]
                check_hash(path, phase[f"{kind}_sha256"])
                inputs.add(path)
            source = REPO / phase["source_config"]
            base_state = source.parent / "state.csv"
            check_hash(base_state, baseline[portable(base_state)])
            inputs.update((source, base_state, source.parent / "curve.csv"))
            if m74.normalized_config(read_json(source), False) != m74.normalized_config(
                    read_json(REPO / phase["config"]), True):
                raise ValueError("non-single-axis config difference")
            rows.append({"device": workflow["device"], "drain_voltage_V": workflow["drain_voltage_V"],
                         "stage": stage, "baseline_state": portable(base_state),
                         "baseline_sha256": sha256(base_state), "candidate_state": phase["state"],
                         "candidate_sha256": phase["state_sha256"], "identity_pass": True})
    if m74.baseline_hash_errors():
        raise ValueError("M65 final curve/state mismatch")
    # M73's frozen geometry, material and Sentaurus fields are needed to
    # reconstruct its node-level predictions, which were not saved in CSV.
    for rel, expected in read_json(f73)["upstream_hashes"].items():
        path = REPO / rel
        check_hash(path, expected)
        inputs.add(path)
    # Source configs were frozen before the M74 intervention.
    old = read_json(m74.FREEZE)["upstream_hashes"]
    for path in list(inputs):
        rel = portable(path)
        if rel in old and path.name == "config.json":
            check_hash(path, old[rel])
    inputs.update((m73.SCRIPT, m74.SCRIPT, m73.STAGES, m73.PAIRS,
                   ROOT / "native_gate_reaction/m72_native_gate_reaction_stage_ledger.csv",
                   ROOT / "tight_convergence_port_burst/m60_tight_convergence_port_burst_report.json"))
    return inputs, rows


def freeze():
    if FREEZE.exists() or CONTRACT.exists():
        raise ValueError("M78 contract already exists; never silently refreeze")
    inputs, identities = historical_inputs()
    contract = {
        "schema": "vela.simplemos.m78_electron_volume_localization_contract.v1",
        "objective": "Locate existing M74 electron-only response by stage and space; separate pair improvement from absolute-error deterioration.",
        "matrix": read_json(m74.CONTRACT)["matrix"],
        "execution": {"new_sentaurus_solves": 0, "new_vela_solves": 0,
                      "fixed_state_linear_observers_only": True, "production_changes": False},
        "comparison": "M65 matched-ni/no-BGN vs M74 poisson_electron_transport_node_volume=true; all other physics and bias paths frozen",
        "reference": "M65 no-BGN Sentaurus endpoint currents for like-physics errors. M60 remains the tight-convergence original-physics reference; do not substitute its BGN currents into the no-BGN ledger.",
        "definitions": {
            "delta": "M74 minus M65 at identical node IDs within each device; never subtract different-device node IDs",
            "stage_increments": "equilibrium delta; drain delta minus equilibrium delta; gate delta minus drain delta",
            "pair_stage_share": "absolute high-minus-low control-node psi (or psi-phin) stage increment divided by sum of absolute increments",
            "control": "M69 frozen Sentaurus gate control node, used unchanged at all three stages",
            "norms": "unweighted node mean, L1, RMS, P95 absolute and maximum; carrier fields restricted to Si; volume-weighted psi L1 also reported",
            "spatial_partition": "Si nodes at depth x - minimum Si/SiO2 interface x > 0.1 um are substrate; remaining y<-0.125 source, y>0.125 drain, others channel; exclusive, exhaustive",
            "interface": "ALL nodes belonging to both Si and SiO2; gate footprint subset abs(y)<=0.125 um reported separately, no restricted profile called full interface",
            "electron_reaction": "Poisson residual sign +q*n*volume in C/m; forcing q*n_baseline*(V_Si-V_all); relaxation q*(n_candidate-n_baseline)*V_Si",
            "gate_reaction": "unconstrained legacy Poisson residual at gate nodes, qualified against existing M72 native reaction before use",
            "poisson_prediction": "M73 raw candidate-minus-legacy fixed-state correction, separately total and electron, using frozen Sentaurus states. Also report electron-only baseline-state forcing response -K_legacy^-1*q*n_baseline*(V_Si-V_all)",
            "pair_prediction": "negative high-minus-low gate-minus-drain delta psi divided by Vt*ln(10); compare both actual electrostatic proxy and M74 actual current reduction",
            "localization": "report forcing support separately from propagated psi support; interface L1 share has all Si nodes as denominator",
            "next_step_gate": "M79 local-support A/B only if at least 6/8 pairs have same-sign current reduction and BOTH devices have >=80% gate-increment psi L1 on Si/SiO2 mixed nodes, stable under 0.05/0.10/0.20 um depth partitions",
        },
        "acceptance": {"cases": 16, "stages": 48, "pairs": 8, "hash_errors": 0,
                       "bias_tolerance_V": 1e-12, "native_gate_replay_tolerance_C_per_m": 1e-20,
                       "poisson_replay_tolerance_V": 1e-7, "m73_ledger_replay_tolerance_V": 1e-10,
                       "balance_tolerance_V": 2e-7, "current_ledger_tolerance_dex": 1e-12},
        "guards": ["No retuning closed topics", "No new self-consistent solves", "Stop on state hash failure",
                   "Pair improvement is not absolute accuracy improvement", "Linear prediction is not a validated correction",
                   "No production-default changes", "No new or changed remote decks"],
    }
    write_json(CONTRACT, contract)
    write_json(FREEZE, {"schema": "vela.simplemos.m78_contract_freeze.v1", "status": "frozen_before_analysis",
                       "contract_sha256": sha256(CONTRACT), "baseline_stage_count": len(identities),
                       "candidate_stage_count": len(identities),
                       "input_hashes": {portable(p): sha256(p) for p in sorted(inputs)}})


def validate():
    frozen = read_json(FREEZE)
    if frozen["status"] != "frozen_before_analysis":
        raise ValueError("M78 contract not frozen")
    check_hash(CONTRACT, frozen["contract_sha256"])
    for rel, expected in frozen["input_hashes"].items():
        check_hash(REPO / rel, expected)
    return read_json(CONTRACT)


def supports(device, geo, depth=0.1):
    mesh = read_json(m73.M8 / "vela" / device / "mesh.json")
    regions = {int(r["id"]): r["material"] for r in mesh["regions"]}
    materials = [set() for _ in range(geo.count)]
    for tri in mesh["triangles"]:
        for node in tri["node_ids"]:
            materials[node].add(regions[tri["region_id"]])
    si = np.array(["Si" in r for r in materials])
    interface = np.array([{"Si", "SiO2"} <= r for r in materials])
    xy = np.array([geo.coords[i] for i in range(geo.count)])
    surface = float(np.min(xy[interface, 0]))
    deep = si & (xy[:, 0]-surface > depth)
    masks = {"source": si & ~deep & (xy[:, 1] < -0.125),
             "channel": si & ~deep & (np.abs(xy[:, 1]) <= 0.125),
             "drain": si & ~deep & (xy[:, 1] > 0.125), "substrate": deep}
    if not np.array_equal(sum(m.astype(int) for m in masks.values()), si.astype(int)):
        raise ValueError("spatial partition is not exhaustive and exclusive")
    masks.update({"all_si": si, "interface": interface,
                  "si_nitride_interface": np.array([{"Si", "Nitride"} <= r for r in materials]),
                  "gate_interface": interface & (np.abs(xy[:, 1]) <= 0.125+1e-12)})
    gate = np.array(next(c["node_ids"] for c in mesh["contacts"] if c["name"] == "gate"), dtype=int)
    return masks, gate, xy


def arrays(path, count):
    rows = read_csv(path)
    if len(rows) != count or {int(r["node_id"]) for r in rows} != set(range(count)):
        raise ValueError(f"invalid node support: {path}")
    rows.sort(key=lambda r: int(r["node_id"]))
    result = {f: np.array([float(r[f]) for r in rows])
              for f in ("psi", "phin", "phip", "electrons_m3", "holes_m3")}
    if not all(np.all(np.isfinite(v)) for v in result.values()):
        raise ValueError("nonfinite state")
    return result


def analyze():
    contract = validate()
    _, identities = historical_inputs()
    write_csv(OUT / "m78_identity_ledger.csv", identities)
    workflows = read_json(m74.PORTABLE_MANIFEST)["workflows"]
    controls, exports = m73.control_nodes(), m73.exports()
    native = {(r["device"], float(r["drain_voltage_V"]), r["stage"]): float(r["vela_native_gate_charge_C_per_m"])
              for r in read_csv(ROOT / "native_gate_reaction/m72_native_gate_reaction_stage_ledger.csv")}
    old73 = {(r["device"], float(r["drain_voltage_V"]), r["stage"], r["variant"]): r
             for r in read_csv(m73.STAGES)}
    cache, data = {}, {}
    stages, spatial, nodes, gates, predictions, sensitivity = [], [], [], [], [], []
    checks = {"maximum_native_gate_replay_error_C_per_m": 0.0,
              "maximum_baseline_poisson_replay_V": 0.0, "maximum_candidate_poisson_replay_V": 0.0,
              "maximum_m73_control_replay_error_V": 0.0, "maximum_charge_response_balance_error_V": 0.0,
              "maximum_bias_error_V": 0.0, "maximum_boltzmann_delta_identity_dex": 0.0}
    for workflow in workflows:
        device, vd = workflow["device"], float(workflow["drain_voltage_V"])
        if device not in cache:
            cache[device] = m73.Geometry(device)
        geo = cache[device]
        masks, gate_ids, xy = supports(device, geo)
        si, interface = masks["all_si"], masks["interface"]
        control = controls[(device, vd)]
        va, vs = geo.volumes["all_cell"], geo.volumes["barycentric_si"]
        net = m73.doping(device, geo.count)
        previous = None
        for stage, phase in zip(STAGES, workflow["phases"], strict=True):
            a = arrays((REPO / phase["source_config"]).parent / "state.csv", geo.count)
            b = arrays(REPO / phase["state"], geo.count)
            ac = read_csv((REPO / phase["source_config"]).parent / "curve.csv")
            bc = read_csv(REPO / phase["curve"])
            for curve in (ac, bc):
                if not curve or any(r["converged"].lower() not in ("1", "true") for r in curve):
                    raise ValueError("unconverged phase curve")
            expected_bias = float(workflow["gate_voltage_V"]) if stage == "gate" else (vd if stage == "drain" else 0.0)
            bias_error = max(abs(float(c[-1]["bias_V"])-expected_bias) for c in (ac, bc))
            checks["maximum_bias_error_V"] = max(checks["maximum_bias_error_V"], bias_error)
            delta = {f: b[f]-a[f] for f in ("psi", "phin", "phip")}
            for out, field in (("log10_n", "electrons_m3"), ("log10_p", "holes_m3")):
                if np.any(a[field][si] <= 0) or np.any(b[field][si] <= 0):
                    raise ValueError("nonpositive Si density")
                delta[out] = np.zeros(geo.count)
                delta[out][si] = np.log10(b[field][si])-np.log10(a[field][si])
            checks["maximum_boltzmann_delta_identity_dex"] = max(checks["maximum_boltzmann_delta_identity_dex"],
                float(np.max(np.abs(delta["log10_n"][si]-(delta["psi"][si]-delta["phin"][si])/m73.VT_LN10))))
            forcing = m73.Q*a["electrons_m3"]*(vs-va)
            relax_n = m73.Q*(b["electrons_m3"]-a["electrons_m3"])*vs
            relax_p = -m73.Q*(b["holes_m3"]-a["holes_m3"])*va
            total_e = forcing+relax_n
            response = {"forcing": geo.solve("legacy", forcing),
                        "electron_relaxation": geo.solve("legacy", relax_n),
                        "hole_relaxation": geo.solve("legacy", relax_p)}
            residual_a = geo.matrices["legacy"]@a["psi"] + m73.Q*(a["electrons_m3"]-a["holes_m3"]-net)*va
            residual_b = geo.matrices["legacy"]@b["psi"] + m73.Q*(b["electrons_m3"]*vs-(b["holes_m3"]+net)*va)
            for label, residual in (("baseline", residual_a), ("candidate", residual_b)):
                name = f"maximum_{label}_poisson_replay_V"
                checks[name] = max(checks[name], float(np.max(np.abs(geo.solve("legacy", residual)))))
            balance = delta["psi"]-sum(response.values())
            checks["maximum_charge_response_balance_error_V"] = max(checks["maximum_charge_response_balance_error_V"], float(np.max(np.abs(balance))))
            qa, qb = float(np.sum(residual_a[gate_ids])), float(np.sum(residual_b[gate_ids]))
            if stage != "equilibrium":
                checks["maximum_native_gate_replay_error_C_per_m"] = max(checks["maximum_native_gate_replay_error_C_per_m"], abs(qa-native[(device, vd, stage)]))
            increment = delta if previous is None else {f: delta[f]-previous["delta"][f] for f in FIELDS}
            forcing_increment = forcing if previous is None else forcing-previous["forcing"]
            psi_l1 = float(np.sum(np.abs(increment["psi"][si])))
            interface_share = float(np.sum(np.abs(increment["psi"][interface])))/max(psi_l1, 1e-300)
            top_count = max(1, math.ceil(int(np.sum(si))*0.1))
            top_share = float(np.sum(np.sort(np.abs(increment["psi"][si]))[-top_count:]))/max(psi_l1, 1e-300)
            key = {"device": device, "drain_voltage_V": vd, "stage": stage,
                   "gate_voltage_V": float(workflow["gate_voltage_V"]) if stage == "gate" else 0.0}
            pred = {}
            if stage != "equilibrium":
                sp, sn, sh, _ = m73.sentaurus_state(exports[(device, vd, stage)], geo)
                corr = {}
                for variant, (matrix, volume) in {k: m73.VARIANTS[k] for k in ("legacy_all_cell", "region_local_barycentric_si")}.items():
                    vol = geo.volumes[volume]
                    corr[variant] = {
                        "total": geo.solve(matrix, geo.matrices[matrix]@sp+m73.Q*(sn-sh-net)*vol),
                        "electron": geo.solve(matrix, m73.Q*sn*vol)}
                    old = old73[(device, vd, stage, variant)]
                    for component in ("total", "electron"):
                        checks["maximum_m73_control_replay_error_V"] = max(checks["maximum_m73_control_replay_error_V"],
                            abs(float(corr[variant][component][control])-float(old[f"{component}_response_V"])))
                pred = {name: corr["region_local_barycentric_si"][name]-corr["legacy_all_cell"][name] for name in ("total", "electron")}
                for name, prediction in {**pred, "baseline_electron_forcing": response["forcing"]}.items():
                    error = prediction-delta["psi"]
                    predictions.append({**key, "prediction": name,
                        "control_predicted_V": float(prediction[control]), "control_actual_V": float(delta["psi"][control]),
                        "control_error_V": float(error[control]),
                        **{f"si_error_{k}_V": v for k, v in norms(error[si]).items()},
                        "si_relative_l2_error": float(np.linalg.norm(error[si])/max(np.linalg.norm(delta["psi"][si]), 1e-300))})
            for scope, fields in (("state_delta", delta), ("stage_increment", increment)):
                for region, mask in masks.items():
                    for field in FIELDS:
                        spatial.append({**key, "scope": scope, "region": region, "field": field,
                                        "node_count": int(np.sum(mask)), **norms(fields[field][mask])})
            row = {**key, "control_node": control,
                   **{f"control_delta_{f}": float(delta[f][control]) for f in FIELDS},
                   **{f"control_increment_{f}": float(increment[f][control]) for f in FIELDS},
                   "delta_gate_charge_C_per_m": qb-qa, "baseline_gate_charge_C_per_m": qa,
                   "candidate_gate_charge_C_per_m": qb, "si_node_count": int(np.sum(si)),
                   "interface_node_count": int(np.sum(interface)), "gate_increment_interface_psi_l1_share": interface_share,
                   "increment_top10pct_si_psi_l1_share": top_share,
                   "increment_interface_volume_weighted_psi_l1_share": float(np.sum(np.abs(increment["psi"][interface])*vs[interface])/max(np.sum(np.abs(increment["psi"][si])*vs[si]), 1e-300)),
                   "forcing_interface_l1_share": float(np.sum(np.abs(forcing[interface]))/max(np.sum(np.abs(forcing[si])), 1e-300)),
                   "forcing_si_nitride_l1_share": float(np.sum(np.abs(forcing[masks["si_nitride_interface"]]))/max(np.sum(np.abs(forcing[si])), 1e-300)),
                   **{f"control_{name}_response_V": float(value[control]) for name, value in response.items()},
                   **{f"control_increment_{name}_response_V": float(value[control])-(previous["row"][f"control_{name}_response_V"] if previous else 0.0)
                      for name, value in response.items()},
                   "forcing_increment_interface_l1_share": float(np.sum(np.abs(forcing_increment[interface]))/max(np.sum(np.abs(forcing_increment[si])), 1e-300))}
            stages.append(row)
            for node in range(geo.count):
                region = next((r for r in ("source", "channel", "drain", "substrate") if masks[r][node]), "insulator")
                nodes.append({**key, "node_id": node, "x_um": xy[node, 0], "y_um": xy[node, 1], "region": region,
                    "si_sio2_interface": bool(interface[node]), "si_nitride_interface": bool(masks["si_nitride_interface"][node]),
                    "gate_footprint_interface": bool(masks["gate_interface"][node]),
                    **{f"delta_{f}": float(delta[f][node]) if si[node] or f == "psi" else "" for f in FIELDS},
                    "increment_psi_V": increment["psi"][node], "electron_forcing_C_per_m": forcing[node],
                    "electron_relaxation_C_per_m": relax_n[node], "electron_total_reaction_delta_C_per_m": total_e[node],
                    "hole_relaxation_C_per_m": relax_p[node],
                    "m73_total_prediction_V": float(pred["total"][node]) if pred else "",
                    "m73_electron_prediction_V": float(pred["electron"][node]) if pred else "",
                    "m73_total_prediction_error_V": float(pred["total"][node]-delta["psi"][node]) if pred else "",
                    "m73_electron_prediction_error_V": float(pred["electron"][node]-delta["psi"][node]) if pred else ""})
            for node in gate_ids:
                gates.append({**key, "node_id": int(node), "x_um": xy[node, 0], "y_um": xy[node, 1],
                              "baseline_reaction_C_per_m": residual_a[node], "candidate_reaction_C_per_m": residual_b[node],
                              "delta_reaction_C_per_m": residual_b[node]-residual_a[node]})
            for depth in (0.05, 0.1, 0.2):
                varied, _, _ = supports(device, geo, depth)
                for region in ("source", "channel", "drain", "substrate"):
                    mask = varied[region]
                    sensitivity.append({**key, "depth_um": depth, "region": region,
                        "increment_psi_l1_share": float(np.sum(np.abs(increment["psi"][mask]))/max(psi_l1, 1e-300)),
                        "forcing_increment_l1_C_per_m": float(np.sum(np.abs(forcing_increment[mask]))),
                        "control_forcing_increment_response_V": float(geo.solve("legacy", np.where(mask, forcing_increment, 0.0))[control])})
            data[(device, vd, stage)] = {"row": row, "delta": delta, "forcing": forcing,
                                         "prediction": pred, "masks": masks}
            previous = data[(device, vd, stage)]
    return finish(contract, checks, data, stages, spatial, nodes, gates, predictions, sensitivity)


def finish(contract, checks, data, stages, spatial, nodes, gates, predictions, sensitivity):
    cases = read_csv(m74.CASES)
    case_map = {(r["device"], float(r["drain_voltage_V"])): r for r in cases}
    old_pairs = {(r["low_device"], r["high_device"], float(r["drain_voltage_V"])): r for r in read_csv(m74.PAIRS)}
    pairs, paired_fields = [], []
    # Pair norm differences refer to scalar physical-region summaries, not
    # correspondence between node IDs of distinct device meshes.
    region_index = {(r["device"], r["drain_voltage_V"], r["stage"], r["scope"], r["region"], r["field"]): r for r in spatial}
    for pair in contract["matrix"]["pairs"]:
        low, high, vd = pair["low_device"], pair["high_device"], float(pair["drain_voltage_V"])
        key = {"low_device": low, "high_device": high, "drain_voltage_V": vd}
        lo, hi = case_map[(low, vd)], case_map[(high, vd)]
        sl, sh = [float(r["endpoint_candidate_over_baseline_log_shift_dex"]) for r in (lo, hi)]
        base_growth = float(hi["baseline_error_dex"])-float(lo["baseline_error_dex"])
        candidate_growth = float(hi["candidate_error_dex"])-float(lo["candidate_error_dex"])
        p = {**key, "gate_voltage_V": pair["diagnostic_gate_voltage_V"], "baseline_pair_growth_dex": base_growth,
             "candidate_pair_growth_dex": candidate_growth, "pair_reduction_dex": sl-sh,
             "low_current_shift_dex": sl, "high_current_shift_dex": sh,
             "low_absolute_error_worsening_dex": -float(lo["absolute_error_improvement_dex"]),
             "high_absolute_error_worsening_dex": -float(hi["absolute_error_improvement_dex"]),
             "pair_accounting_identity_error_dex": (base_growth-candidate_growth)-(sl-sh)}
        for field in ("psi", "barrier"):
            inc = []
            for stage in STAGES:
                h, l = [data[(d, vd, stage)]["row"] for d in (high, low)]
                value = h["control_increment_psi"]-l["control_increment_psi"]
                if field == "barrier":
                    value -= h["control_increment_phin"]-l["control_increment_phin"]
                inc.append(value)
                p[f"pair_{stage}_{field}_increment_V"] = value
            denominator = max(sum(abs(v) for v in inc), 1e-300)
            for stage, value in zip(STAGES, inc, strict=True):
                p[f"pair_{stage}_{field}_absolute_share"] = abs(value)/denominator
        for name in ("total", "electron"):
            changes = []
            for device in (low, high):
                control = data[(device, vd, "gate")]["row"]["control_node"]
                changes.append(float(data[(device, vd, "gate")]["prediction"][name][control]-data[(device, vd, "drain")]["prediction"][name][control]))
            predicted = -(changes[1]-changes[0])/m73.VT_LN10
            actual_psi = -p["pair_gate_psi_increment_V"]/m73.VT_LN10
            p[f"m73_{name}_predicted_reduction_dex"] = predicted
            p[f"m73_{name}_vs_actual_psi_error_dex"] = predicted-actual_psi
            p[f"m73_{name}_vs_current_error_dex"] = predicted-(sl-sh)
        p["actual_gate_psi_reduction_proxy_dex"] = -p["pair_gate_psi_increment_V"]/m73.VT_LN10
        p["actual_gate_barrier_reduction_proxy_dex"] = -p["pair_gate_barrier_increment_V"]/m73.VT_LN10
        p["actual_gate_phin_contribution_dex"] = p["actual_gate_barrier_reduction_proxy_dex"]-p["actual_gate_psi_reduction_proxy_dex"]
        for component in ("forcing", "electron_relaxation", "hole_relaxation"):
            h, l = [data[(d, vd, "gate")]["row"] for d in (high, low)]
            p[f"pair_gate_{component}_reduction_proxy_dex"] = -(h[f"control_increment_{component}_response_V"]-l[f"control_increment_{component}_response_V"])/m73.VT_LN10
        p["m73_electron_state_operator_difference_dex"] = p["m73_electron_predicted_reduction_dex"]-p["pair_gate_forcing_reduction_proxy_dex"]
        p["m73_electron_missing_feedback_difference_dex"] = p["pair_gate_forcing_reduction_proxy_dex"]-p["actual_gate_psi_reduction_proxy_dex"]
        error = p["m73_electron_vs_actual_psi_error_dex"]
        p["m73_electron_missing_feedback_error_fraction"] = p["m73_electron_missing_feedback_difference_dex"]/error if abs(error) > 1e-15 else 0.0
        p["m74_frozen_prediction_difference_dex"] = p["m73_total_predicted_reduction_dex"]-float(old_pairs[(low, high, vd)]["m73_predicted_pair_reduction_dex"])
        for label, device in (("low", low), ("high", high)):
            r = data[(device, vd, "gate")]["row"]
            p[f"{label}_interface_psi_l1_share"] = r["gate_increment_interface_psi_l1_share"]
            p[f"{label}_forcing_interface_l1_share"] = r["forcing_interface_l1_share"]
            p[f"{label}_forcing_increment_interface_l1_share"] = r["forcing_increment_interface_l1_share"]
            p[f"{label}_interface_volume_weighted_psi_l1_share"] = r["increment_interface_volume_weighted_psi_l1_share"]
            p[f"{label}_control_gate_delta_psi_V"] = r["control_delta_psi"]
            p[f"{label}_control_gate_increment_psi_V"] = r["control_increment_psi"]
            p[f"{label}_gate_charge_increment_response_C_per_m"] = r["delta_gate_charge_C_per_m"]-data[(device, vd, "drain")]["row"]["delta_gate_charge_C_per_m"]
        p["local_mixed_node_response_qualified"] = min(p["low_interface_psi_l1_share"], p["high_interface_psi_l1_share"]) >= 0.8
        pairs.append(p)
        for stage in STAGES:
            for scope in ("state_delta", "stage_increment"):
                for region in ("source", "channel", "drain", "substrate", "all_si", "interface", "gate_interface"):
                    for field in FIELDS:
                        a, b = [region_index[(d, vd, stage, scope, region, field)] for d in (low, high)]
                        paired_fields.append({**key, "stage": stage, "scope": scope, "region": region, "field": field,
                            **{f"high_minus_low_{n}": b[n]-a[n] for n in ("mean", "l1", "rms", "p95_abs", "max_abs")}})
    acceptance = {
        "case_and_stage_count": len(stages) == 48 and len(cases) == 16,
        "pair_count": len(pairs) == 8, "historical_hash_identity": True, "single_axis_config_identity": True,
        "bias_identity": checks["maximum_bias_error_V"] <= contract["acceptance"]["bias_tolerance_V"],
        "native_gate_replay": checks["maximum_native_gate_replay_error_C_per_m"] <= 1e-20,
        "baseline_poisson_replay": checks["maximum_baseline_poisson_replay_V"] <= 1e-7,
        "candidate_poisson_replay": checks["maximum_candidate_poisson_replay_V"] <= 1e-7,
        "charge_response_balance": checks["maximum_charge_response_balance_error_V"] <= 2e-7,
        "m73_replay": checks["maximum_m73_control_replay_error_V"] <= 1e-10,
        "pair_accounting": max(abs(r["pair_accounting_identity_error_dex"]) for r in pairs) <= 1e-12,
        "m74_prediction_replay": max(abs(r["m74_frozen_prediction_difference_dex"]) for r in pairs) <= 1e-10,
        "pair_charge_response_identity": max(abs(sum(r[f"pair_gate_{c}_reduction_proxy_dex"] for c in ("forcing", "electron_relaxation", "hole_relaxation"))-r["actual_gate_psi_reduction_proxy_dex"]) for r in pairs) <= 2e-7/m73.VT_LN10,
    }
    acceptance["all_checks_pass"] = all(acceptance.values())
    summaries = {name: distribution([p[name] for p in pairs]) for name in (
        "low_current_shift_dex", "high_current_shift_dex", "pair_reduction_dex",
        "baseline_pair_growth_dex", "candidate_pair_growth_dex", "pair_gate_psi_absolute_share", "pair_gate_barrier_absolute_share",
        "low_interface_psi_l1_share", "high_interface_psi_l1_share",
        "low_forcing_interface_l1_share", "high_forcing_interface_l1_share",
        "low_forcing_increment_interface_l1_share", "high_forcing_increment_interface_l1_share",
        "low_interface_volume_weighted_psi_l1_share", "high_interface_volume_weighted_psi_l1_share",
        "actual_gate_psi_reduction_proxy_dex", "actual_gate_barrier_reduction_proxy_dex", "actual_gate_phin_contribution_dex",
        "m73_electron_state_operator_difference_dex", "m73_electron_missing_feedback_difference_dex", "m73_electron_missing_feedback_error_fraction",
        "pair_gate_forcing_reduction_proxy_dex", "pair_gate_electron_relaxation_reduction_proxy_dex", "pair_gate_hole_relaxation_reduction_proxy_dex")}
    for label, devices in (("low", {"n17", "n18", "n19", "n20"}), ("high", {"n21", "n22", "n23", "n24"})):
        for region in ("source", "channel", "drain", "substrate"):
            rows = [r for r in sensitivity if r["device"] in devices and r["stage"] == "gate" and r["depth_um"] == 0.1 and r["region"] == region]
            summaries[f"{label}_{region}_gate_increment_psi_l1_share"] = distribution([r["increment_psi_l1_share"] for r in rows])
    for name in ("total", "electron"):
        for comparison in ("actual_psi", "current"):
            column = f"m73_{name}_vs_{comparison}_error_dex"
            summaries[f"absolute_{column}"] = distribution([abs(p[column]) for p in pairs])
    summaries["same_sign_improvement_pair_count"] = sum(p["pair_reduction_dex"] > 0 for p in pairs)
    summaries["localized_mixed_node_pair_count"] = sum(p["local_mixed_node_response_qualified"] for p in pairs)
    summaries["all_case_absolute_errors_worsen"] = all(float(r["absolute_error_improvement_dex"]) < 0 for r in cases)
    summaries["gate_si_nitride_forcing_l1_share"] = distribution([r["forcing_si_nitride_l1_share"] for r in stages if r["stage"] == "gate"])
    summaries["maximum_control_barrier_proxy_current_error_dex"] = max(abs(p["actual_gate_barrier_reduction_proxy_dex"]-p["pair_reduction_dex"]) for p in pairs)
    proceed = summaries["same_sign_improvement_pair_count"] >= 6 and summaries["localized_mixed_node_pair_count"] >= 6
    classification = ("gate_response_localized_m79_candidate" if proceed else "gate_response_requires_distributed_support_no_m79_local_ab") if acceptance["all_checks_pass"] else "observer_or_identity_failure"
    report = {"schema": "vela.simplemos.m78_report.v1", "status": "qualified" if acceptance["all_checks_pass"] else "failed",
              "classification": classification,
              "execution": contract["execution"], "acceptance": acceptance, "checks": checks, "summary": summaries,
              "interpretation": "Pair reduction equals low-current increase minus high-current increase. All endpoint errors increase; geometric forcing localization does not prove localized self-consistent potential response.",
              "m79_local_support_ab_qualified": proceed and acceptance["all_checks_pass"]}
    outputs = {"m78_stage_ledger.csv": stages, "m78_spatial_field_ledger.csv": spatial,
               "m78_pair_field_ledger.csv": paired_fields, "m78_pair_ledger.csv": pairs,
               "m78_interface_node_ledger.csv": [r for r in nodes if r["si_sio2_interface"]],
               "m78_gate_node_reaction_ledger.csv": gates, "m78_prediction_ledger.csv": predictions,
               "m78_partition_sensitivity_ledger.csv": sensitivity}
    for name, rows in outputs.items():
        write_csv(OUT / name, rows)
    write_csv(LOCAL / "m78_all_node_ledger.csv", nodes)
    write_json(OUT / "m78_report.json", report)
    lines = ["# SimpleMOS M78 电子电荷体积响应的阶段与空间定位", "",
             "M78 只读复用 M65/M74 的 48 对阶段状态；没有新增器件求解或修改生产默认值。", "",
             "配对差的改善来自低 NWell 电流上升更多，16 个诊断端点的绝对误差全部恶化。",
             "no-BGN 绝对误差沿用 M65 同物理 Sentaurus 参考；M60 保留为原物理收紧收敛参考，不混入此账本。", "",
             "| 配对 | Vd | 低/高电流增量 (dex) | 配对差降幅 (dex) | gate 电势份额 | 低/高界面响应份额 |",
             "|---|---:|---:|---:|---:|---:|"]
    for p in pairs:
        lines.append(f"| {p['low_device']}/{p['high_device']} | {p['drain_voltage_V']:g} | {p['low_current_shift_dex']:.6f} / {p['high_current_shift_dex']:.6f} | {p['pair_reduction_dex']:.6f} | {p['pair_gate_psi_absolute_share']:.2%} | {p['low_interface_psi_l1_share']:.2%} / {p['high_interface_psi_l1_share']:.2%} |")
    lines += ["", "## 新的定位结论", "",
              f"1. gate 阶段承担配对电势响应的中位 {summaries['pair_gate_psi_absolute_share']['median']:.4%}。这是电子体积干预响应的阶段份额，与 M69 原始跨求解器残差的约 75.20% 份额属于不同问题。",
              f"2. 低/高 NWell 的电流上升中位分别为 {summaries['low_current_shift_dex']['median']:.6f} / {summaries['high_current_shift_dex']['median']:.6f} dex。配对增长中位由 {summaries['baseline_pair_growth_dex']['median']:.6f} 降至 {summaries['candidate_pair_growth_dex']['median']:.6f} dex，但这来自低 NWell 误差追近高 NWell，不能作为精度改进。",
              "3. 电势响应不局限于少数 Si/SiO2 混合节点。下表给出 gate 增量响应在 0.10 um 深度分区中的份额；分区互斥且覆盖全部 Si。0.05/0.10/0.20 um 深度敏感性账本全部保存，主响应仍在浅沟道和源漏邻域。", "",
              "| 空间区域 | 低 NWell 中位份额 | 高 NWell 中位份额 |", "|---|---:|---:|"]
    for region in ("source", "channel", "drain", "substrate"):
        lines.append(f"| {region} | {summaries[f'low_{region}_gate_increment_psi_l1_share']['median']:.2%} | {summaries[f'high_{region}_gate_increment_psi_l1_share']['median']:.2%} |")
    lines += ["", "各列是分别取中位数，不能要求四个中位数严格相加为 100%。", "",
              f"4. 无载流子反馈的直接电子体积项，其配对 gate 响应代理中位为 {summaries['pair_gate_forcing_reduction_proxy_dex']['median']:.6f} dex；电子浓度自洽调整贡献 {summaries['pair_gate_electron_relaxation_reduction_proxy_dex']['median']:.6f} dex，空穴调整贡献 {summaries['pair_gate_hole_relaxation_reduction_proxy_dex']['median']:.6f} dex。电子调整抵消了很大一部分直接驱动。分项在每个工况精确加和，不能将各项中位数当作加和恒等式。这里的空穴项是电子单轴干预引发的被动响应，不是 M75 空穴体积开关。",
              "", "| 配对 | Vd | 直接电子体积项 | 电子自洽调整 | 空穴自洽调整 | 实际电势代理 | 实际电流差降幅 |",
              "|---|---:|---:|---:|---:|---:|---:|"]
    for p in pairs:
        lines.append(f"| {p['low_device']}/{p['high_device']} | {p['drain_voltage_V']:g} | {p['pair_gate_forcing_reduction_proxy_dex']:.6f} | {p['pair_gate_electron_relaxation_reduction_proxy_dex']:.6f} | {p['pair_gate_hole_relaxation_reduction_proxy_dex']:.6f} | {p['actual_gate_psi_reduction_proxy_dex']:.6f} | {p['pair_reduction_dex']:.6f} |")
    lines += ["", "上表代理单位均为 dex（电势除以 Vt ln 10）。电势代理与电流差仍不完全相等；本任务未声称电流连续性闭合。", "",
              "5. M73 逐节点线性观察器已从冻结输入重建，并与原控制点账本核对。M73 total 预测对实际电势配对响应的绝对误差中位/P95/最大如下；electron-only 预测仍有相近误差，因此不能只靠去掉 total 中的其它分项来解决预测偏强。", "",
              "| 统计量 | 中位 | P95 | 最大 |", "|---|---:|---:|---:|"]
    for name in ("absolute_m73_total_vs_actual_psi_error_dex", "absolute_m73_electron_vs_actual_psi_error_dex",
                 "absolute_m73_total_vs_current_error_dex", "low_interface_psi_l1_share", "high_interface_psi_l1_share",
                 "low_interface_volume_weighted_psi_l1_share", "high_interface_volume_weighted_psi_l1_share"):
        v = summaries[name]
        lines.append(f"| {name} | {v['median']:.7g} | {v['p95']:.7g} | {v['maximum']:.7g} |")
    lines += ["", "界面份额行单位为比例，其余误差行单位为 dex。节点级预测误差保存在完整节点 CSV，Si 区范数保存在 prediction CSV。", "",
              f"M73 electron 预测与实际电势响应的差，可逐对拆为“冻结 Sentaurus 状态/算子与 M65 基态驱动之差”及“直接驱动减去自洽反馈后响应之差”。后一项占该预测误差的中位 {summaries['m73_electron_missing_feedback_error_fraction']['median']:.2%}；这是 M73→M74 的响应误差记账，不是对原始跨求解器 Id–Vg 偏差的闭合率。", "",
              "6. M69 固定单点的 psi-phin 代理不能作为所有控制工况的电流解释：n19/n23、Vd=1 V 甚至给出相反方向；n20/n24、Vd=1 V 显著放大。电势分量与直接电流账本需单独保留，不据此重开已闭环的准费米实现排查。", "",
              "## 范围与下一步边界", "",
              "本次核验的是 16 个已冻结诊断端点及其三阶段路径，不是重新求解或重验 816 点生产 Id–Vg 全矩阵。",
              "完整材料几何还含 Si/Nitride 接面。电子体积开关作用于所有含绝缘体体积的 Si 节点，并不限于 Si/SiO2；完整节点账本增加了 Nitride 支撑标签。直接电荷总量与 gate 增量、以及经 Poisson 传播后的电势响应须分开解释。",
              "按预冻结的 80% 混合节点响应门槛，0/8 配对满足局部化条件。因此本轮停止“少数 Si/SiO2 节点局部 A/B”路线，不启动 M79 器件求解，不启用 M74–M77 生产开关。",
              "高 NWell 的绝对 Id–Vg 偏差仍未闭环。后续若另开任务，应先研究已保存电荷反馈及分布式电势响应的可比性，重新冻结观察器合同，不能把本次配对差缩小直接升级为修复。", "",
              "## 验收与复现", "",
              f"M65/M74 各 48 份状态的历史哈希通过；配置唯一比较轴、阶段偏压及 8 配对身份通过。Poisson 分项响应最大加和误差 {checks['maximum_charge_response_balance_error_V']:.3e} V，原生栅反力回放最大误差 {checks['maximum_native_gate_replay_error_C_per_m']:.3e} C/m，M73 控制点回放最大误差 {checks['maximum_m73_control_replay_error_V']:.3e} V。"]
    lines += ["", f"分类：`{report['classification']}`；验收：`{acceptance['all_checks_pass']}`。",
              "", "界面份额为完整 Si/SiO2 混合节点上的 gate 增量 |delta psi| 之和 / 全部 Si 节点相同量之和。",
              "电子电荷反力采用 Poisson 残差正号 +q*n*V，单位 C/m；它不是端口电流。",
              "M73 total 与 electron 预测分开保存；原 total 候选包含空穴/掺杂/介电分项，不能当作电子单轴的定量预测。",
              "", "复核：`python scripts/run_simplemos_m78_electron_volume_localization.py --verify`。",
              "完整节点账本位于 ignored build-release/m78_electron_volume_localization；界面、栅反力和各级汇总位于 reference_tcad/simplemos_sentaurus2022/electron_volume_localization。", ""]
    DOC.write_text("\n".join(lines), encoding="utf-8", newline="\n")
    validate()
    artifacts = [OUT / name for name in outputs] + [OUT / "m78_identity_ledger.csv", OUT / "m78_report.json", DOC, LOCAL / "m78_all_node_ledger.csv"]
    write_json(EVIDENCE, {"schema": "vela.simplemos.m78_evidence.v1", "status": report["status"],
                         "contract_sha256": sha256(CONTRACT), "implementation_sha256": sha256(SCRIPT),
                         "test_sha256": sha256(REPO / "tests/regression/test_simplemos_m78_electron_volume_localization.py"),
                         "artifacts": {portable(p): sha256(p) for p in artifacts}})
    return report


def verify():
    validate()
    evidence = read_json(EVIDENCE)
    check_hash(CONTRACT, evidence["contract_sha256"])
    check_hash(SCRIPT, evidence["implementation_sha256"])
    check_hash(REPO / "tests/regression/test_simplemos_m78_electron_volume_localization.py", evidence["test_sha256"])
    for rel, expected in evidence["artifacts"].items():
        check_hash(REPO / rel, expected)
    report = read_json(OUT / "m78_report.json")
    if not report["acceptance"]["all_checks_pass"]:
        raise ValueError(f"M78 observer qualification failed: {report['acceptance']}")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--freeze-contract", action="store_true")
    group.add_argument("--analyze", action="store_true")
    group.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.freeze_contract:
        freeze()
        print("M78 frozen before analysis; 48 baseline + 48 candidate state hashes verified")
    else:
        import json
        report = verify() if args.verify else analyze()
        print(json.dumps(report, indent=2))
        if not report["acceptance"]["all_checks_pass"]:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
