"""P1 read-only replay. Unit-valued references are not accuracy tolerances.

No metric produced here is used by a solver. Requires h5py/numpy on the VM.
"""
from pathlib import Path
from fractions import Fraction
from decimal import Decimal, localcontext
import argparse
import csv
import hashlib
import json
import math


def missing(reason):
    return {"status": "not_evaluated", "reason": reason}


def metric(value, **extra):
    if not math.isfinite(float(value)):
        return missing("nonfinite diagnostic value")
    return {"status": "not_calibrated", "value": float(value), **extra}


def combine(states):
    if not states:
        return "not_evaluated"
    for state in ("fail", "not_calibrated", "not_evaluated"):
        if state in states:
            return state
    return "pass" if all(s == "pass" for s in states) else "not_evaluated"


def legacy_gate(evaluation):
    """Mirror the production !enforced || satisfied predicate, not its audit."""
    if not all(isinstance(evaluation.get(k), bool) for k in ("enforced", "satisfied")):
        return missing("production enforced/satisfied outcome not captured")
    return {"status": "pass" if not evaluation["enforced"] or evaluation["satisfied"] else "fail",
            "enforced": evaluation["enforced"], "satisfied": evaluation["satisfied"],
            "identity": "production predicate: !enforced || satisfied"}


def exact(x):
    return Fraction.from_float(float(x))


def decimal(x):
    with localcontext() as ctx:
        ctx.prec = 100
        return str(Decimal(x.numerator) / Decimal(x.denominator))


def load_state(path):
    import h5py
    with h5py.File(path) as f:
        meta = json.loads(f.attrs["metadata_json"])
        fields = {k: f["fields"][k][()] for k in f["fields"]}
    return meta, fields, compose_state(meta, fields)


def compose_state(meta, fields):
    """Reconstruct the saved representation, without regenerating populations."""
    if "packed_psi" not in fields:
        # SLOT population-memory checkpoints can be referenced but unpacked.
        # Do not invent a packed scale or a missing low coordinate.
        state = [exact(x) for x in fields["psi"]]
        for name, carrier in [("phin", "electron"), ("phip", "hole")]:
            ref, inc = carrier+"_qf_reference_V", carrier+"_qf_increment_V"
            if ref in fields and inc in fields:
                state.extend(exact(a)+exact(b) for a,b in zip(fields[ref],fields[inc]))
            else:
                state.extend(exact(x) for x in fields[name])
        return state
    scale = exact(meta["packed_potential_scale_V"])
    state = []
    for key, ref in [
        ("packed_psi", None),
        ("packed_electron_qf_increment", "electron_qf_reference_V"),
        ("packed_hole_qf_increment", "hole_qf_reference_V"),
    ]:
        hi = fields[key]
        lo = fields.get(key + "_low", [0.] * len(hi))
        state.extend((exact(h) + exact(l)) * scale +
                     (exact(fields[ref][i]) if ref else Fraction(0))
                     for i, (h, l) in enumerate(zip(hi, lo)))
    return state


def summarize(values, allowed, materials):
    indices = [i for i in range(len(values)) if allowed[i]]
    if not indices:
        return missing("no free active components")
    i = max(indices, key=lambda j: abs(values[j]))
    answer = metric(abs(values[i]), node=i, signed_value=float(values[i]),
                    count=len(indices), rms=math.sqrt(math.fsum(float(values[j])**2 for j in indices)/len(indices)))
    answer["per_material"] = {}
    for material in sorted({m for j in indices for m in materials[j]}):
        ids = [j for j in indices if material in materials[j]]
        j = max(ids, key=lambda k: abs(values[k]))
        answer["per_material"][material] = metric(abs(values[j]), node=j, count=len(ids))
    return answer


def eta_components(delta_V, vt):
    """Fixed T, DOS and band-edge offsets; includes psi in eta, never QF alone."""
    n = len(delta_V) // 3
    vt = exact(vt)
    return [
        [delta_V[i] / vt for i in range(n)],
        [(delta_V[i] - delta_V[n+i]) / vt for i in range(n)],
        [(delta_V[2*n+i] - delta_V[i]) / vt for i in range(n)],
    ]


def transition(expected, actual, allowed):
    ids = [i for i, ok in enumerate(allowed) if ok]
    changed = [i for i in ids if actual[i] != expected[i]]
    reduced = [i for i in ids if expected[i] and abs(actual[i]) < abs(expected[i])]
    zeroed = [i for i in ids if expected[i] and not actual[i]]
    worst = max(changed, key=lambda i: abs(actual[i]-expected[i]), default=None)
    relative = max((float(1-abs(actual[i]/expected[i])) for i in reduced), default=0.)
    return {"status": "not_calibrated", "changed_components": len(changed),
            "reduced_components": len(reduced), "zeroed_components": len(zeroed),
            "zeroed_indices": zeroed, "max_relative_reduction": relative,
            "max_abs_difference_V": decimal(abs(actual[worst]-expected[worst])) if worst is not None else "0",
            "max_index": worst,
            "zero_request_nonzero_actual": sum(not expected[i] and bool(actual[i]) for i in ids)}


def material_map(mesh, n):
    materials = [set() for _ in range(n)]
    triangles = mesh.get("triangles", mesh.get("cells", []))
    for region in mesh["regions"]:
        for ci in region["cell_ids"]:
            cell = triangles[ci]
            nodes = cell if isinstance(cell, list) else cell.get("node_ids", cell.get("nodes"))
            for node in nodes:
                materials[node].add(region["material"])
    return materials


def analyze_run(folder, contract, output_name="p1_replay"):
    get = lambda p: json.loads(Path(p).read_text())
    mesh = get(contract["template"]["mesh_file"])
    result = get(folder/"result.json") if (folder/"result.json").exists() else {}
    snapshots = sorted((folder / "progress").glob("iteration_*.json"),
                       key=lambda p: int(p.stem.split("_")[-1]))
    if Path(output_name).name != output_name or output_name in ("", ".", ".."):
        raise ValueError("output_name must be a single directory name")
    output = folder / output_name
    output.mkdir(exist_ok=True)
    reports, previous, previous_peak = [], None, None
    for path in snapshots:
        doc = get(path)
        state_file = path.with_suffix(".h5")
        meta, fields, state = load_state(state_file)
        n = len(state)//3
        materials = material_map(mesh, n)
        s = doc.get("p1_shadow")
        report = {"iteration": doc["iteration"], "input_sha256": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [path, state_file]},
            "arithmetic": "exact rational composition of binary64 high/low/reference/scale; 100-digit decimal output",
            "affects_acceptance": False}
        if s is None:
            report["E_adv"] = missing("no P1 direction capture")
            report["E_act"] = missing("no bound active/constraint mask and thermal constant")
        else:
            constrained = set(s["constrained_indices"])
            # Insulating-region carrier identities are included in the production
            # boundary map. Material filtering additionally prevents accidental use.
            semiconductor = [bool(m & {"Si", "Silicon", "PolySilicon"}) for m in materials]
            masks = [[i not in constrained for i in range(n)],
                     [semiconductor[i] and n+i not in constrained for i in range(n)],
                     [semiconductor[i] and 2*n+i not in constrained for i in range(n)]]
            flat_mask = sum(masks, [])
            report["coverage"] = {"free_counts": [sum(x) for x in masks],
                "all_materials": sorted(set.union(*materials)),
                "density_cutoff": None, "constraint_count": len(constrained)}
            report["linear_quality"] = s["linear_quality"]
            report["alpha"] = s["alpha"]
            report["E_phys"] = s["s0_physical_direction"]
            residual = s["residual_internal"]
            report["R_free_internal"] = {name: summarize(residual[b*n:(b+1)*n], masks[b], materials)
                                        for b, name in enumerate(["psi", "electron", "hole"])}
            report["R_constraint_internal"] = metric(max((abs(residual[i]) for i in constrained), default=0.),
                                                     count=len(constrained))
            report["physical_reference_scaling"] = s["fixed_reference_scales"]
            report["R_circuit"] = missing("fixed-device Newton solve; external circuit equation not present")
            report["legacy_row_gate_passed"] = missing("capture does not include original predicate result; use final result ledger")
            if previous is None or "s1_direction_internal" not in s:
                for k in ["E_adv", "E_act", "T_policy", "T_repr"]:
                    report[k] = missing("no previous accepted update")
                report["kappa"] = missing("no previous update norm")
            else:
                pmeta, pstate = previous
                if pmeta.get("split_mesh_fingerprint") != meta.get("split_mesh_fingerprint"):
                    raise ValueError("state mesh identity changed")
                if pmeta["packed_potential_scale_V"] != meta["packed_potential_scale_V"]:
                    raise ValueError("coordinate scale changed")
                scale, alpha = exact(s["potential_scale_V"]), exact(s["alpha"]["value"])
                raw = [exact(x)*scale for x in s["s1_direction_internal"]]
                processed = [exact(x)*scale*alpha for x in s["s3_direction_before_alpha_internal"]]
                actual = [y-x for x, y in zip(pstate, state)]
                adv_eta, act_eta = eta_components(raw, s["thermal_voltage_V"]), eta_components(actual, s["thermal_voltage_V"])
                report["E_adv"] = {name: summarize(adv_eta[b], masks[b], materials) for b, name in enumerate(["psi", "electron_eta", "hole_eta"])}
                report["E_act"] = {name: summarize(act_eta[b], masks[b], materials) for b, name in enumerate(["psi", "electron_eta", "hole_eta"])}
                for name in ("E_adv", "E_act"):
                    for block, entry in report[name].items():
                        if "node" not in entry:
                            continue
                        i = entry["node"]
                        entry["state_at_max"] = {
                            "psi_V": decimal(state[i]),
                            "electron_QF_V": decimal(state[n+i]),
                            "hole_QF_V": decimal(state[2*n+i]),
                            "electrons_cm3": float(fields["electrons_m3"][i])/1e6,
                            "holes_cm3": float(fields["holes_m3"][i])/1e6,
                            "materials": sorted(materials[i]),
                            "state_identity": "accepted candidate, not direction base state",
                        }
                report["eta_contract"] = "fixed isothermal DOS/BGN offsets in this same-config QF solve; delta_eta_n=(delta_psi-delta_phi_n)/Vt, holes opposite"
                report["delta_log_DOS"] = metric(0., basis="frozen isothermal material/doping configuration; no parameter update")
                report["nonlinear_map"] = {"identity": True, "coordinate": "QF", "not_truncation": True}
                report["T_policy"] = transition(raw, processed, flat_mask)
                report["T_repr"] = transition(processed, actual, flat_mask)
                peak = max(float(abs(v)) for b in range(3) for i, v in enumerate(adv_eta[b]) if masks[b][i])
                if previous_peak is None:
                    report["kappa"] = missing("no previous direction")
                elif previous_peak == 0:
                    report["kappa"] = {"status": "undefined", "reason": "zero previous norm", "both_zero": peak == 0}
                else:
                    report["kappa"] = metric(peak/previous_peak, interpretable=False,
                        reason="branch continuity not captured; no contraction inference")
                previous_peak = peak
                with (output / f"update_{doc['iteration']}.csv").open("w", newline="") as f:
                    w = csv.writer(f);w.writerow(["index","node","block","participates","s1_V","s3_V","s4_V","policy_delta_V","representation_delta_V"])
                    for i in range(3*n):
                        w.writerow([i,i%n,i//n,flat_mask[i]]+[decimal(v) for v in
                            [raw[i],processed[i],actual[i],processed[i]-raw[i],actual[i]-processed[i]]])
        # Missing production guard outcomes must not be inferred from a small
        # residual or from the fact that a state file was written.
        guards = {name: missing("original predicate outcome not captured") for name in
                  ["inverse_mapping", "boundary_identity", "global_closure",
                   "contact_QF_branch", "numerical_protections", "file_consistency"]}
        finite = all(math.isfinite(float(x)) for x in state)
        guards["finite_archived_state"] = {"status": "pass" if finite else "fail",
                                          "scope": "saved state only"}
        report["G_retained"] = guards
        if s:
            guard_aliases = {"global_continuity_closure": "global_closure",
                             "contact_majority_qf_branch": "contact_QF_branch",
                             "archive_consistency": "file_consistency"}
            for key, value in s.get("retained_guards", {}).items():
                guards[guard_aliases.get(key, key)] = value
            for key in ("solved_equation_quality", "legacy_row_gate_passed", "global_ledger"):
                if key in s:
                    report[key] = s[key]
            if "continuity_residual_to_A_per_um" in s:
                factor = s["continuity_residual_to_A_per_um"]
                report["continuity_residual_to_A_per_um"] = factor
                report["R_free_A_per_um"] = {}
                for block, name in ((1, "electron"), (2, "hole")):
                    values = [v * factor for v in s["residual_internal"][block*n:(block+1)*n]]
                    report["R_free_A_per_um"][name] = summarize(values, masks[block], materials)
                report["physical_normalization"] = missing("unit conversion verified; fixed references not yet selected")
            if "poisson_residual_to_C_per_m" in s:
                report["R_free_C_per_m"] = summarize(
                    [v*s["poisson_residual_to_C_per_m"] for v in s["residual_internal"][:n]],
                    masks[0], materials)
                constraint_values = [s["residual_internal"][i]*s["constraint_residual_to_V"] for i in constrained]
                report["R_constraint_V"] = metric(max(map(abs, constraint_values), default=0.),
                    count=len(constraint_values), rms=math.sqrt(math.fsum(v*v for v in constraint_values)/len(constraint_values))
                    if constraint_values else 0.) if constraint_values else missing("no constraint rows")
                report["physical_normalization"] = s["fixed_reference_scales"]
        report["shadow_numerically_converged"] = {"status": combine(
            ["not_calibrated"] + [x["status"] for x in guards.values()])}
        report["shadow_device_audit_passed"] = missing("reference tolerances and complete retained-guard outcomes not captured")
        report["shadow_continuation_eligible"] = {"status": combine(
            [report["shadow_numerically_converged"]["status"], "not_evaluated"])}
        reports.append(report)
        previous = (meta, state)
    covered_solve = any("solved_equation_quality" in x and x["iteration"] > 0 for x in reports)
    covered_units = any("R_free_C_per_m" in x for x in reports)
    not_covered = ["independent physical direction", "complete retained-guard labels",
                   "per-region source/flux audit", "independent circuit and path audit"]
    if not covered_solve:
        not_covered.append("closure-replaced returned-coordinate solve quality")
    if not covered_units:
        not_covered.append("physical fixed residual references")
    final = {"schema":"vela.p1_shadow.replay.v1","snapshots":reports,
        "legacy_result": result, "observer_does_not_override_legacy":True,
        "legacy_row_gate_passed": legacy_gate(result.get("carrier_row_convergence", {})),
        "legacy_closure_gate_passed": legacy_gate(result.get("global_continuity_closure", {})),
        "final_device_raw_records": result.get("convergence_report", {}).get("device_checks", {}),
        "final_device_audit_status": missing("record availability does not establish independent accuracy tolerances"),
        "not_covered":not_covered,
        "legacy_full_qualification":missing("Newton result is not full outer circuit qualification")}
    (output/"summary.json").write_text(json.dumps(final,indent=2,allow_nan=False))
    return final


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("root",type=Path)
    parser.add_argument("labels",nargs="+")
    parser.add_argument("--output-name",default="p1_replay")
    args=parser.parse_args()
    c=json.loads((args.root/"contract.json").read_text())
    for label in args.labels:
        answer=analyze_run(args.root/label,c,args.output_name)
        print(json.dumps({"label":label,"snapshots":len(answer["snapshots"]),
                          "last":answer["snapshots"][-1] if answer["snapshots"] else None}))


if __name__ == "__main__":
    main()
