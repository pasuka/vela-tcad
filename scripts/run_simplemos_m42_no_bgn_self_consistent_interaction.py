#!/usr/bin/env python3
"""Run the SimpleMOS M42 no-BGN self-consistent interaction audit."""

from __future__ import annotations

import argparse
import copy
import math
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import run_simplemos_m10_fixed_state_replay as m10
import run_simplemos_m31_minority_poisson_perturbation as m31
import run_simplemos_m40_true_no_bgn_factorial as m40
import run_simplemos_m41_equal_ni_flux_ablation as m41


REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "build-release/vela_example_runner.exe"
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m42_no_bgn_self_consistent_interaction_contract_v1.json")
M41_OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
              / "m41_equal_ni_flux_ablation")
M41_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "simplemos_m41_equal_ni_flux_ablation_evidence.json")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m42_no_bgn_self_consistent_interaction")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "no_bgn_self_consistent_interaction")
MODES = m41.MODES
CONTACTS = {
    "dominant": "dominant_signed_contact_mean",
    "node_local": "legacy_node_local",
}
SRH_VALUES = {"srh_on": True, "srh_off": False}
Q = 1.602176634e-19


def set_operator(deck: dict[str, Any], mode: str, srh: bool,
                 contact: str) -> None:
    m41.configure(deck, srh, MODES[mode])
    deck["solver"]["contact_boundary_reconstruction"] = CONTACTS[contact]


def target_row(curve: Path) -> dict[str, str]:
    return min(m40.read_csv(curve),
               key=lambda row: abs(float(row["bias_V"]) - 0.05))


def self_consistent_workflow(mode: str, srh_name: str, contact: str,
                             force: bool) -> dict[str, Any]:
    srh = SRH_VALUES[srh_name]
    cell = "bgn_off_srh_on" if srh else "bgn_off_srh_off"
    if contact == "dominant":
        upstream = M41_OUTPUT / "self_consistent" / cell / mode / "20_gate"
        if not (upstream / "state.csv").is_file():
            raise FileNotFoundError(f"missing M41 state: {upstream}")
        row = target_row(upstream / "curve.csv")
        return {
            "mode": mode,
            "srh_name": srh_name,
            "srh": srh,
            "contact": contact,
            "contact_boundary_reconstruction": CONTACTS[contact],
            "reported_current_A_per_um": abs(float(
                row["current_total_A_per_um"])),
            "state": m40.portable(upstream / "state.csv"),
            "config": m40.portable(upstream / "config.json"),
            "converged": True,
            "reused_m41": True,
            "reused": True,
        }

    root = OUTPUT / "self_consistent" / mode / srh_name / contact
    if force and root.exists():
        shutil.rmtree(root)
    final = root / "20_gate"
    if (not force and (final / "state.csv").is_file()
            and (final / "curve.csv").is_file()):
        row = target_row(final / "curve.csv")
        return {
            "mode": mode, "srh_name": srh_name, "srh": srh,
            "contact": contact,
            "contact_boundary_reconstruction": CONTACTS[contact],
            "reported_current_A_per_um": abs(float(
                row["current_total_A_per_um"])),
            "state": m40.portable(final / "state.csv"),
            "config": m40.portable(final / "config.json"),
            "converged": True, "reused_m41": False, "reused": True,
        }

    previous: Path | None = None
    phases = (("00_equilibrium", "00_equilibrium.json"),
              ("10_drain", "10_drain_ramp.json"),
              ("20_gate", "20_gate_sweep.json"))
    converged = True
    for label, source in phases:
        phase = root / label
        phase.mkdir(parents=True, exist_ok=True)
        deck = copy.deepcopy(m40.read_json(m40.VELA_BASE / source))
        set_operator(deck, mode, srh, contact)
        deck["output_csv"] = str((phase / "curve.csv").resolve())
        deck["log_file"] = str((phase / "run.log").resolve())
        deck["sweep"]["write_state_file"] = str(
            (phase / "state.csv").resolve())
        if previous is not None:
            deck["sweep"]["initial_state_file"] = str(previous.resolve())
        if label == "20_gate":
            deck["sweep"].update({
                "start": 0.0, "stop": 0.05, "step": 0.05,
                "bias_points": [0.0, 0.05], "max_step": 0.05,
            })
        deck["simplemos_m42"] = {
            "mode": mode, "srh": srh,
            "contact_boundary_reconstruction": CONTACTS[contact],
            "default_model_changed": False,
        }
        config = phase / "config.json"
        m40.write_json(config, deck)
        status = m10.execute_runner(config, RUNNER)
        converged = converged and bool(status.get("converged"))
        if not status.get("converged"):
            raise RuntimeError(
                f"M42 workflow failed: {mode} {srh_name} {contact} {label}")
        previous = phase / "state.csv"
    row = target_row(final / "curve.csv")
    return {
        "mode": mode, "srh_name": srh_name, "srh": srh,
        "contact": contact,
        "contact_boundary_reconstruction": CONTACTS[contact],
        "reported_current_A_per_um": abs(float(
            row["current_total_A_per_um"])),
        "state": m40.portable(final / "state.csv"),
        "config": m40.portable(final / "config.json"),
        "converged": converged, "reused_m41": False, "reused": False,
    }


def contact_sets(config: Path) -> tuple[set[int], set[int]]:
    deck = m40.read_json(config)
    mesh = m40.read_json(Path(deck["mesh_file"]))
    all_contacts = {int(node) for item in mesh["contacts"]
                    for node in item["node_ids"]}
    drain = next({int(node) for node in item["node_ids"]}
                 for item in mesh["contacts"]
                 if item["name"].lower() == "drain")
    return all_contacts, drain


def drain_cut_condition(rows: list[dict[str, str]], nodes: set[int]) -> float:
    contributions: list[float] = []
    for row in rows:
        n0, n1 = int(row["node0"]), int(row["node1"])
        at0, at1 = n0 in nodes, n1 in nodes
        if at0 == at1:
            continue
        sign = 1.0 if at0 else -1.0
        contributions.extend([
            -Q * sign * float(
                row["electron_particle_line_flux_per_m_s"]) * 1.0e-6,
            Q * sign * float(
                row["hole_particle_line_flux_per_m_s"]) * 1.0e-6,
        ])
    return sum(abs(value) for value in contributions) / max(
        abs(sum(contributions)), 1.0e-300)


def probe(state_row: dict[str, Any], op_mode: str, op_srh_name: str,
          op_contact: str, root: Path, force: bool) -> dict[str, Any]:
    root.mkdir(parents=True, exist_ok=True)
    state = REPO / state_row["state"]
    source_config = REPO / state_row["config"]
    outputs: dict[str, Path] = {}
    for label, simulation_type in (
            ("sg", "sg_edge_flux_probe"),
            ("terms", "newton_carrier_term_probe")):
        output = root / f"{label}.csv"
        config = root / f"{label}.json"
        deck = m10.probe_config(source_config, state, output, 0.05, 0.05,
                                simulation_type)
        set_operator(deck, op_mode, SRH_VALUES[op_srh_name], op_contact)
        if label == "terms":
            deck["carrier_term_probe"] = {"solved_equation_terms": False}
        deck["simplemos_m42_probe"] = {
            "state_mode": state_row["mode"],
            "state_srh": state_row["srh_name"],
            "state_contact": state_row["contact"],
            "operator_mode": op_mode,
            "operator_srh": op_srh_name,
            "operator_contact": op_contact,
        }
        m40.write_json(config, deck)
        if force or not output.is_file():
            status = m10.execute_runner(config, RUNNER)
            if not status.get("converged"):
                raise RuntimeError(f"M42 probe failed: {config}")
        outputs[label] = output

    contacts, drain = contact_sets(source_config)
    sg = m40.read_csv(outputs["sg"])
    cut = m10.drain_cut_current(sg, drain)
    terms = m40.read_csv(outputs["terms"])
    free = [row for row in terms if int(row["node_id"]) not in contacts]
    boundary = [row for row in terms if int(row["node_id"]) in contacts]

    def l2(rows: list[dict[str, str]], column: str) -> float:
        return math.sqrt(sum(float(row[column]) ** 2 for row in rows))

    return {
        "state_mode": state_row["mode"],
        "state_srh": state_row["srh_name"],
        "state_contact": state_row["contact"],
        "operator_mode": op_mode,
        "operator_srh": op_srh_name,
        "operator_contact": op_contact,
        "drain_cut_total_A_per_um": cut["total_A_per_um"],
        "drain_cut_electron_A_per_um": cut["electron_A_per_um"],
        "drain_cut_condition": drain_cut_condition(sg, drain),
        "free_electron_residual_l2": l2(free, "electron_residual"),
        "free_hole_residual_l2": l2(free, "hole_residual"),
        "free_electron_flux_l2": l2(free, "electron_flux"),
        "free_electron_recombination_l2": l2(
            free, "electron_recombination"),
        "contact_electron_boundary_l2": l2(
            boundary, "electron_boundary"),
        "contact_hole_boundary_l2": l2(boundary, "hole_boundary"),
        "sg": m40.portable(outputs["sg"]),
        "terms": m40.portable(outputs["terms"]),
    }


def state_rows(path: Path) -> tuple[list[str], dict[int, dict[str, str]]]:
    return m31.state_map(path)


def hybrid_probes(states: dict[tuple[str, str, str], dict[str, Any]],
                  force: bool) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    specifications = {
        "srh_on_baseline": (),
        "absolute_qf_only_roundtrip": None,
        "psi_from_srh_off": ("psi",),
        "electron_qf_density_from_srh_off": (
            "phin", "electrons_m3", "electron_qf_increment_V",
            "electron_qf_reference_V"),
        "hole_qf_density_from_srh_off": (
            "phip", "holes_m3", "hole_qf_increment_V",
            "hole_qf_reference_V"),
        "all_carriers_from_srh_off": (
            "phin", "electrons_m3", "electron_qf_increment_V",
            "electron_qf_reference_V", "phip", "holes_m3",
            "hole_qf_increment_V", "hole_qf_reference_V"),
        "full_srh_off_state": (
            "psi", "phin", "electrons_m3", "electron_qf_increment_V",
            "electron_qf_reference_V", "phip", "holes_m3",
            "hole_qf_increment_V", "hole_qf_reference_V"),
    }
    for mode in MODES:
        on = states[(mode, "srh_on", "dominant")]
        off = states[(mode, "srh_off", "dominant")]
        fields, baseline = state_rows(REPO / on["state"])
        _, replacement = state_rows(REPO / off["state"])
        for variant, replace_fields in specifications.items():
            root = OUTPUT / "qf_substitution" / mode / variant
            state = root / "state.csv"
            if replace_fields is None:
                core_fields = ["node_id", "psi", "phin", "phip",
                               "electrons_m3", "holes_m3"]
                m31.write_state(state, core_fields, baseline)
            elif replace_fields:
                m31.hybrid_state(fields, baseline, replacement,
                                 replace_fields, None, state)
            else:
                state.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(REPO / on["state"], state)
            hybrid = dict(on)
            hybrid["state"] = m40.portable(state)
            result = probe(hybrid, mode, "srh_on", "dominant",
                           root / "probe", True)
            results.append({
                "mode": mode,
                "variant": variant,
                "drain_cut_total_A_per_um": result[
                    "drain_cut_total_A_per_um"],
                "drain_cut_condition": result["drain_cut_condition"],
                "state": m40.portable(state),
                "sg": result["sg"],
            })
    return results


def maximum_state_delta(left: Path, right: Path) -> dict[str, float]:
    _, a = state_rows(left)
    _, b = state_rows(right)
    result: dict[str, float] = {}
    for field in ("psi", "phin", "phip"):
        result[f"maximum_abs_{field}_delta_V"] = max(
            abs(float(a[node][field]) - float(b[node][field])) for node in a)
    for field, label in (("electrons_m3", "logn"), ("holes_m3", "logp")):
        result[f"maximum_abs_{label}_delta_dex"] = max(
            abs(math.log10(max(float(a[node][field]), 1e-300))
                - math.log10(max(float(b[node][field]), 1e-300)))
            for node in a)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    contract = m40.read_json(CONTRACT)
    if contract.get("schema") != (
            "vela.simplemos.sdevice.m42_no_bgn_self_consistent_interaction.v1"):
        raise ValueError("unexpected M42 contract schema")
    if not M41_EVIDENCE.is_file():
        raise FileNotFoundError("M42 requires frozen M41 evidence")

    workflow_jobs = [(mode, srh, contact, args.force)
                     for mode in MODES for srh in SRH_VALUES
                     for contact in CONTACTS]
    with ThreadPoolExecutor(max_workers=min(args.jobs, len(workflow_jobs))) as pool:
        workflows = list(pool.map(
            lambda values: self_consistent_workflow(*values), workflow_jobs))
    states = {(row["mode"], row["srh_name"], row["contact"]): row
              for row in workflows}

    diagonal_jobs = []
    for row in workflows:
        root = (OUTPUT / "diagonal" / row["mode"] / row["srh_name"]
                / row["contact"])
        diagonal_jobs.append((row, row["mode"], row["srh_name"],
                              row["contact"], root, args.force))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        diagonal = list(pool.map(lambda values: probe(*values), diagonal_jobs))
    diagonal_map = {(row["state_mode"], row["state_srh"],
                     row["state_contact"]): row for row in diagonal}

    self_consistent_rows = []
    for workflow in workflows:
        key = (workflow["mode"], workflow["srh_name"], workflow["contact"])
        direct = diagonal_map[key]
        reported = workflow["reported_current_A_per_um"]
        cut = abs(direct["drain_cut_total_A_per_um"])
        self_consistent_rows.append({
            **workflow,
            "operator_consistent_cut_A_per_um": cut,
            "cut_vs_reported_dex": abs(math.log10(
                max(cut, 1e-300) / max(reported, 1e-300))),
            "free_electron_residual_l2": direct[
                "free_electron_residual_l2"],
            "free_hole_residual_l2": direct["free_hole_residual_l2"],
            "drain_cut_condition": direct["drain_cut_condition"],
        })
    m40.write_csv(PORTABLE / "m42_self_consistent_matrix.csv",
                  self_consistent_rows)

    default_states = [states[(mode, srh, "dominant")]
                      for mode in MODES for srh in SRH_VALUES]
    cross_jobs = []
    for state in default_states:
        for op_mode in MODES:
            for op_srh in SRH_VALUES:
                for op_contact in CONTACTS:
                    root = (OUTPUT / "frozen_operator_matrix"
                            / f"state_{state['mode']}_{state['srh_name']}"
                            / f"op_{op_mode}_{op_srh}_{op_contact}")
                    cross_jobs.append((state, op_mode, op_srh, op_contact,
                                       root, args.force))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        frozen_matrix = list(pool.map(lambda values: probe(*values), cross_jobs))
    m40.write_csv(PORTABLE / "m42_frozen_operator_matrix.csv", frozen_matrix)

    substitutions = hybrid_probes(states, args.force)
    for mode in MODES:
        rows = [row for row in substitutions if row["mode"] == mode]
        by_variant = {row["variant"]: row for row in rows}
        baseline = abs(by_variant["srh_on_baseline"][
            "drain_cut_total_A_per_um"])
        full = abs(by_variant["full_srh_off_state"][
            "drain_cut_total_A_per_um"])
        full_shift = math.log10(max(full, 1e-300) / max(baseline, 1e-300))
        for row in rows:
            current = abs(row["drain_cut_total_A_per_um"])
            shift = math.log10(max(current, 1e-300) / max(baseline, 1e-300))
            row["log_current_shift_from_srh_on_dex"] = shift
            row["fraction_of_full_srh_state_shift"] = (
                shift / full_shift if abs(full_shift) > 1e-15 else 0.0)
    m40.write_csv(PORTABLE / "m42_qf_state_substitution.csv", substitutions)

    contact_rows = []
    by_self = {(row["mode"], row["srh_name"], row["contact"]): row
               for row in self_consistent_rows}
    for mode in MODES:
        for srh in SRH_VALUES:
            dominant = by_self[(mode, srh, "dominant")]
            local = by_self[(mode, srh, "node_local")]
            contact_rows.append({
                "mode": mode, "srh_name": srh,
                "operator_cut_log_shift_dex": math.log10(
                    max(local["operator_consistent_cut_A_per_um"], 1e-300)
                    / max(dominant["operator_consistent_cut_A_per_um"], 1e-300)),
                "reported_current_log_shift_dex": math.log10(
                    max(local["reported_current_A_per_um"], 1e-300)
                    / max(dominant["reported_current_A_per_um"], 1e-300)),
                **maximum_state_delta(REPO / dominant["state"],
                                      REPO / local["state"]),
            })
    m40.write_csv(PORTABLE / "m42_contact_reconstruction_response.csv",
                  contact_rows)

    primary = [row for row in frozen_matrix
               if row["state_contact"] == "dominant"]
    srh_direct_current_delta = 0.0
    contact_direct_current_delta = 0.0
    transport_direct_shifts: list[float] = []
    srh_recombination_changes: list[float] = []
    for state in default_states:
        matching = [row for row in primary
                    if row["state_mode"] == state["mode"]
                    and row["state_srh"] == state["srh_name"]
                    and row["operator_mode"] == state["mode"]]
        on = next(row for row in matching
                  if row["operator_srh"] == "srh_on"
                  and row["operator_contact"] == "dominant")
        off = next(row for row in matching
                   if row["operator_srh"] == "srh_off"
                   and row["operator_contact"] == "dominant")
        local = next(row for row in matching
                     if row["operator_srh"] == state["srh_name"]
                     and row["operator_contact"] == "node_local")
        diagonal_row = next(row for row in matching
                            if row["operator_srh"] == state["srh_name"]
                            and row["operator_contact"] == "dominant")
        transport = next(row for row in primary
                         if row["state_mode"] == state["mode"]
                         and row["state_srh"] == state["srh_name"]
                         and row["operator_mode"] != state["mode"]
                         and row["operator_srh"] == state["srh_name"]
                         and row["operator_contact"] == "dominant")
        transport_direct_shifts.append(abs(math.log10(
            max(abs(transport["drain_cut_total_A_per_um"]), 1e-300)
            / max(abs(diagonal_row["drain_cut_total_A_per_um"]), 1e-300))))
        srh_recombination_changes.append(abs(
            on["free_electron_recombination_l2"]
            - off["free_electron_recombination_l2"]))
        srh_direct_current_delta = max(
            srh_direct_current_delta,
            abs(on["drain_cut_total_A_per_um"]
                - off["drain_cut_total_A_per_um"]))
        contact_direct_current_delta = max(
            contact_direct_current_delta,
            abs(local["drain_cut_total_A_per_um"]
                - diagonal_row["drain_cut_total_A_per_um"]))

    legacy_gaps = [row["cut_vs_reported_dex"]
                   for row in self_consistent_rows if row["mode"] == "legacy"]
    compensated_gaps = [row["cut_vs_reported_dex"]
                        for row in self_consistent_rows
                        if row["mode"] == "compensated"]
    qf_by = {(row["mode"], row["variant"]): row for row in substitutions}
    electron_fractions = {
        mode: qf_by[(mode, "electron_qf_density_from_srh_off")][
            "fraction_of_full_srh_state_shift"] for mode in MODES}
    compensated_absolute_roundtrip_loss = abs(math.log10(
        max(abs(qf_by[("compensated", "absolute_qf_only_roundtrip")][
            "drain_cut_total_A_per_um"]), 1e-300)
        / max(abs(qf_by[("compensated", "srh_on_baseline")][
            "drain_cut_total_A_per_um"]), 1e-300)))
    legacy_residuals = [row["free_electron_residual_l2"]
                        for row in self_consistent_rows
                        if row["mode"] == "legacy"]
    compensated_residuals = [row["free_electron_residual_l2"]
                             for row in self_consistent_rows
                             if row["mode"] == "compensated"]
    checks = {
        "self_consistent_state_count": len(workflows) == 8,
        "all_self_consistent_workflows_converged": all(
            row["converged"] for row in workflows),
        "frozen_operator_pair_count": len(frozen_matrix) == 32,
        "qf_substitution_count": len(substitutions) == 14,
        "compensated_cut_reported_closure": max(compensated_gaps) < 0.01,
        "legacy_cut_reported_mismatch_reproduced": min(legacy_gaps) > 1.0,
        "srh_frozen_cut_is_null": srh_direct_current_delta < 1e-30,
        "contact_frozen_cut_is_null": contact_direct_current_delta < 1e-30,
        "qf_reference_increment_required": (
            compensated_absolute_roundtrip_loss > 2.0),
    }
    report = {
        "schema": "vela.simplemos.sdevice.m42_no_bgn_self_consistent_interaction_report.v1",
        "status": "complete" if all(checks.values()) else "failed",
        "execution": {
            "default_model_changed": False,
            "new_sentaurus_execution": False,
            "self_consistent_state_count": len(workflows),
            "frozen_operator_pair_count": len(frozen_matrix),
            "qf_substitution_count": len(substitutions),
        },
        "findings": {
            "maximum_compensated_cut_reported_gap_dex": max(compensated_gaps),
            "minimum_legacy_cut_reported_gap_dex": min(legacy_gaps),
            "maximum_srh_direct_frozen_cut_delta_A_per_um": (
                srh_direct_current_delta),
            "maximum_contact_direct_frozen_cut_delta_A_per_um": (
                contact_direct_current_delta),
            "maximum_contact_self_consistent_cut_shift_abs_dex": max(
                abs(row["operator_cut_log_shift_dex"])
                for row in contact_rows),
            "maximum_contact_state_phin_delta_V": max(
                row["maximum_abs_phin_delta_V"] for row in contact_rows),
            "electron_qf_density_fraction_of_srh_shift": electron_fractions,
            "maximum_transport_direct_frozen_cut_shift_abs_dex": max(
                transport_direct_shifts),
            "maximum_srh_recombination_operator_change_l2": max(
                srh_recombination_changes),
            "minimum_legacy_free_electron_residual_l2": min(
                legacy_residuals),
            "maximum_compensated_free_electron_residual_l2": max(
                compensated_residuals),
            "legacy_to_compensated_residual_separation_min_ratio": (
                min(legacy_residuals) / max(compensated_residuals)),
            "compensated_absolute_qf_roundtrip_loss_dex": (
                compensated_absolute_roundtrip_loss),
        },
        "causal_interpretation": {
            "srh_terminal_effect_is_state_feedback_not_direct_operator": (
                srh_direct_current_delta < 1e-30),
            "contact_reconstruction_is_null_for_simplemos": max(
                abs(row["operator_cut_log_shift_dex"])
                for row in contact_rows) < 1e-8,
            "electron_qf_density_is_primary_srh_feedback_path": all(
                abs(value) > 0.8 for value in electron_fractions.values()),
            "legacy_state_terminal_extraction_is_operator_inconsistent": (
                min(legacy_gaps) > 1.0),
            "qf_reference_increment_representation_is_required": (
                compensated_absolute_roundtrip_loss > 2.0),
        },
        "acceptance": {**checks, "all_checks_pass": all(checks.values())},
        "artifacts": {
            "self_consistent_matrix": m40.portable(
                PORTABLE / "m42_self_consistent_matrix.csv"),
            "frozen_operator_matrix": m40.portable(
                PORTABLE / "m42_frozen_operator_matrix.csv"),
            "qf_state_substitution": m40.portable(
                PORTABLE / "m42_qf_state_substitution.csv"),
            "contact_response": m40.portable(
                PORTABLE / "m42_contact_reconstruction_response.csv"),
        },
        "claim_policy": contract["claim_policy"],
    }
    report_path = PORTABLE / "m42_no_bgn_self_consistent_interaction_report.json"
    m40.write_json(report_path, report)
    artifact_paths = [report_path,
                      PORTABLE / "m42_self_consistent_matrix.csv",
                      PORTABLE / "m42_frozen_operator_matrix.csv",
                      PORTABLE / "m42_qf_state_substitution.csv",
                      PORTABLE / "m42_contact_reconstruction_response.csv"]
    source_paths = [CONTRACT, Path(__file__).resolve(), M41_EVIDENCE,
                    REPO / "src/equation/CoupledDDAssembler.cpp",
                    REPO / "src/solver/NewtonSolver.cpp"]
    evidence = {
        "schema": "vela.simplemos.sdevice.m42_no_bgn_self_consistent_interaction_evidence.v1",
        "status": "frozen" if report["status"] == "complete" else "failed",
        "artifacts": [{"path": m40.portable(path),
                       "sha256": m10.sha256(path)} for path in artifact_paths],
        "source_hashes": {m40.portable(path): m10.sha256(path)
                          for path in source_paths},
        "default_model_changed": False,
        "acceptance": report["acceptance"],
    }
    m40.write_json(REPO / "reference_tcad/simplemos_sentaurus2022"
                   / "simplemos_m42_no_bgn_self_consistent_interaction_evidence.json",
                   evidence)
    print(report["findings"])
    print(report["causal_interpretation"])
    return 0 if report["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
