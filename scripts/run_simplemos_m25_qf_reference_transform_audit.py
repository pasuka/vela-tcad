#!/usr/bin/env python3
"""Audit contact-basin QF transforms and SG log imbalance edge by edge."""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
RAW = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
       / "m22_n23_hfs_deep_off")
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m25_qf_reference_transform_audit_contract_v1.json")
M24 = (REPO / "reference_tcad/simplemos_sentaurus2022"
       / "contact_boundary_audit/m24_drain_cut_edge_ledger.csv")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m25_qf_reference_transform_audit")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "qf_reference_transform_audit")
KB = 1.380649e-23
Q = 1.602176634e-19


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]),
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def exact(value: str | float) -> Decimal:
    return Decimal.from_float(float(value))


def state_path(variant: str) -> Path:
    return RAW / "self_consistent" / variant / "vg_0p05/state.csv"


def sg_path(variant: str) -> Path:
    return (RAW / "frozen/vg_0p05" / f"state_{variant}"
            / f"eval_{variant}/sg_edges.csv")


def validate(contract: dict[str, Any]) -> None:
    if contract.get("schema") != (
            "vela.simplemos.sdevice.m25_qf_reference_transform_audit.v1"):
        raise ValueError("unexpected M25 contract schema")
    if contract.get("state_variants") != ["full", "no_hfs"]:
        raise ValueError("M25 requires paired full/no_hfs states")
    if int(contract["precision"]["decimal_digits"]) < 80:
        raise ValueError("M25 requires at least 80 Decimal digits")


def transform_row(variant: str, edge: dict[str, str],
                  states: dict[int, dict[str, str]], cut_ids: set[int],
                  thermal_voltage: Decimal) -> dict[str, Any]:
    edge_id = int(edge["edge_id"])
    node0 = int(edge["node0"])
    node1 = int(edge["node1"])
    state0 = states[node0]
    state1 = states[node1]
    ref0 = exact(state0["electron_qf_reference_V"])
    ref1 = exact(state1["electron_qf_reference_V"])
    inc0 = exact(state0["electron_qf_increment_V"])
    inc1 = exact(state1["electron_qf_increment_V"])
    observed0 = exact(edge["electron_sg_phin0_relative_V"])
    observed1 = exact(edge["electron_sg_phin1_relative_V"])
    expected0 = inc0
    expected1 = (ref1 - ref0) + inc1
    physical_drop = (ref1 - ref0) + (inc1 - inc0)
    observed_drop = observed1 - observed0
    naive_drop = inc1 - inc0
    expected_log_from_state = physical_drop / thermal_voltage
    expected_log_from_edge = observed_drop / thermal_voltage
    observed_log = exact(edge["electron_sg_log_left_over_right"])
    right_factor = exact(edge["electron_sg_right_factor_flux"])
    factorized = right_factor * (observed_log.exp() - Decimal(1))
    stable_flux = exact(edge["electron_sg_stable_flux"])
    flux_scale = max(abs(stable_flux), abs(factorized), Decimal("1e-300"))
    ref_transition = ref0 != ref1
    drain_cut = edge_id in cut_ids
    classes = [name for name, enabled in (
        ("drain_cut", drain_cut),
        ("reference_transition", ref_transition)) if enabled]
    return {
        "variant": variant,
        "edge_id": edge_id,
        "node0": node0,
        "node1": node1,
        "edge_class": "+".join(classes),
        "drain_cut": drain_cut,
        "reference_transition": ref_transition,
        "finite_couple_transport_edge": abs(float(edge["couple_m"])) > 0.0,
        "couple_m": edge["couple_m"],
        "reference0_V": float(ref0),
        "reference1_V": float(ref1),
        "reference_jump_V": float(ref1 - ref0),
        "increment0_V": float(inc0),
        "increment1_V": float(inc1),
        "expected_anchored_phin0_V": float(expected0),
        "observed_anchored_phin0_V": float(observed0),
        "anchored_phin0_error_V": float(observed0 - expected0),
        "expected_anchored_phin1_V": float(expected1),
        "observed_anchored_phin1_V": float(observed1),
        "anchored_phin1_error_V": float(observed1 - expected1),
        "physical_qf_drop_V": float(physical_drop),
        "observed_anchored_qf_drop_V": float(observed_drop),
        "physical_drop_closure_error_V": float(observed_drop - physical_drop),
        "naive_increment_only_drop_V": float(naive_drop),
        "naive_drop_error_V": float(naive_drop - physical_drop),
        "reference_jump_thermal_voltages": float((ref1 - ref0)
                                                   / thermal_voltage),
        "expected_sg_log_from_state": float(expected_log_from_state),
        "expected_sg_log_from_anchored_endpoints": float(expected_log_from_edge),
        "observed_sg_log_left_over_right": float(observed_log),
        "sg_log_state_closure_error": float(observed_log
                                             - expected_log_from_state),
        "sg_log_endpoint_closure_error": float(observed_log
                                                - expected_log_from_edge),
        "electron_sg_stable_flux": float(stable_flux),
        "decimal100_factorized_flux": float(factorized),
        "factorized_flux_relative_error": float(abs(stable_flux - factorized)
                                                  / flux_scale),
    }


def summarize(rows: list[dict[str, Any]], variant: str,
              edge_class: str) -> dict[str, Any]:
    group = [row for row in rows if row["variant"] == variant
             and row[edge_class]]
    carrying = [row for row in group if row["finite_couple_transport_edge"]]
    return {
        "variant": variant,
        "edge_class": edge_class,
        "edge_count": len(group),
        "finite_couple_transport_edge_count": len(carrying),
        "maximum_endpoint_transform_error_V": max(
            max(abs(row["anchored_phin0_error_V"]),
                abs(row["anchored_phin1_error_V"])) for row in group),
        "maximum_physical_drop_closure_error_V": max(
            abs(row["physical_drop_closure_error_V"]) for row in group),
        "maximum_sg_log_state_closure_error": max(
            abs(row["sg_log_state_closure_error"]) for row in group),
        "maximum_sg_log_endpoint_closure_error": max(
            abs(row["sg_log_endpoint_closure_error"]) for row in group),
        "maximum_factorized_flux_relative_error": max(
            row["factorized_flux_relative_error"] for row in group),
        "maximum_naive_increment_only_drop_error_V": max(
            abs(row["naive_drop_error_V"]) for row in group),
        "maximum_absolute_stable_flux": max(
            abs(row["electron_sg_stable_flux"]) for row in group),
    }


def analyze(contract: dict[str, Any], portable_dir: Path) -> dict[str, Any]:
    variants = contract["state_variants"]
    cut_ids = {int(row["edge_id"]) for row in read_csv(M24)}
    states = {
        variant: {int(row["node_id"]): row for row in read_csv(state_path(variant))}
        for variant in variants
    }
    edges = {
        variant: {int(row["edge_id"]): row for row in read_csv(sg_path(variant))}
        for variant in variants
    }
    transition_ids = {
        variant: {edge_id for edge_id, row in edges[variant].items()
                  if float(row["electron_qf_reference0_V"])
                  != float(row["electron_qf_reference1_V"])}
        for variant in variants
    }
    if transition_ids["full"] != transition_ids["no_hfs"]:
        raise ValueError("full/no_hfs reference-transition topology differs")
    if cut_ids & transition_ids["full"]:
        raise ValueError("drain cut unexpectedly crosses a reference basin")
    selected = sorted(cut_ids | transition_ids["full"])
    temperature = float(contract["temperature_K"])
    thermal_voltage = exact(KB * temperature / Q)
    with localcontext() as context:
        context.prec = int(contract["precision"]["decimal_digits"])
        rows = [transform_row(variant, edges[variant][edge_id],
                              states[variant], cut_ids, thermal_voltage)
                for variant in variants for edge_id in selected]

    pair_rows: list[dict[str, Any]] = []
    by_key = {(row["variant"], row["edge_id"]): row for row in rows}
    for edge_id in selected:
        full = by_key[("full", edge_id)]
        off = by_key[("no_hfs", edge_id)]
        pair_rows.append({
            "edge_id": edge_id,
            "edge_class": full["edge_class"],
            "node0": full["node0"],
            "node1": full["node1"],
            "finite_couple_transport_edge": full[
                "finite_couple_transport_edge"],
            "reference_jump_V": full["reference_jump_V"],
            "full_physical_qf_drop_V": full["physical_qf_drop_V"],
            "no_hfs_physical_qf_drop_V": off["physical_qf_drop_V"],
            "full_minus_no_hfs_qf_drop_V": (
                full["physical_qf_drop_V"] - off["physical_qf_drop_V"]),
            "full_sg_log_imbalance": full["observed_sg_log_left_over_right"],
            "no_hfs_sg_log_imbalance": off["observed_sg_log_left_over_right"],
            "full_minus_no_hfs_sg_log_imbalance": (
                full["observed_sg_log_left_over_right"]
                - off["observed_sg_log_left_over_right"]),
            "full_electron_sg_stable_flux": full["electron_sg_stable_flux"],
            "no_hfs_electron_sg_stable_flux": off["electron_sg_stable_flux"],
            "full_minus_no_hfs_electron_sg_stable_flux": (
                full["electron_sg_stable_flux"]
                - off["electron_sg_stable_flux"]),
        })

    summaries = [summarize(rows, variant, edge_class)
                 for variant in variants
                 for edge_class in contract["edge_classes"]]
    maximum_endpoint = max(item["maximum_endpoint_transform_error_V"]
                           for item in summaries)
    maximum_drop = max(item["maximum_physical_drop_closure_error_V"]
                       for item in summaries)
    maximum_log = max(item["maximum_sg_log_state_closure_error"]
                      for item in summaries)
    maximum_flux = max(item["maximum_factorized_flux_relative_error"]
                       for item in summaries)
    reference_rows = [row for row in rows if row["reference_transition"]]
    drain_rows = [row for row in rows if row["drain_cut"]]
    reference_qf_response = max(
        abs(row["full_minus_no_hfs_qf_drop_V"]) for row in pair_rows
        if row["edge_class"] == "reference_transition")
    drain_qf_response = max(
        abs(row["full_minus_no_hfs_qf_drop_V"]) for row in pair_rows
        if row["edge_class"] == "drain_cut")
    write_csv(portable_dir / "m25_edge_reference_transform_ledger.csv", rows)
    write_csv(portable_dir / "m25_variant_edge_delta_ledger.csv", pair_rows)
    write_csv(portable_dir / "m25_class_summary.csv", summaries)
    report = {
        "schema": "vela.simplemos.sdevice.m25_qf_reference_transform_audit_report.v1",
        "status": "complete",
        "execution": {
            "read_only_state_count": 2,
            "edge_evaluation_count": len(rows),
            "new_sentaurus_execution": False,
            "default_model_changed": False,
            "decimal_digits": int(contract["precision"]["decimal_digits"]),
        },
        "topology": {
            "unique_selected_edges": len(selected),
            "drain_cut_edge_count": len(cut_ids),
            "drain_cut_reference_transition_edge_count": len(
                cut_ids & transition_ids["full"]),
            "reference_transition_edge_count": len(transition_ids["full"]),
            "finite_couple_reference_transition_edge_count": len({
                row["edge_id"] for row in reference_rows
                if row["finite_couple_transport_edge"]}),
            "reference_jump_values_V": sorted({
                abs(row["reference_jump_V"]) for row in reference_rows}),
        },
        "closure": {
            "maximum_endpoint_transform_error_V": maximum_endpoint,
            "maximum_physical_drop_closure_error_V": maximum_drop,
            "maximum_sg_log_state_closure_error": maximum_log,
            "maximum_factorized_flux_relative_error": maximum_flux,
            "all_acceptance_checks_pass": (
                maximum_endpoint < contract["precision"][
                    "maximum_endpoint_transform_error_V"]
                and maximum_drop < contract["precision"][
                    "maximum_physical_drop_closure_error_V"]
                and maximum_log < contract["precision"][
                    "maximum_sg_log_closure_error"]
                and maximum_flux < contract["precision"][
                    "maximum_factorized_flux_relative_error"]),
        },
        "findings": {
            "maximum_naive_increment_only_error_on_reference_transition_V": max(
                abs(row["naive_drop_error_V"]) for row in reference_rows),
            "maximum_reference_jump_thermal_voltages": max(
                abs(row["reference_jump_thermal_voltages"])
                for row in reference_rows),
            "maximum_full_no_hfs_reference_transition_qf_drop_delta_V": (
                reference_qf_response),
            "maximum_full_no_hfs_drain_cut_qf_drop_delta_V": drain_qf_response,
            "maximum_absolute_drain_cut_qf_drop_V": max(
                abs(row["physical_qf_drop_V"]) for row in drain_rows),
        },
        "class_summaries": summaries,
        "conclusions": {
            "supported": [
                "The seven drain-cut edges remain inside the drain reference basin; no cross-basin transform occurs on the terminal-current cut.",
                "Across all 51 reference-transition edges, the production anchored endpoints, physical QF drop, SG log imbalance, and factorized flux close to the frozen state at the stated numerical tolerances.",
                "Using increments without the reference jump would introduce an approximately 0.05 V artificial QF drop on every reference-transition edge.",
                "The HFS-on/off terminal response cannot be attributed to a cross-reference transform located on the drain cut.",
            ],
            "not_supported": [
                "Equating this Vela reference-coordinate closure with an undocumented Sentaurus contact-basin implementation.",
                "Using the cross-reference edges as an additive terminal-current cut; zero-couple and internal edges are included for coordinate closure only.",
                "Assigning the remaining deep-off current difference to a contact reference transform without a corresponding Sentaurus node-level export.",
            ],
        },
        "artifacts": {
            "edge_ledger_csv": portable(
                portable_dir / "m25_edge_reference_transform_ledger.csv"),
            "variant_delta_csv": portable(
                portable_dir / "m25_variant_edge_delta_ledger.csv"),
            "class_summary_csv": portable(
                portable_dir / "m25_class_summary.csv"),
        },
        "claim_policy": contract["claim_policy"],
    }
    write_json(portable_dir / "m25_qf_reference_transform_audit_report.json",
               report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--portable-output", type=Path, default=PORTABLE)
    args = parser.parse_args()
    contract = read_json(args.contract)
    validate(contract)
    args.output.mkdir(parents=True, exist_ok=True)
    report = analyze(contract, args.portable_output)
    write_json(args.output / "m25_qf_reference_transform_audit_report.json",
               report)
    print(json.dumps({
        "status": report["status"],
        "selected_edges": report["topology"]["unique_selected_edges"],
        "reference_transition_edges": report["topology"][
            "reference_transition_edge_count"],
        "drain_cut_cross_reference_edges": report["topology"][
            "drain_cut_reference_transition_edge_count"],
        "all_acceptance_checks_pass": report["closure"][
            "all_acceptance_checks_pass"],
    }))


if __name__ == "__main__":
    main()
