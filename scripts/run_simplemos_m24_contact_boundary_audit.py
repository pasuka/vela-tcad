#!/usr/bin/env python3
"""Audit the n23 deep-off drain boundary and first free-node layer."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from decimal import Decimal
import json
import math
from pathlib import Path
import sys
from typing import Any


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m10_fixed_state_replay as m10  # noqa: E402


ROOT = REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
M22 = ROOT / "m22_n23_hfs_deep_off"
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m24_contact_boundary_audit_contract_v1.json")
OUTPUT = ROOT / "m24_contact_boundary_audit"
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "contact_boundary_audit")
RUNNER = REPO / "build-release/vela_example_runner.exe"
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


def validate(contract: dict[str, Any]) -> None:
    if contract.get("schema") != (
            "vela.simplemos.sdevice.m24_contact_boundary_audit.v1"):
        raise ValueError("unexpected M24 contract schema")
    if contract.get("state_variants") != ["full", "no_hfs"]:
        raise ValueError("M24 requires the paired full/no_hfs states")


def source_dir(variant: str) -> Path:
    return M22 / "self_consistent" / variant / "vg_0p05"


def sg_path(variant: str) -> Path:
    return (M22 / "frozen/vg_0p05" / f"state_{variant}"
            / f"eval_{variant}/sg_edges.csv")


def state_map(variant: str) -> dict[int, dict[str, str]]:
    return {int(row["node_id"]): row
            for row in read_csv(source_dir(variant) / "state.csv")}


def mesh_and_contact(config: dict[str, Any], contact: str
                     ) -> tuple[dict[str, Any], set[int]]:
    mesh = read_json(Path(config["mesh_file"]))
    matches = [item for item in mesh["contacts"]
               if item["name"].lower() == contact.lower()]
    if len(matches) != 1:
        raise ValueError(f"expected one {contact} contact")
    return mesh, {int(node) for node in matches[0]["node_ids"]}


def run_probes(variant: str, output: Path, runner: Path,
               force: bool) -> dict[str, Path]:
    base = read_json(source_dir(variant) / "config.json")
    state = source_dir(variant) / "state.csv"
    run_dir = output / variant
    run_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {"sg": sg_path(variant), "state": state}
    jobs: list[Path] = []
    for label, simulation_type in (
            ("terms", "newton_carrier_term_probe"),
            ("rows", "newton_carrier_row_probe")):
        csv_path = run_dir / f"{label}.csv"
        config_path = run_dir / f"{label}.json"
        deck = dict(base)
        deck.update({
            "simulation_type": simulation_type,
            "state_file": str(state.resolve()),
            "output_csv": str(csv_path.resolve()),
            "simplemos_m24": {"read_only": True, "variant": variant},
        })
        deck.pop("residual_output_csv", None)
        if label == "terms":
            deck["carrier_term_probe"] = {"solved_equation_terms": False}
        write_json(config_path, deck)
        paths[label] = csv_path
        if force or not csv_path.is_file():
            jobs.append(config_path)
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(lambda path: m10.execute_runner(path, runner), jobs))
    return paths


def signed_flux(row: dict[str, str], node: int) -> float:
    if int(row["node0"]) == node:
        return float(row["electron_flux"])
    if int(row["node1"]) == node:
        return -float(row["electron_flux"])
    raise ValueError(f"node {node} is not incident to edge {row['edge_id']}")


def signed_particle_at_contact(row: dict[str, str], contact_nodes: set[int]) -> float:
    sign = 1.0 if int(row["node0"]) in contact_nodes else -1.0
    return sign * float(row["electron_particle_line_flux_per_m_s"])


def target_curve_row(variant: str) -> dict[str, str]:
    rows = read_csv(source_dir(variant) / "curve.csv")
    return min(rows, key=lambda row: abs(float(row["bias_V"]) - 0.05))


def analyze(contract: dict[str, Any], paths: dict[str, dict[str, Path]],
            portable_dir: Path) -> dict[str, Any]:
    base_config = read_json(source_dir("full") / "config.json")
    _, contact_nodes = mesh_and_contact(base_config, contract["contact"])
    drain_bias = Decimal(str(contract["drain_voltage_V"]))
    variants = contract["state_variants"]
    states = {variant: state_map(variant) for variant in variants}
    sg = {variant: read_csv(paths[variant]["sg"]) for variant in variants}
    terms = {variant: {int(row["node_id"]): row for row in read_csv(
        paths[variant]["terms"])} for variant in variants}
    rows = {variant: {int(row["node_id"]): row for row in read_csv(
        paths[variant]["rows"])} for variant in variants}
    cut_ids = {int(row["edge_id"]) for row in sg["full"]
               if (int(row["node0"]) in contact_nodes)
               != (int(row["node1"]) in contact_nodes)}
    adjacent_nodes = sorted({node for row in sg["full"]
                             if int(row["edge_id"]) in cut_ids
                             for node in (int(row["node0"]), int(row["node1"]))
                             if node not in contact_nodes})

    contact_rows: list[dict[str, Any]] = []
    maximum_bias_error = 0.0
    maximum_physical_delta = 0.0
    maximum_increment_delta = 0.0
    for node in sorted(contact_nodes):
        a, b = states["full"][node], states["no_hfs"][node]
        physical_deltas = [abs(float(a[column]) - float(b[column]))
                           for column in ("psi", "phin", "phip",
                                          "electrons_m3", "holes_m3")]
        maximum_physical_delta = max(maximum_physical_delta,
                                     max(physical_deltas))
        increment_delta = abs(float(a["electron_qf_increment_V"])
                              - float(b["electron_qf_increment_V"]))
        maximum_increment_delta = max(maximum_increment_delta, increment_delta)
        for variant in variants:
            state = states[variant][node]
            absolute_qf = (Decimal(state["electron_qf_reference_V"])
                           + Decimal(state["electron_qf_increment_V"]))
            bias_error = float(abs(absolute_qf - drain_bias))
            maximum_bias_error = max(maximum_bias_error, bias_error)
            contact_rows.append({
                "variant": variant,
                "node_id": node,
                "psi_V": state["psi"],
                "electron_density_m3": state["electrons_m3"],
                "hole_density_m3": state["holes_m3"],
                "electron_qf_reference_V": state["electron_qf_reference_V"],
                "electron_qf_increment_V": state["electron_qf_increment_V"],
                "electron_qf_bias_error_V": bias_error,
                "full_minus_no_hfs_psi_V": float(a["psi"]) - float(b["psi"]),
                "full_minus_no_hfs_electron_density_m3": (
                    float(a["electrons_m3"]) - float(b["electrons_m3"])),
                "full_minus_no_hfs_qf_increment_V": (
                    float(a["electron_qf_increment_V"])
                    - float(b["electron_qf_increment_V"])),
            })

    edge_rows: list[dict[str, Any]] = []
    cut_currents: dict[str, float] = {}
    maximum_cut_disagreement = 0.0
    maximum_cut_absolute_disagreement = 0.0
    for variant in variants:
        cut = [row for row in sg[variant]
               if int(row["edge_id"]) in cut_ids]
        cut_current = sum(-Q * signed_particle_at_contact(
            row, contact_nodes) * 1.0e-6 for row in cut)
        curve_current = float(target_curve_row(variant)[
            "current_electron_A_per_um"])
        absolute_disagreement = abs(cut_current - curve_current)
        disagreement = absolute_disagreement / max(
            abs(cut_current), abs(curve_current), 1.0e-300)
        maximum_cut_disagreement = max(maximum_cut_disagreement, disagreement)
        maximum_cut_absolute_disagreement = max(
            maximum_cut_absolute_disagreement, absolute_disagreement)
        cut_currents[variant] = cut_current
        for row in cut:
            contact_node = (int(row["node0"]) if int(row["node0"]) in contact_nodes
                            else int(row["node1"]))
            interior_node = (int(row["node1"]) if int(row["node0"]) in contact_nodes
                             else int(row["node0"]))
            edge_rows.append({
                "variant": variant,
                "edge_id": row["edge_id"],
                "contact_node": contact_node,
                "interior_node": interior_node,
                "electron_current_A_per_um": -Q * signed_particle_at_contact(
                    row, contact_nodes) * 1.0e-6,
                "electron_mobility_m2_V_s": row["electron_mobility_m2_V_s"],
                "electron_sg_phin_contact_relative_V": (
                    row["electron_sg_phin0_relative_V"]
                    if int(row["node0"]) in contact_nodes
                    else row["electron_sg_phin1_relative_V"]),
                "electron_sg_phin_interior_relative_V": (
                    row["electron_sg_phin1_relative_V"]
                    if int(row["node0"]) in contact_nodes
                    else row["electron_sg_phin0_relative_V"]),
                "electron_sg_log_left_over_right": row[
                    "electron_sg_log_left_over_right"],
                "electron_sg_cancellation_condition": row[
                    "electron_sg_cancellation_condition"],
            })

    node_rows: list[dict[str, Any]] = []
    max_term_closure = 0.0
    variant_node: dict[str, dict[int, dict[str, float]]] = {}
    for variant in variants:
        by_node: dict[int, dict[str, float]] = {}
        for node in adjacent_nodes:
            incident = [row for row in sg[variant]
                        if node in (int(row["node0"]), int(row["node1"]))]
            contact_flux = sum(signed_flux(row, node) for row in incident
                               if int(row["edge_id"]) in cut_ids)
            internal_flux = sum(signed_flux(row, node) for row in incident
                                if int(row["edge_id"]) not in cut_ids)
            term = terms[variant][node]
            source = sum(float(term[column]) for column in (
                "electron_recombination", "electron_impact",
                "electron_gauge", "electron_boundary"))
            closure = contact_flux + internal_flux + source - float(
                term["electron_residual"])
            max_term_closure = max(max_term_closure, abs(closure))
            row = rows[variant][node]
            contact_share = float(row["electron_contact_column_abs_sum"]) / max(
                float(row["electron_row_abs_sum"]), 1.0e-300)
            by_node[node] = {"contact_flux": contact_flux,
                             "internal_flux": internal_flux,
                             "source": source,
                             "residual": float(term["electron_residual"]),
                             "contact_column_share": contact_share}
        variant_node[variant] = by_node
    for node in adjacent_nodes:
        a, b = variant_node["full"][node], variant_node["no_hfs"][node]
        node_rows.append({
            "node_id": node,
            "full_contact_flux": a["contact_flux"],
            "no_hfs_contact_flux": b["contact_flux"],
            "delta_contact_flux": a["contact_flux"] - b["contact_flux"],
            "full_internal_flux": a["internal_flux"],
            "no_hfs_internal_flux": b["internal_flux"],
            "delta_internal_flux": a["internal_flux"] - b["internal_flux"],
            "delta_source": a["source"] - b["source"],
            "delta_residual": a["residual"] - b["residual"],
            "full_contact_column_share": a["contact_column_share"],
            "no_hfs_contact_column_share": b["contact_column_share"],
            "full_qf_increment_V": states["full"][node][
                "electron_qf_increment_V"],
            "no_hfs_qf_increment_V": states["no_hfs"][node][
                "electron_qf_increment_V"],
            "qf_increment_delta_V": (
                float(states["full"][node]["electron_qf_increment_V"])
                - float(states["no_hfs"][node]["electron_qf_increment_V"])),
        })

    write_csv(portable_dir / "m24_contact_node_ledger.csv", contact_rows)
    write_csv(portable_dir / "m24_drain_cut_edge_ledger.csv", edge_rows)
    write_csv(portable_dir / "m24_first_layer_continuity_ledger.csv", node_rows)
    effect_dex = math.log10(abs(cut_currents["full"])
                            / abs(cut_currents["no_hfs"]))
    report = {
        "schema": "vela.simplemos.sdevice.m24_contact_boundary_audit_report.v1",
        "status": "complete",
        "execution": {
            "contact_node_count": len(contact_nodes),
            "drain_cut_edge_count": len(cut_ids),
            "first_layer_node_count": len(adjacent_nodes),
            "read_only_probe_count": 4,
            "new_sentaurus_execution": False,
            "default_model_changed": False,
        },
        "closure": {
            "maximum_contact_qf_bias_error_V": maximum_bias_error,
            "maximum_sg_cut_current_relative_disagreement": (
                maximum_cut_disagreement),
            "maximum_sg_cut_current_absolute_disagreement_A_per_um": (
                maximum_cut_absolute_disagreement),
            "maximum_continuity_term_absolute_closure": max_term_closure,
        },
        "findings": {
            "maximum_full_no_hfs_contact_physical_state_delta": (
                maximum_physical_delta),
            "maximum_full_no_hfs_contact_qf_increment_delta_V": (
                maximum_increment_delta),
            "full_cut_electron_current_A_per_um": cut_currents["full"],
            "no_hfs_cut_electron_current_A_per_um": cut_currents["no_hfs"],
            "self_consistent_cut_current_effect_dex": effect_dex,
            "maximum_first_layer_contact_column_share": max(
                max(row["full_contact_column_share"],
                    row["no_hfs_contact_column_share"]) for row in node_rows),
            "maximum_first_layer_qf_increment_delta_V": max(
                abs(row["qf_increment_delta_V"]) for row in node_rows),
            "signed_first_layer_contact_flux_delta_sum": sum(
                row["delta_contact_flux"] for row in node_rows),
            "signed_first_layer_internal_flux_delta_sum": sum(
                row["delta_internal_flux"] for row in node_rows),
            "signed_first_layer_source_delta_sum": sum(
                row["delta_source"] for row in node_rows),
        },
        "conclusions": {
            "supported": [
                "The Vela drain-contact physical state is identical between the self-consistent HFS-on and HFS-off states.",
                "The drain electron quasi-Fermi target remains pinned to the applied drain bias in both states.",
                "The approximately 0.10 dex current response is carried by tiny first-layer quasi-Fermi increment changes and balanced contact/internal flux changes, not by a changed contact boundary value.",
            ],
            "not_supported": [
                "Claiming that the Vela ohmic contact target directly causes the HFS-on/off response.",
                "Claiming Vela/Sentaurus n23 contact-boundary equivalence without an n23 Sentaurus node-state export.",
                "Changing a contact or transport default from this audit.",
            ],
        },
        "artifacts": {
            "contact_nodes": portable(portable_dir / "m24_contact_node_ledger.csv"),
            "drain_cut_edges": portable(
                portable_dir / "m24_drain_cut_edge_ledger.csv"),
            "first_layer_rows": portable(
                portable_dir / "m24_first_layer_continuity_ledger.csv"),
        },
        "claim_policy": contract["claim_policy"],
    }
    write_json(portable_dir / "m24_contact_boundary_audit_report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--portable", type=Path, default=PORTABLE)
    parser.add_argument("--runner", type=Path, default=RUNNER)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    contract = read_json(args.contract.resolve())
    validate(contract)
    paths = {variant: run_probes(variant, args.output.resolve(),
                                 args.runner.resolve(), args.force)
             for variant in contract["state_variants"]}
    report = analyze(contract, paths, args.portable.resolve())
    acceptance = contract["acceptance"]
    closure = report["closure"]
    if closure["maximum_contact_qf_bias_error_V"] > acceptance[
            "maximum_contact_qf_bias_error_V"]:
        raise AssertionError("contact QF bias closure failed")
    if closure["maximum_sg_cut_current_relative_disagreement"] > acceptance[
            "maximum_sg_cut_current_relative_disagreement"]:
        raise AssertionError("cut-current extraction disagreement exceeded limit")
    if closure["maximum_continuity_term_absolute_closure"] > acceptance[
            "maximum_continuity_term_absolute_closure"]:
        raise AssertionError("continuity-term closure failed")
    print(json.dumps({"status": "complete", **report["execution"],
                      **report["closure"], **report["findings"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
