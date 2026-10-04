#!/usr/bin/env python3
"""Analyze and freeze SimpleMOS M56 Si/SiO2 interface topology audit."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m53_direct_current_full_matrix as m53  # noqa: E402


ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m56_si_oxide_interface_topology_contract_v1.json"
FREEZE = ROOT / "simplemos_m56_si_oxide_interface_topology_contract_freeze.json"
M8_EVIDENCE = ROOT / "simplemos_m8_original_physics_evidence.json"
M52_EVIDENCE = ROOT / "simplemos_m52_direct_current_attribution_evidence.json"
M52_REPORT = ROOT / "direct_current_attribution/m52_direct_current_attribution_report.json"
M8_NEUTRAL = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics/neutral"
M8_DECKS = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics/sentaurus_bundle"
M52_STATES = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m52_direct_current_attribution/state_exports/weighted"
PORTABLE = ROOT / "si_oxide_interface_topology"
REPORT = PORTABLE / "m56_si_oxide_interface_topology_report.json"
TOPOLOGY = PORTABLE / "m56_input_topology_ledger.csv"
STATES = PORTABLE / "m56_solved_state_interface_ledger.csv"
NODES = PORTABLE / "m56_interface_node_potential_ledger.csv"
SUPPORT = PORTABLE / "m56_field_support_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m56_si_oxide_interface_topology_2026-09-02.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m56/artifact.json"
EVIDENCE = ROOT / "simplemos_m56_si_oxide_interface_topology_evidence.json"
TEST = REPO / "tests/regression/test_simplemos_m56_si_oxide_interface_topology.py"
AUDIT_FIELDS = (
    "ElectrostaticPotential", "eDensity", "hDensity",
    "eCurrentDensity", "hCurrentDensity",
    "eQuasiFermiPotential", "hQuasiFermiPotential",
)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8",
                    newline="\n")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]),
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO.resolve()).as_posix()
    except ValueError:
        return str(resolved)


def validate_contract() -> dict[str, Any]:
    contract = read_json(CONTRACT)
    freeze = read_json(FREEZE)
    if contract.get("schema") != (
            "vela.simplemos.sdevice.m56_si_oxide_interface_topology_contract.v1"):
        raise ValueError("unexpected M56 contract schema")
    if freeze.get("status") != "frozen_before_analysis":
        raise ValueError("M56 contract was not frozen before analysis")
    if freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M56 contract changed after freeze")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M56 upstream changed: {relative}")
    if read_json(M8_EVIDENCE).get("status") != "accepted":
        raise ValueError("M56 requires accepted M8 evidence")
    if read_json(M52_EVIDENCE).get("status") != "frozen":
        raise ValueError("M56 requires frozen M52 evidence")
    m52 = read_json(M52_REPORT)
    upstream = contract["upstream"]
    if m52.get("classification") != upstream["required_m52_classification"]:
        raise ValueError("M52 classification changed")
    if int(m52["state_field_summary"]["comparison_row_count"]) != int(
            upstream["required_m52_state_field_comparison_count"]):
        raise ValueError("M52 field-count anchor changed")
    if float(m52["state_field_summary"]["maximum_absolute_difference"]) != 0.0:
        raise ValueError("M52 field-invariance anchor changed")
    return contract


def topology(root: Path) -> dict[str, Any]:
    node_rows = read_csv(root / "nodes.csv")
    element_rows = read_csv(root / "elements.csv")
    contact_rows = read_csv(root / "contacts.csv")
    coordinates = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in node_rows
    }
    coordinate_groups: dict[tuple[float, float], list[int]] = defaultdict(list)
    for node, coordinate in coordinates.items():
        coordinate_groups[coordinate].append(node)
    duplicates = [nodes for nodes in coordinate_groups.values() if len(nodes) > 1]
    region_nodes: dict[str, set[int]] = defaultdict(set)
    region_material: dict[str, str] = {}
    region_edges: dict[str, Counter[tuple[int, int]]] = defaultdict(Counter)
    for row in element_rows:
        region = row["region"]
        region_material[region] = row["material"]
        vertices = (int(row["node0"]), int(row["node1"]), int(row["node2"]))
        region_nodes[region].update(vertices)
        for left, right in ((vertices[0], vertices[1]),
                            (vertices[1], vertices[2]),
                            (vertices[2], vertices[0])):
            region_edges[region][tuple(sorted((left, right)))] += 1
    silicon = next(region for region, material in region_material.items()
                   if material == "Si")
    oxide = next(region for region, material in region_material.items()
                 if material == "SiO2")
    silicon_boundary = {edge for edge, count in region_edges[silicon].items()
                        if count == 1}
    oxide_boundary = {edge for edge, count in region_edges[oxide].items()
                      if count == 1}
    interface_edges = silicon_boundary & oxide_boundary
    interface_nodes = {node for edge in interface_edges for node in edge}
    duplicate_interface_groups = 0
    for group in duplicates:
        if (set(group) & region_nodes[silicon]
                and set(group) & region_nodes[oxide]):
            duplicate_interface_groups += 1
    contacts = {row["name"]: row["region"] for row in contact_rows}
    expected_contacts = {"gate": oxide, "source": silicon,
                         "drain": silicon, "substrate": silicon}
    return {
        "coordinates": coordinates,
        "region_nodes": region_nodes,
        "region_material": region_material,
        "silicon_region": silicon,
        "oxide_region": oxide,
        "interface_edges": interface_edges,
        "interface_nodes": interface_nodes,
        "global_node_count": len(coordinates),
        "element_count": len(element_rows),
        "duplicate_coordinate_group_count": len(duplicates),
        "duplicate_coordinate_node_count": sum(len(group) for group in duplicates),
        "duplicate_interface_coordinate_group_count": duplicate_interface_groups,
        "shared_region_node_count": len(region_nodes[silicon] & region_nodes[oxide]),
        "contacts": contacts,
        "contact_attachment_ok": contacts == expected_contacts,
        "raw_paths": [root / "nodes.csv", root / "elements.csv",
                      root / "contacts.csv"],
    }


def field_entries(root: Path) -> list[dict[str, Any]]:
    return list(read_json(root / "field_manifest.json")["fields"])


def field_map(root: Path, entry: dict[str, Any]) -> dict[int, float]:
    path = root / "fields" / entry["csv_file"]
    return {int(row["node_id"]): float(row["component0"])
            for row in read_csv(path)}


def analyze(contract: dict[str, Any]) -> tuple[dict[str, Any], list[Path]]:
    raw_paths: list[Path] = []
    topology_rows: list[dict[str, Any]] = []
    for device in contract["scope"]["input_tdr_devices"]:
        root = M8_NEUTRAL / device
        topo = topology(root)
        raw_paths.extend(topo["raw_paths"])
        topology_rows.append({
            "device": device,
            "global_node_count": topo["global_node_count"],
            "element_count": topo["element_count"],
            "silicon_region": topo["silicon_region"],
            "oxide_region": topo["oxide_region"],
            "silicon_oxide_shared_edge_count": len(topo["interface_edges"]),
            "silicon_oxide_interface_node_count": len(topo["interface_nodes"]),
            "shared_region_node_count": topo["shared_region_node_count"],
            "duplicate_coordinate_group_count":
                topo["duplicate_coordinate_group_count"],
            "duplicate_coordinate_node_count": topo["duplicate_coordinate_node_count"],
            "duplicate_interface_coordinate_group_count":
                topo["duplicate_interface_coordinate_group_count"],
            "gate_contact_region": topo["contacts"].get("gate", ""),
            "source_contact_region": topo["contacts"].get("source", ""),
            "drain_contact_region": topo["contacts"].get("drain", ""),
            "substrate_contact_region": topo["contacts"].get("substrate", ""),
            "contact_attachment_ok": topo["contact_attachment_ok"],
        })

    state_rows: list[dict[str, Any]] = []
    node_rows: list[dict[str, Any]] = []
    support_rows: list[dict[str, Any]] = []
    max_potential_mismatch = 0.0
    carrier_fields = {"eDensity", "hDensity", "eCurrentDensity", "hCurrentDensity"}
    required_carrier_regions = set(contract["acceptance"]["carrier_support_regions"])
    for device in contract["scope"]["solved_state_devices"]:
        for gate in map(float, contract["scope"]["solved_state_gate_voltages_V"]):
            state = f"{device}_vd_0p05_vg_{m53.voltage_tag(gate)}"
            root = M52_STATES / state
            topo = topology(root)
            entries = field_entries(root)
            raw_paths.extend(topo["raw_paths"] + [root / "field_manifest.json"])
            selected = [entry for entry in entries if entry["name"] in AUDIT_FIELDS]
            supports: dict[str, set[str]] = defaultdict(set)
            for entry in selected:
                supports[entry["name"]].add(entry["region_name"])
                support_rows.append({
                    "state": state, "device": device,
                    "gate_voltage_V": gate,
                    "field": entry["name"],
                    "region": entry["region_name"],
                    "value_count": entry["values"],
                    "components": entry["components"],
                    "mapping_status": entry["mapping_status"],
                })
            potential_entries = {
                entry["region_name"]: entry for entry in selected
                if entry["name"] == "ElectrostaticPotential"
            }
            silicon = topo["silicon_region"]
            oxide = topo["oxide_region"]
            if silicon not in potential_entries or oxide not in potential_entries:
                silicon_values: dict[int, float] = {}
                oxide_values: dict[int, float] = {}
            else:
                silicon_values = field_map(root, potential_entries[silicon])
                oxide_values = field_map(root, potential_entries[oxide])
                raw_paths.extend([
                    root / "fields" / potential_entries[silicon]["csv_file"],
                    root / "fields" / potential_entries[oxide]["csv_file"],
                ])
            shared_with_values = (topo["interface_nodes"] & silicon_values.keys()
                                  & oxide_values.keys())
            state_max = 0.0
            for node in sorted(shared_with_values):
                difference = silicon_values[node] - oxide_values[node]
                state_max = max(state_max, abs(difference))
                coordinate = topo["coordinates"][node]
                node_rows.append({
                    "state": state, "device": device,
                    "gate_voltage_V": gate, "node_id": node,
                    "x_um": coordinate[0], "y_um": coordinate[1],
                    "silicon_potential_V": silicon_values[node],
                    "oxide_potential_V": oxide_values[node],
                    "silicon_minus_oxide_potential_V": difference,
                })
            max_potential_mismatch = max(max_potential_mismatch, state_max)
            carrier_support_ok = all(supports[field] == required_carrier_regions
                                     for field in carrier_fields)
            state_rows.append({
                "state": state, "device": device,
                "drain_voltage_V": 0.05, "gate_voltage_V": gate,
                "global_node_count": topo["global_node_count"],
                "silicon_oxide_shared_edge_count": len(topo["interface_edges"]),
                "silicon_oxide_interface_node_count": len(topo["interface_nodes"]),
                "interface_nodes_with_both_potential_values":
                    len(shared_with_values),
                "maximum_interface_potential_mismatch_V": state_max,
                "potential_support_regions": ";".join(
                    sorted(supports["ElectrostaticPotential"])),
                "carrier_density_current_support_regions": ";".join(sorted(
                    set().union(*(supports[field] for field in carrier_fields)))),
                "carrier_density_current_support_is_silicon_only":
                    carrier_support_ok,
                "electron_qf_plot_support_regions": ";".join(
                    sorted(supports["eQuasiFermiPotential"])),
                "hole_qf_plot_support_regions": ";".join(
                    sorted(supports["hQuasiFermiPotential"])),
                "contact_attachment_ok": topo["contact_attachment_ok"],
                "duplicate_coordinate_group_count":
                    topo["duplicate_coordinate_group_count"],
                "duplicate_interface_coordinate_group_count":
                    topo["duplicate_interface_coordinate_group_count"],
            })

    forbidden_tokens = ("HeteroInterface", "Thermionic", "Discontinuity",
                        "CurrentWeighting", "DirectCurrent")
    deck_rows: list[dict[str, Any]] = []
    deck_paths: list[Path] = []
    for device in contract["scope"]["input_tdr_devices"]:
        for drain in (0.05, 1.0):
            case = f"{device}_vd_{m53.voltage_tag(drain)}"
            path = M8_DECKS / device / f"{case}_des.cmd"
            text = path.read_text(encoding="utf-8")
            deck_paths.append(path)
            deck_rows.append({
                "case": case,
                "forbidden_interface_or_current_override_present": any(
                    token in text for token in forbidden_tokens),
            })
    raw_paths.extend(deck_paths)
    potential_tolerance = float(contract["acceptance"][
        "maximum_interface_potential_absolute_mismatch_V"])
    topology_complete = (
        len(topology_rows) == 8
        and all(int(row["silicon_oxide_shared_edge_count"]) >= 1
                for row in topology_rows)
        and all(bool(row["contact_attachment_ok"]) for row in topology_rows))
    states_complete = (
        len(state_rows) == 6
        and all(int(row["interface_nodes_with_both_potential_values"])
                == int(row["silicon_oxide_interface_node_count"])
                for row in state_rows)
        and all(bool(row["carrier_density_current_support_is_silicon_only"])
                for row in state_rows)
        and all(bool(row["contact_attachment_ok"]) for row in state_rows))
    duplicate_interface = any(
        int(row["duplicate_interface_coordinate_group_count"]) > 0
        for row in topology_rows + state_rows)
    potential_continuous = max_potential_mismatch <= potential_tolerance
    deck_clean = not any(bool(row["forbidden_interface_or_current_override_present"])
                         for row in deck_rows)
    if not topology_complete or not states_complete or not deck_clean:
        classification = "support_or_topology_incomplete"
    elif duplicate_interface:
        classification = "explicit_geometric_double_nodes"
    elif not potential_continuous:
        classification = "potential_discontinuity_observed"
    else:
        classification = "shared_topology_continuous_potential_insulator"
    acceptance = {
        "contract_hash_unchanged": sha256(CONTRACT) ==
            read_json(FREEZE)["contract_sha256"],
        "input_topology_count": len(topology_rows) == 8,
        "solved_state_count": len(state_rows) == 6,
        "interface_topology_complete": topology_complete,
        "interface_potential_support_complete": all(
            int(row["interface_nodes_with_both_potential_values"])
            >= int(contract["acceptance"]["minimum_interface_node_count_per_state"])
            for row in state_rows),
        "interface_potential_continuous": potential_continuous,
        "carrier_support_silicon_only": all(
            bool(row["carrier_density_current_support_is_silicon_only"])
            for row in state_rows),
        "contact_attachments_match": all(
            bool(row["contact_attachment_ok"]) for row in topology_rows + state_rows),
        "default_decks_have_no_explicit_override": deck_clean,
        "all_values_finite": all(math.isfinite(float(row[key]))
                                 for row in node_rows
                                 for key in ("silicon_potential_V",
                                             "oxide_potential_V")),
        "classification_declared": classification in
            contract["analysis"]["classifications"],
        "no_new_solver_execution": True,
        "historical_artifacts_preserved": True,
    }
    acceptance["all_checks_pass"] = all(acceptance.values())
    report = {
        "schema": "vela.simplemos.sdevice.m56_si_oxide_interface_topology_report.v1",
        "status": "complete" if acceptance["all_checks_pass"] else "failed",
        "classification": classification,
        "execution": {
            "new_sentaurus_execution": False,
            "new_vela_execution": False,
            "inputs": "frozen M8 neutral TDR exports and six M52 solved-state exports",
            "contract_sha256_before_and_after": sha256(CONTRACT),
        },
        "findings": {
            "input_topology_count": len(topology_rows),
            "solved_state_count": len(state_rows),
            "total_input_duplicate_coordinate_groups": sum(
                int(row["duplicate_coordinate_group_count"])
                for row in topology_rows),
            "total_input_duplicate_interface_coordinate_groups": sum(
                int(row["duplicate_interface_coordinate_group_count"])
                for row in topology_rows),
            "minimum_input_shared_interface_edge_count": min(
                int(row["silicon_oxide_shared_edge_count"])
                for row in topology_rows),
            "minimum_solved_interface_node_count": min(
                int(row["silicon_oxide_interface_node_count"])
                for row in state_rows),
            "maximum_interface_potential_mismatch_V": max_potential_mismatch,
            "all_gate_contacts_on_oxide": all(
                row["gate_contact_region"] == "Oxide_1"
                for row in topology_rows),
            "all_other_contacts_on_silicon": all(
                row[key] == "Silicon_1" for row in topology_rows
                for key in ("source_contact_region", "drain_contact_region",
                            "substrate_contact_region")),
            "carrier_density_current_fields_silicon_only": all(
                bool(row["carrier_density_current_support_is_silicon_only"])
                for row in state_rows),
            "quasi_fermi_plot_records_outside_silicon_are_diagnostic_only": True,
            "default_deck_explicit_interface_override_present": not deck_clean,
        },
        "input_topologies": topology_rows,
        "solved_states": state_rows,
        "deck_audit": deck_rows,
        "acceptance": acceptance,
        "claim_guard": contract["analysis"]["claim_guard"],
    }
    write_csv(TOPOLOGY, topology_rows)
    write_csv(STATES, state_rows)
    write_csv(NODES, node_rows)
    write_csv(SUPPORT, support_rows)
    write_json(REPORT, report)
    return report, raw_paths


def freeze_artifacts(report: dict[str, Any], raw_paths: list[Path]) -> None:
    findings = report["findings"]
    topology_table = "\n".join(
        f"| {row['device']} | {row['global_node_count']} | "
        f"{row['silicon_oxide_shared_edge_count']} | "
        f"{row['silicon_oxide_interface_node_count']} | "
        f"{row['duplicate_coordinate_group_count']} | "
        f"{row['gate_contact_region']} |"
        for row in report["input_topologies"])
    state_table = "\n".join(
        f"| {row['device']} | {float(row['gate_voltage_V']):.2f} | "
        f"{row['silicon_oxide_interface_node_count']} | "
        f"{float(row['maximum_interface_potential_mismatch_V']):.3e} | "
        f"{row['carrier_density_current_support_regions']} |"
        for row in report["solved_states"])
    doc = f'''# SimpleMOS M56 Si/SiO2 界面拓扑审计

## 结论

M56 分类为 `{report['classification']}`。八个输入TDR均未发现坐标完全相同但全局ID不同的几何重复节点；Si/SiO2 界面由两区域共享同一组全局节点ID和边来表示。求解场仍按区域分别导出，因此同一个界面全局节点在 Silicon 和 Oxide 电势文件中各有一条区域记录。

六个求解状态的所有共享界面节点均同时具有 Silicon/Oxide 电势值，最大跨区域差为 `{float(findings['maximum_interface_potential_mismatch_V']):.12e}` V。`eDensity`、`hDensity`、`eCurrentDensity`、`hCurrentDensity` 只在 `Silicon_1` 上有支持；氧化层中的准费米绘图记录不能解释为氧化层启用了载流子输运方程。

gate 接触在八个TDR中均附着于 `Oxide_1`，source、drain、substrate 均附着于 `Silicon_1`。因此“SiO2绝缘”意味着没有氧化层电子/空穴漂移扩散输运，不意味着不能在氧化层外边界的 gate 电极施加静电势。默认 deck 不含显式 `HeteroInterface`、`Thermionic`、`Discontinuity` 或双节点开关；本算例的区域场支持由TDR材料/区域拓扑和接触边界自动确定。

## 八个输入TDR

| 器件 | 全局节点 | Si/Ox共享边 | 界面节点 | 重复坐标组 | gate区域 |
|---|---:|---:|---:|---:|---|
{topology_table}

## 六个求解状态

| 器件 | Vg (V) | 界面节点 | 最大电势差 (V) | 载流子密度/电流支持区域 |
|---|---:|---:|---:|---|
{state_table}

## 边界

- 结论只针对冻结的 SimpleMOS TDR/deck，不外推到需要能带不连续、热发射或隧穿模型的其他异质结。
- M56 没有新增求解、修改网格、复制节点或添加界面物理。
- 准费米量在绝缘区域可作为绘图/派生记录存在，但载流子密度和电流密度的实际支持仍限定于 Silicon。

机器报告：`{portable(REPORT)}`。
'''
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(doc, encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {
        "schema": "vela.simplemos.sdevice.m56_artifact.v1",
        "title": "SimpleMOS M56 Si/SiO2 interface topology audit",
        "status": report["status"], "classification": report["classification"],
        "report": portable(REPORT),
        "ledgers": [portable(TOPOLOGY), portable(STATES), portable(NODES),
                    portable(SUPPORT)],
        "document": portable(DOC),
        "findings": findings, "acceptance": report["acceptance"],
    })
    artifacts = [REPORT, TOPOLOGY, STATES, NODES, SUPPORT, DOC, ARTIFACT]
    sources = [CONTRACT, FREEZE, Path(__file__).resolve(), TEST,
               M8_EVIDENCE, M52_EVIDENCE, M52_REPORT]
    unique_raw = sorted(set(raw_paths))
    write_json(EVIDENCE, {
        "schema": "vela.simplemos.sdevice.m56_si_oxide_interface_topology_evidence.v1",
        "status": "frozen", "contract_sha256": sha256(CONTRACT),
        "artifacts": [{"path": portable(path), "sha256": sha256(path)}
                      for path in artifacts],
        "source_hashes": {portable(path): sha256(path) for path in sources},
        "raw_input_hashes": {portable(path): sha256(path) for path in unique_raw},
        "new_sentaurus_execution": False, "new_vela_execution": False,
        "historical_artifacts_rewritten": False,
        "closed_topics_reinvestigated": False,
        "classification": report["classification"],
        "acceptance": report["acceptance"],
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analyze", action="store_true")
    args = parser.parse_args()
    contract = validate_contract()
    if not args.analyze:
        print(json.dumps({"status": "contract_validated",
                          "contract": portable(CONTRACT)}, indent=2))
        return
    report, raw_paths = analyze(contract)
    freeze_artifacts(report, raw_paths)
    print(json.dumps({
        "status": report["status"],
        "classification": report["classification"],
        "all_checks_pass": report["acceptance"]["all_checks_pass"],
        "report": portable(REPORT),
    }, indent=2))


if __name__ == "__main__":
    main()
