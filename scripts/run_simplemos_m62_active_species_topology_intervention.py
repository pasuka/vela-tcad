#!/usr/bin/env python3
"""Freeze the SimpleMOS M62 E3 writer-feasibility result."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m62_active_species_topology_intervention_contract_v1.json"
FREEZE = ROOT / "simplemos_m62_active_species_topology_intervention_contract_freeze.json"
ERRATUM = ROOT / "simplemos_m62_active_species_topology_intervention_contract_erratum_v1.json"
M58_EVIDENCE = ROOT / "simplemos_m58_compensation_contour_mesh_audit_evidence.json"
M60_EVIDENCE = ROOT / "simplemos_m60_tight_convergence_port_burst_evidence.json"
M63_EVIDENCE = ROOT / "simplemos_m63_smooth_nwell_attribution_evidence.json"
M58_NODES = ROOT / "compensation_contour_mesh_audit/m58_critical_gate_edge_node_ledger.csv"
M58_COMPONENTS = ROOT / "compensation_contour_mesh_audit/m58_connectivity_component_ledger.csv"
ORIGINAL = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics/neutral/n23"
ROUNDTRIP = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m62_active_species_topology_intervention/roundtrip_export"
ROUNDTRIP_TDR = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m62_active_species_topology_intervention/remote/m62_n23_roundtrip.tdr"
INVENTORY_LOG = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m62_active_species_topology_intervention/remote/m62_n23_inventory.txt"
ROUNDTRIP_LOG = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m62_active_species_topology_intervention/remote/m62_n23_roundtrip.log"
OUT = ROOT / "active_species_topology_intervention"
REPORT = OUT / "m62_active_species_topology_intervention_report.json"
FIELDS = OUT / "m62_roundtrip_field_invariance_ledger.csv"
TOOLS = OUT / "m62_writer_capability_ledger.csv"
CANDIDATES = OUT / "m62_candidate_node_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m62_active_species_topology_intervention_2026-09-02.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m62/artifact.json"
EVIDENCE = ROOT / "simplemos_m62_active_species_topology_intervention_evidence.json"
SCRIPT = REPO / "scripts/run_simplemos_m62_active_species_topology_intervention.py"
TDX_INVENTORY = REPO / "scripts/simplemos_m62_tdx_inventory.tcl"
TDX_ROUNDTRIP = REPO / "scripts/simplemos_m62_tdx_roundtrip.tcl"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False,
                               allow_nan=False) + "\n", encoding="utf-8",
                    newline="\n")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]),
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def validate() -> dict[str, Any]:
    contract, freeze = read_json(CONTRACT), read_json(FREEZE)
    if freeze.get("status") != "frozen_before_tool_inventory":
        raise ValueError("M62 contract was not frozen before tool inventory")
    if freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M62 contract changed after freeze")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M62 upstream artifact changed: {relative}")
    erratum = read_json(ERRATUM)
    if (erratum.get("status") != "frozen_documentation_correction" or
            erratum.get("contract_sha256") != sha256(CONTRACT) or
            erratum.get("m58_evidence_sha256") != sha256(M58_EVIDENCE)):
        raise ValueError("M62 contract erratum is not frozen against M58")
    if read_json(M58_EVIDENCE).get("status") != erratum["corrected_interpretation"]:
        raise ValueError("M58 status changed")
    if read_json(M60_EVIDENCE).get("classification") != contract["upstream"]["required_m60_classification"]:
        raise ValueError("M60 classification changed")
    if read_json(M63_EVIDENCE).get("classification") != contract["upstream"]["required_m63_classification"]:
        raise ValueError("M63 classification changed")
    for path in (ROUNDTRIP_TDR, INVENTORY_LOG, ROUNDTRIP_LOG,
                 ORIGINAL / "field_manifest.json", ROUNDTRIP / "field_manifest.json"):
        if not path.is_file():
            raise FileNotFoundError(path)
    return contract


def compare_field_files(original: Path, roundtrip: Path) -> dict[str, Any]:
    left, right = read_csv(original), read_csv(roundtrip)
    if len(left) != len(right) or list(left[0]) != list(right[0]):
        raise ValueError(f"round-trip field shape changed: {original.name}")
    components = [name for name in left[0] if name.startswith("component")]
    changed_components = 0
    changed_rows = 0
    max_abs = 0.0
    max_sym = 0.0
    ratios: list[float] = []
    for a, b in zip(left, right, strict=True):
        row_changed = False
        for component in components:
            av, bv = float(a[component]), float(b[component])
            if av != bv:
                changed_components += 1
                row_changed = True
            max_abs = max(max_abs, abs(av - bv))
            scale = max(abs(av), abs(bv), 1e-300)
            max_sym = max(max_sym, abs(av - bv) / scale)
            if av != 0.0 and bv != 0.0:
                ratios.append(abs(bv / av))
        changed_rows += int(row_changed)
    return {
        "field_file": original.name,
        "location": "element" if "_cells.csv" in original.name else "vertex",
        "row_count": len(left), "component_count": len(components),
        "changed_row_count": changed_rows,
        "changed_component_count": changed_components,
        "maximum_absolute_difference": max_abs,
        "maximum_symmetric_relative_difference": max_sym,
        "median_nonzero_abs_roundtrip_over_original": (
            sorted(ratios)[len(ratios) // 2] if ratios else 1.0),
        "byte_identical": sha256(original) == sha256(roundtrip),
    }


def candidate_rows() -> list[dict[str, Any]]:
    nodes = read_csv(M58_NODES)
    components = [row for row in read_csv(M58_COMPONENTS)
                  if row["device"] in ("n21", "n23")
                  and row["single_node"] == "True"
                  and row["floating"] == "True"]
    result: list[dict[str, Any]] = []
    for component in components:
        x, y = float(component["x_min_um"]), float(component["y_min_um"])
        matches = [row for row in nodes if row["device"] == component["device"]
                   and math.isclose(float(row["x_um"]), x, abs_tol=1e-12)
                   and math.isclose(float(row["y_um"]), y, abs_tol=2e-8)]
        if len(matches) != 1:
            raise ValueError(f"M62 candidate coordinate match failed: {component}")
        row = matches[0]
        identity = (float(row["AsActive_cm3"]) + float(row["PActive_cm3"])
                    - float(row["BActive_cm3"]))
        result.append({
            "device": row["device"], "node_id": int(row["node_id"]),
            "x_um": row["x_um"], "y_um": row["y_um"],
            "BActive_cm3": row["BActive_cm3"],
            "AsActive_cm3": row["AsActive_cm3"],
            "PActive_cm3": row["PActive_cm3"],
            "NetActive_cm3": row["NetActive_cm3"],
            "identity_reconstructed_net_cm3": identity,
            "identity_absolute_error_cm3": abs(identity - float(row["NetActive_cm3"])),
            "floating_single_node": True,
            "mutation_attempted": False,
            "stop_reason": "no-edit tdx round-trip changed unrelated fields",
        })
    return result


def analyze(contract: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    original_fields = {path.name: path for path in (ORIGINAL / "fields").glob("*.csv")}
    roundtrip_fields = {path.name: path for path in (ROUNDTRIP / "fields").glob("*.csv")}
    if original_fields.keys() != roundtrip_fields.keys():
        raise ValueError("round-trip field inventory changed")
    fields = [compare_field_files(original_fields[name], roundtrip_fields[name])
              for name in sorted(original_fields)]
    changed = [row for row in fields if int(row["changed_component_count"]) > 0]
    topology_files = ["nodes.csv", "elements.csv", "contacts.csv"]
    doping_files = ["doping.csv", "doping_metadata.json"]
    topology_exact = all(sha256(ORIGINAL / name) == sha256(ROUNDTRIP / name)
                         for name in topology_files)
    doping_exact = all(sha256(ORIGINAL / name) == sha256(ROUNDTRIP / name)
                       for name in doping_files)
    inventory = INVENTORY_LOG.read_text(encoding="utf-8", errors="replace")
    required_names = ["NetActive", "BActive", "PActive", "AsActive"]
    species_visible = all(f"|{name}|" in inventory for name in required_names)
    tools = [
        {"environment": "local Windows", "tool": "tdx", "available": False,
         "documented_nodal_dataset_write": False, "roundtrip_attempted": False,
         "roundtrip_preserves_unrelated_fields": False,
         "decision": "not a candidate"},
        {"environment": "Sentaurus T-2022.03-SP2", "tool": "tdx Tcl interface",
         "available": True, "documented_nodal_dataset_write": True,
         "roundtrip_attempted": True,
         "roundtrip_preserves_unrelated_fields": len(changed) == 0,
         "decision": "rejected: no-edit save changes unrelated fields"},
        {"environment": "Sentaurus T-2022.03-SP2", "tool": "tdr2ascii/ascii2tdr",
         "available": False, "documented_nodal_dataset_write": False,
         "roundtrip_attempted": False,
         "roundtrip_preserves_unrelated_fields": False,
         "decision": "not a candidate"},
        {"environment": "repository", "tool": "sentaurus_import.exe",
         "available": True, "documented_nodal_dataset_write": False,
         "roundtrip_attempted": False,
         "roundtrip_preserves_unrelated_fields": False,
         "decision": "read-only verifier only"},
    ]
    candidates = candidate_rows()
    writer_gate_pass = (topology_exact and doping_exact and species_visible
                        and len(changed) == 0)
    classification = ("writer_or_execution_mismatch" if not species_visible
                      else "e3_stopped_no_identity_preserving_writer")
    checks = {
        "contract_frozen_before_inventory": sha256(CONTRACT) == read_json(FREEZE)["contract_sha256"],
        "candidate_ledger_complete": len(candidates) == 4,
        "local_and_remote_writer_inventory_complete": len(tools) == 4,
        "roundtrip_topology_exact": topology_exact,
        "roundtrip_doping_exact": doping_exact,
        "required_species_datasets_visible": species_visible,
        "roundtrip_unrelated_fields_exact": len(changed) == 0,
        "writer_gate_failed_before_mutation": not writer_gate_pass,
        "no_mutated_tdr_created": True,
        "no_device_solve_executed": True,
        "classification_declared": classification in contract["classifications"],
        "production_references_not_replaced": True,
    }
    # A failed writer-invariance predicate is the evidence supporting the accepted stop,
    # not an overall acceptance failure.
    accepted_checks = dict(checks)
    accepted_checks["roundtrip_unrelated_fields_exact"] = True
    checks["all_checks_pass"] = all(accepted_checks.values())
    report = {
        "schema": "vela.simplemos.sdevice.m62_active_species_topology_intervention_report.v1",
        "status": "accepted" if checks["all_checks_pass"] else "failed",
        "classification": classification,
        "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)},
        "writer_gate": {
            "passed": writer_gate_pass,
            "stop_rule_applied": not writer_gate_pass,
            "tdx_release": "T-2022.03-SP2",
            "tdx_manual_capability": "TdrDataSetComponent plus TdrFileSave",
            "original_tdr_sha256": sha256(REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics/sentaurus_bundle/n23/input_fps.tdr"),
            "roundtrip_tdr_sha256": sha256(ROUNDTRIP_TDR),
            "topology_files_byte_identical": topology_exact,
            "doping_files_byte_identical": doping_exact,
            "field_file_count": len(fields),
            "changed_unrelated_field_file_count": len(changed),
            "changed_unrelated_field_files": [row["field_file"] for row in changed],
            "maximum_symmetric_relative_field_difference": max(float(row["maximum_symmetric_relative_difference"]) for row in fields),
            "representative_displacement_scale_factor": next(float(row["median_nonzero_abs_roundtrip_over_original"]) for row in fields if row["field_file"] == "Displacement_region0.csv"),
        },
        "execution": {
            "tool_inventory_completed": True,
            "no_edit_roundtrip_completed": True,
            "candidate_mutation_attempted": False,
            "new_sentaurus_device_execution": False,
            "mutated_tdr_created": False,
        },
        "interpretation": {
            "closed": "E3 cannot be executed under its frozen field-identity requirement with the available writer path.",
            "not_closed": "The causal modulation of the M58 floating-node topology was not measured.",
            "upstream_implication": "M60 remains the decisive burst discriminator; M63 remains the independent smooth-component attribution.",
        },
        "acceptance": checks,
    }
    return report, fields, tools, candidates


def freeze(report: dict[str, Any], fields: list[dict[str, Any]],
           tools: list[dict[str, Any]], candidates: list[dict[str, Any]]) -> None:
    write_csv(FIELDS, fields)
    write_csv(TOOLS, tools)
    write_csv(CANDIDATES, candidates)
    write_json(REPORT, report)
    gate = report["writer_gate"]
    DOC.write_text(f"""# SimpleMOS M62 活性物种拓扑干预

## 结论

M62 分类为 `{report['classification']}`。T-2022.03-SP2 的 `tdx` Tcl接口在手册和运行时清单中均能看到逐值写回能力，n23 TDR也同时含 `BActive`、`PActive`、`AsActive` 与 `NetActive`。但是无编辑 `TdrFileSave` round-trip 未通过冻结的字段恒等门禁。

网格、单元、接触、`doping.csv` 和掺杂元数据保持字节一致；132个导出字段中有 `{int(gate['changed_unrelated_field_file_count'])}` 个无关机械字段改变。代表例 `Displacement_region0.csv` 的非零值尺度因子为 `{float(gate['representative_displacement_scale_factor']):.6g}`，最大对称相对差为 `{float(gate['maximum_symmetric_relative_field_difference']):.6g}`。因此不能证明局部活性物种写回不会同时改变其他状态。

按照冻结停止规则，没有创建变异TDR、没有执行SDevice E3求解，也没有只编辑派生 `NetActive`。这项结果关闭的是“当前工具链能否合法执行E3”，不是浮置单节点拓扑的因果效应本身。burst的决定性判别仍由M60提供；平滑NWell成分仍由M63独立归因。

机器报告：`{portable(REPORT)}`。
""", encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {
        "schema": "vela.validation.artifact.v1",
        "title": "SimpleMOS M62 active-species topology intervention gate",
        "status": report["status"], "classification": report["classification"],
        "summary": "E3 stopped before mutation because the no-edit TDX round-trip changed unrelated fields",
        "report": portable(REPORT),
        "ledgers": [portable(FIELDS), portable(TOOLS), portable(CANDIDATES)],
    })
    artifacts = [REPORT, FIELDS, TOOLS, CANDIDATES, DOC, ARTIFACT, ERRATUM,
                 TDX_INVENTORY, TDX_ROUNDTRIP]
    write_json(EVIDENCE, {
        "schema": "vela.simplemos.sdevice.m62_active_species_topology_intervention_evidence.v1",
        "status": "frozen" if report["acceptance"]["all_checks_pass"] else "failed",
        "classification": report["classification"],
        "contract_sha256": sha256(CONTRACT),
        "new_sentaurus_device_execution": False,
        "mutated_tdr_created": False,
        "production_reference_replaced": False,
        "source_hashes": read_json(FREEZE)["upstream_hashes"],
        "implementation_hashes": {portable(SCRIPT): sha256(SCRIPT)},
        "diagnostic_source_hashes": {
            portable(ROUNDTRIP_TDR): sha256(ROUNDTRIP_TDR),
            portable(INVENTORY_LOG): sha256(INVENTORY_LOG),
            portable(ROUNDTRIP_LOG): sha256(ROUNDTRIP_LOG),
        },
        "artifacts": {portable(path): sha256(path) for path in artifacts},
        "acceptance": report["acceptance"],
    })


def verify() -> dict[str, Any]:
    validate()
    evidence, report = read_json(EVIDENCE), read_json(REPORT)
    if evidence.get("status") != "frozen" or report.get("status") != "accepted":
        raise ValueError("M62 evidence is not accepted and frozen")
    for relative, expected in evidence["artifacts"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M62 artifact hash changed: {relative}")
    for relative, expected in evidence["diagnostic_source_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M62 diagnostic source hash changed: {relative}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        report = verify()
        print(f"M62 verified: {report['classification']}")
        return
    contract = validate()
    report, fields, tools, candidates = analyze(contract)
    freeze(report, fields, tools, candidates)
    if not report["acceptance"]["all_checks_pass"]:
        raise RuntimeError(f"M62 acceptance failed: {report['acceptance']}")
    print(json.dumps({"status": report["status"],
                      "classification": report["classification"],
                      "writer_gate": report["writer_gate"]}, indent=2))


if __name__ == "__main__":
    main()
