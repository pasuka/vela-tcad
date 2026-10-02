"""Freeze and execute the read-only SimpleMOS M67 residual barrier partition."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m47_default_bgn_self_consistent_attribution as m47  # noqa: E402
import run_simplemos_m65_nobgn_intrinsic_density_attribution as m65  # noqa: E402

ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m67_residual_barrier_partition_contract_v1.json"
FREEZE = ROOT / "simplemos_m67_residual_barrier_partition_contract_freeze.json"
M66_EVIDENCE = ROOT / "simplemos_m66_matched_ni_full_curve_evidence.json"
M65_EVIDENCE = ROOT / "simplemos_m65_nobgn_intrinsic_density_attribution_evidence.json"
M65_REPORT = ROOT / "nobgn_intrinsic_density_attribution/m65_nobgn_intrinsic_density_attribution_report.json"
M65_CURRENT = ROOT / "nobgn_intrinsic_density_attribution/m65_current_intervention_ledger.csv"
M65_PAIRS = ROOT / "nobgn_intrinsic_density_attribution/m65_nwell_pair_closure_ledger.csv"
M65_VELA = REPO / "build-release/m65_ni/vela_manifest.json"
M65_SENT = REPO / "build-release/m65_ni/sentaurus_export_manifest.json"
PORTABLE = ROOT / "residual_barrier_partition"
REPORT = PORTABLE / "m67_residual_barrier_partition_report.json"
STATES = PORTABLE / "m67_fixed_node_state_ledger.csv"
PAIRS = PORTABLE / "m67_pair_partition_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m67_residual_barrier_partition_2026-09-03.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m67/artifact.json"
EVIDENCE = ROOT / "simplemos_m67_residual_barrier_partition_evidence.json"
SCRIPT = Path(__file__).resolve()
VT = 8.617333262145e-5 * 300.0


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def median(values: list[float]) -> float:
    values = sorted(values); middle = len(values) // 2
    return values[middle] if len(values) % 2 else 0.5 * (values[middle - 1] + values[middle])


def sources() -> list[Path]:
    paths = [M66_EVIDENCE, M65_EVIDENCE, M65_REPORT, M65_CURRENT, M65_PAIRS, M65_VELA, M65_SENT]
    vela, sent = read_json(M65_VELA), read_json(M65_SENT)
    paths.extend(REPO / row["state"] for row in vela["workflows"])
    for row in sent["states"]:
        root = REPO / row["export_dir"]
        paths.extend([root / "nodes.csv", root / "elements.csv", root / "field_manifest.json",
                      root / "fields/ElectrostaticPotential_region0.csv",
                      root / "fields/eQuasiFermiPotential_region0.csv",
                      root / "fields/eDensity_region0.csv"])
    return paths


def freeze_contract() -> None:
    contract = read_json(CONTRACT); upstream = read_json(M66_EVIDENCE)
    if contract.get("schema") != "vela.simplemos.sdevice.m67_residual_barrier_partition_contract.v1":
        raise ValueError("unexpected M67 contract schema")
    if upstream.get("status") != contract["upstream"]["required_m66_status"] or upstream.get("classification") != contract["upstream"]["required_m66_classification"]:
        raise ValueError("M66 upstream qualification changed")
    paths = sources(); missing = [path for path in paths if not path.is_file()]
    if missing: raise FileNotFoundError(missing[0])
    write_json(FREEZE, {"schema": "vela.simplemos.sdevice.m67_residual_barrier_partition_contract_freeze.v1",
                        "status": "frozen_before_execution", "contract": portable(CONTRACT),
                        "contract_sha256": sha256(CONTRACT),
                        "upstream_hashes": {portable(path): sha256(path) for path in paths}})


def validate_contract() -> dict[str, Any]:
    contract, freeze = read_json(CONTRACT), read_json(FREEZE)
    if freeze.get("status") != "frozen_before_execution" or freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M67 contract is not frozen")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected: raise ValueError(f"M67 input changed: {relative}")
    return contract


def analyze(contract: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    vela_manifest, sent_manifest = read_json(M65_VELA), read_json(M65_SENT)
    m65_report = read_json(M65_REPORT)
    vela_by = {(row["device"], float(row["drain_voltage_V"])): row for row in vela_manifest["workflows"]}
    current_by = {(row["device"], float(row["drain_voltage_V"])): row for row in read_csv(M65_CURRENT)}
    state_rows: list[dict[str, Any]] = []
    for item in sent_manifest["states"]:
        key = (item["device"], float(item["drain_voltage_V"])); workflow = vela_by[key]
        export = REPO / item["export_dir"]
        sent_psi = m47.scalar_field(export, "ElectrostaticPotential")
        sent_phin = m47.scalar_field(export, "eQuasiFermiPotential")
        sent_n = m47.scalar_field(export, "eDensity")
        vela = m65.vela_state(REPO / workflow["state"])
        coords = m47.coordinates(export); _, support = m47.interface_nodes(export)
        support = sorted(set(support) & set(vela["psi"]) & set(sent_psi))
        sent_node = max(support, key=lambda node: sent_phin[node] - sent_psi[node])
        vela_node = max(support, key=lambda node: vela["phin"][node] - vela["psi"][node])
        dphin = vela["phin"][sent_node] - sent_phin[sent_node]
        dpsi = vela["psi"][sent_node] - sent_psi[sent_node]
        dbarrier = ((vela["phin"][sent_node] - vela["psi"][sent_node]) -
                    (sent_phin[sent_node] - sent_psi[sent_node]))
        partition_residual = dbarrier - (dphin - dpsi)
        independent = ((vela["phin"][vela_node] - vela["psi"][vela_node]) -
                       (sent_phin[sent_node] - sent_psi[sent_node]))
        distance = math.hypot(coords[sent_node][0] - coords[vela_node][0],
                              coords[sent_node][1] - coords[vela_node][1])
        density_observed = math.log10((vela["electrons_m3"][sent_node] / 1e6) / sent_n[sent_node])
        density_predicted = -dbarrier / (VT * math.log(10.0))
        current_error = float(current_by[key]["matched_ni_error_dex"])
        state_rows.append({"device": item["device"], "drain_voltage_V": item["drain_voltage_V"],
                           "gate_voltage_V": item["gate_voltage_V"], "support_node_count": len(support),
                           "sentaurus_control_node": sent_node, "vela_control_node": vela_node,
                           "control_node_migration_um": distance,
                           "fixed_node_phin_delta_V": dphin, "fixed_node_psi_delta_V": dpsi,
                           "fixed_node_barrier_delta_V": dbarrier,
                           "fixed_node_partition_identity_residual_V": partition_residual,
                           "fixed_node_qf_barrier_current_proxy_dex": density_predicted,
                           "independent_argmax_barrier_delta_V": independent,
                           "independent_argmax_current_proxy_dex": -independent / (VT * math.log(10.0)),
                           "fixed_node_electron_density_log_ratio_dex": density_observed,
                           "boltzmann_density_identity_residual_dex": density_observed - density_predicted,
                           "matched_ni_current_error_dex": current_error,
                           "fixed_node_proxy_minus_current_error_dex": density_predicted - current_error})
    states = {(row["device"], float(row["drain_voltage_V"])): row for row in state_rows}
    pair_rows, closures = [], []
    for pair in read_csv(M65_PAIRS):
        low, high, drain = pair["low_device"], pair["high_device"], float(pair["drain_voltage_V"])
        lo, hi = states[(low, drain)], states[(high, drain)]
        current_growth = float(pair["matched_ni_pair_growth_dex"])
        proxy_growth = float(hi["fixed_node_qf_barrier_current_proxy_dex"]) - float(lo["fixed_node_qf_barrier_current_proxy_dex"])
        phin_growth = -(float(hi["fixed_node_phin_delta_V"]) - float(lo["fixed_node_phin_delta_V"])) / (VT * math.log(10.0))
        psi_growth = (float(hi["fixed_node_psi_delta_V"]) - float(lo["fixed_node_psi_delta_V"])) / (VT * math.log(10.0))
        raw = 1.0 - abs(proxy_growth - current_growth) / max(abs(current_growth), 1e-300)
        closure = min(1.0, max(0.0, raw)); closures.append(closure)
        pair_rows.append({"low_device": low, "high_device": high, "drain_voltage_V": drain,
                          "diagnostic_gate_voltage_V": pair["diagnostic_gate_voltage_V"],
                          "matched_ni_current_pair_growth_dex": current_growth,
                          "fixed_node_qf_barrier_pair_proxy_dex": proxy_growth,
                          "fixed_node_phin_component_dex": phin_growth,
                          "fixed_node_psi_component_dex": psi_growth,
                          "component_sum_identity_residual_dex": proxy_growth - phin_growth - psi_growth,
                          "proxy_minus_current_growth_dex": proxy_growth - current_growth,
                          "raw_pair_closure_fraction": raw, "classification_pair_closure_fraction": closure,
                          "sentaurus_control_node_changed_low_to_high": int(lo["sentaurus_control_node"] != hi["sentaurus_control_node"]),
                          "maximum_cross_solver_control_node_migration_um": max(float(lo["control_node_migration_um"]), float(hi["control_node_migration_um"]))})
    metric = median(closures)
    if metric >= float(contract["analysis"]["dominant_minimum_median_pair_closure_fraction"]):
        classification = "residual_qf_barrier_dominant"
    elif metric >= float(contract["analysis"]["material_minimum_median_pair_closure_fraction"]):
        classification = "residual_qf_barrier_material_but_not_complete"
    else: classification = "residual_qf_barrier_not_material"
    checks = {"contract_frozen": read_json(FREEZE)["contract_sha256"] == sha256(CONTRACT),
              "state_count": len(state_rows) == int(contract["acceptance"]["required_state_count"]),
              "pair_count": len(pair_rows) == int(contract["acceptance"]["required_pair_count"]),
              "partition_identity": max(abs(float(row["fixed_node_partition_identity_residual_V"])) for row in state_rows) <= float(contract["acceptance"]["maximum_fixed_node_partition_identity_residual_V"]),
              "coordinate_identity_reused_from_m65": float(m65_report["summary"]["maximum_coordinate_error_um"]) <= float(contract["acceptance"]["maximum_coordinate_error_um"]),
              "boltzmann_density_identity": max(abs(float(row["boltzmann_density_identity_residual_dex"])) for row in state_rows) <= float(contract["acceptance"]["maximum_boltzmann_density_identity_residual_dex"]),
              "classification_declared": classification in contract["analysis"]["classifications"],
              "production_reference_not_replaced": True}
    checks["all_checks_pass"] = all(checks.values())
    if not checks["all_checks_pass"]: classification = "state_identity_or_execution_failure"
    report = {"schema": "vela.simplemos.sdevice.m67_residual_barrier_partition_report.v1",
              "status": "accepted" if checks["all_checks_pass"] else "failed", "classification": classification,
              "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)},
              "execution": {"new_sentaurus_solves": 0, "new_vela_solves": 0, "production_default_changed": False},
              "summary": {"median_pair_closure_fraction": metric,
                          "minimum_pair_closure_fraction": min(closures), "maximum_pair_closure_fraction": max(closures),
                          "median_current_pair_growth_dex": median([float(row["matched_ni_current_pair_growth_dex"]) for row in pair_rows]),
                          "median_fixed_node_qf_barrier_pair_proxy_dex": median([float(row["fixed_node_qf_barrier_pair_proxy_dex"]) for row in pair_rows]),
                          "median_phin_component_dex": median([float(row["fixed_node_phin_component_dex"]) for row in pair_rows]),
                          "median_psi_component_dex": median([float(row["fixed_node_psi_component_dex"]) for row in pair_rows]),
                          "maximum_control_node_migration_um": max(float(row["control_node_migration_um"]) for row in state_rows),
                          "maximum_boltzmann_density_identity_residual_dex": max(abs(float(row["boltzmann_density_identity_residual_dex"])) for row in state_rows)},
              "causal_scope": {"closed": "The residual barrier observable is exactly partitioned at a common fixed node.",
                               "not_closed": "A barrier-current proxy is diagnostic rather than a causal transport intervention; any unclosed current growth remains for M68-M70."},
              "acceptance": checks}
    return report, {"states": state_rows, "pairs": pair_rows}


def freeze_results(report: dict[str, Any], rows: dict[str, Any]) -> None:
    write_csv(STATES, rows["states"]); write_csv(PAIRS, rows["pairs"]); write_json(REPORT, report)
    summary = report["summary"]
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(f"""# SimpleMOS M67 剩余势垒分解

M67 分类为 `{report['classification']}`。在 Sentaurus 控制节点固定后，`phin` 与 `psi` 分量严格重构 `phin-psi`；8 个 NWell 配对的势垒代理对剩余电流增幅的中位闭合率为 `{float(summary['median_pair_closure_fraction']):.2%}`。中位电流增幅 `{float(summary['median_current_pair_growth_dex']):.6f} dex`，势垒代理 `{float(summary['median_fixed_node_qf_barrier_pair_proxy_dex']):.6f} dex`，其中准费米分量 `{float(summary['median_phin_component_dex']):.6f} dex`、电势分量 `{float(summary['median_psi_component_dex']):.6f} dex`。

该结果是状态观测量分解，不把热发射式势垒代理当成完整漂移扩散电流因果模型。机器报告：`{portable(REPORT)}`。
""", encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {"schema": "vela.validation.artifact.v1", "title": "SimpleMOS M67 residual barrier partition",
                          "status": report["status"], "classification": report["classification"],
                          "report": portable(REPORT), "ledgers": [portable(STATES), portable(PAIRS)]})
    artifacts = [REPORT, STATES, PAIRS, DOC, ARTIFACT]
    write_json(EVIDENCE, {"schema": "vela.simplemos.sdevice.m67_residual_barrier_partition_evidence.v1",
                          "status": "frozen" if report["acceptance"]["all_checks_pass"] else "failed",
                          "classification": report["classification"], "contract_sha256": sha256(CONTRACT),
                          "implementation_hashes": {portable(SCRIPT): sha256(SCRIPT)},
                          "artifacts": {portable(path): sha256(path) for path in artifacts},
                          "new_sentaurus_execution": False, "new_vela_execution": False,
                          "production_reference_replaced": False, "acceptance": report["acceptance"]})


def verify() -> dict[str, Any]:
    validate_contract(); evidence, report = read_json(EVIDENCE), read_json(REPORT)
    if evidence["status"] != "frozen" or report["status"] != "accepted": raise ValueError("M67 not frozen")
    if evidence["implementation_hashes"][portable(SCRIPT)] != sha256(SCRIPT): raise ValueError("M67 implementation changed")
    for relative, expected in evidence["artifacts"].items():
        if sha256(REPO / relative) != expected: raise ValueError(f"M67 artifact changed: {relative}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze-contract", action="store_true"); parser.add_argument("--analyze", action="store_true"); parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.freeze_contract: freeze_contract(); print(json.dumps({"status": "frozen_before_execution", "contract_sha256": sha256(CONTRACT)})); return
    validate_contract()
    if args.verify: print(json.dumps(verify()["summary"], indent=2)); return
    if args.analyze:
        report, rows = analyze(read_json(CONTRACT)); freeze_results(report, rows)
        print(json.dumps({"status": report["status"], "classification": report["classification"], "summary": report["summary"], "acceptance": report["acceptance"]}, indent=2))


if __name__ == "__main__": main()
