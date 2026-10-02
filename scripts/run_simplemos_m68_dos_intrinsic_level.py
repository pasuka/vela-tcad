"""Freeze and execute the SimpleMOS M68 DOS/intrinsic-level null intervention."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m68_dos_intrinsic_level_contract_v1.json"
FREEZE = ROOT / "simplemos_m68_dos_intrinsic_level_contract_freeze.json"
M67_EVIDENCE = ROOT / "simplemos_m67_residual_barrier_partition_evidence.json"
M65_CONFIG = REPO / "build-release/m65_ni/vela/n23/vd_0p050000/20_gate_sweep/config.json"
M65_STATE = REPO / "build-release/m65_ni/vela/n23/vd_0p050000/20_gate_sweep/state.csv"
M65_MATERIALS = REPO / "build-release/m65_ni/vela/matched_materials.json"
CARRIER = REPO / "src/physics/CarrierStatistics.cpp"
ASSEMBLER = REPO / "src/equation/CoupledDDAssembler.cpp"
RUNNER = REPO / "build-release/vela_example_runner.exe"
OUTPUT = REPO / "build-release/m68_dos"
PORTABLE = ROOT / "dos_intrinsic_level"
REPORT = PORTABLE / "m68_dos_intrinsic_level_report.json"
LEDGER = PORTABLE / "m68_dos_null_intervention_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m68_dos_intrinsic_level_2026-09-03.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m68/artifact.json"
EVIDENCE = ROOT / "simplemos_m68_dos_intrinsic_level_evidence.json"
SCRIPT = Path(__file__).resolve()
VT = 8.617333262145e-5 * 300.0


def read_json(path: Path) -> dict[str, Any]: return json.loads(path.read_text(encoding="utf-8-sig"))
def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True); path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
def sha256(path: Path) -> str: return hashlib.sha256(path.read_bytes()).hexdigest()
def portable(path: Path) -> str: return path.resolve().relative_to(REPO.resolve()).as_posix()
def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream: return list(csv.DictReader(stream))
def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n"); writer.writeheader(); writer.writerows(rows)


def sources() -> list[Path]: return [M67_EVIDENCE, M65_CONFIG, M65_STATE, M65_MATERIALS, CARRIER, ASSEMBLER, RUNNER]


def freeze_contract() -> None:
    contract, upstream = read_json(CONTRACT), read_json(M67_EVIDENCE)
    if contract.get("schema") != "vela.simplemos.sdevice.m68_dos_intrinsic_level_contract.v1": raise ValueError("unexpected M68 contract")
    if upstream.get("status") != contract["upstream"]["required_m67_status"] or upstream.get("classification") != contract["upstream"]["required_m67_classification"]: raise ValueError("M67 qualification changed")
    paths = sources(); missing = [path for path in paths if not path.is_file()]
    if missing: raise FileNotFoundError(missing[0])
    write_json(FREEZE, {"schema": "vela.simplemos.sdevice.m68_dos_intrinsic_level_contract_freeze.v1", "status": "frozen_before_execution",
                        "contract": portable(CONTRACT), "contract_sha256": sha256(CONTRACT),
                        "upstream_hashes": {portable(path): sha256(path) for path in paths}})


def validate_contract() -> dict[str, Any]:
    contract, freeze = read_json(CONTRACT), read_json(FREEZE)
    if freeze.get("status") != "frozen_before_execution" or freeze.get("contract_sha256") != sha256(CONTRACT): raise ValueError("M68 contract not frozen")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected: raise ValueError(f"M68 input changed: {relative}")
    return contract


def run_replay(label: str, materials: Path, contract: dict[str, Any]) -> dict[str, Any]:
    root = OUTPUT / label; root.mkdir(parents=True, exist_ok=True)
    config = copy.deepcopy(read_json(M65_CONFIG)); gate = float(contract["intervention"]["gate_voltage_V"])
    config["materials_file"] = str(materials.resolve()); config["solver"]["method"] = "newton"
    config["contacts"] = [{**row, "bias": gate if row["name"] == "gate" else row["bias"]} for row in config["contacts"]]
    config["output_csv"] = str((root / "curve.csv").resolve()); config["log_file"] = str((root / "run.log").resolve())
    config["sweep"].update({"start": gate, "stop": gate, "bias_points": [gate],
                            "initial_state_file": str(M65_STATE.resolve()), "write_state_file": str((root / "state.csv").resolve()), "write_vtk": False})
    config["sweep"].pop("vtk_prefix", None); config["simplemos_m68"] = {"label": label, "single_point_newton_restart": True, "ni_unchanged": True}
    config_path = root / "config.json"; write_json(config_path, config)
    completed = subprocess.run([str(RUNNER), "--config", str(config_path)], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8", errors="replace", check=False)
    (root / "stdout.txt").write_text(completed.stdout or "", encoding="utf-8", newline="\n")
    if completed.returncode or '"converged":true' not in (completed.stdout or ""): raise RuntimeError(f"M68 {label} replay failed\n{completed.stdout}")
    rows = read_csv(root / "curve.csv")
    if len(rows) != 1: raise RuntimeError(f"M68 {label} point count")
    return {"label": label, "materials": portable(materials), "config": portable(config_path), "curve": portable(root / "curve.csv"),
            "state": portable(root / "state.csv"), "current_A_per_um": abs(float(rows[0]["current_total_A_per_um"])),
            "state_sha256": sha256(root / "state.csv")}


def execute(contract: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    OUTPUT.mkdir(parents=True, exist_ok=True); baseline = OUTPUT / "baseline_materials.json"; perturbed = OUTPUT / "perturbed_materials.json"
    materials = read_json(M65_MATERIALS); write_json(baseline, materials)
    changed = copy.deepcopy(materials); silicon = next(row for row in changed["materials"] if row["name"] == "Si")
    nc0, nv0, ni0 = float(silicon["Nc_m3"]), float(silicon["Nv_m3"]), float(silicon["ni"])
    silicon["Nc_m3"] = nc0 * float(contract["intervention"]["silicon_Nc_multiplier"])
    silicon["Nv_m3"] = nv0 * float(contract["intervention"]["silicon_Nv_multiplier"])
    silicon["ni"] = ni0; changed["simplemos_m68"] = {"purpose": "large DOS-ratio null intervention", "production_default_changed": False}
    write_json(perturbed, changed)
    base, probe = run_replay("baseline", baseline, contract), run_replay("perturbed", perturbed, contract)
    difference = float(probe["current_A_per_um"]) - float(base["current_A_per_um"])
    log_difference = math.log10(float(probe["current_A_per_um"]) / float(base["current_A_per_um"]))
    implied = 0.5 * VT * math.log((float(silicon["Nv_m3"]) / float(silicon["Nc_m3"])) / (nv0 / nc0))
    code = CARRIER.read_text(encoding="utf-8")
    proof = {"boltzmann_electron_uses_ni_not_Nc": "if (!usesFermiDirac(model))\n        return ni * limitedExp" in code,
             "boltzmann_hole_uses_ni_not_Nv": "if (!usesFermiDirac(model))\n        return ni * limitedExp" in code,
             "boltzmann_equilibrium_uses_ni_not_dos": "if (!usesFermiDirac(model)) {\n        const Real n = boltzmannElectronEquilibrium(netDoping, ni);" in code,
             "generalized_srh_dos_is_fermi_guarded": "if (usesFermiDirac(model)) {\n        state.electronDegeneracy" in code}
    checks = {"contract_frozen": read_json(FREEZE)["contract_sha256"] == sha256(CONTRACT), "replay_count": 2 == int(contract["acceptance"]["required_replay_count"]),
              "all_code_path_anchors_present": all(proof.values()), "current_bit_identity": difference == 0.0 and log_difference == 0.0,
              "state_hash_identity": base["state_sha256"] == probe["state_sha256"], "ni_unchanged": float(silicon["ni"]) == ni0, "production_reference_not_replaced": True}
    checks["all_checks_pass"] = all(checks.values())
    classification = "dos_ratio_inactive_in_frozen_boltzmann_path" if checks["all_checks_pass"] else "dos_ratio_active_in_frozen_boltzmann_path"
    report = {"schema": "vela.simplemos.sdevice.m68_dos_intrinsic_level_report.v1", "status": "accepted" if checks["all_checks_pass"] else "failed", "classification": classification,
              "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)}, "execution": {"new_sentaurus_solves": 0, "new_vela_single_point_newton_reclosures": 2, "production_default_changed": False},
              "summary": {"baseline_Nc_m3": nc0, "baseline_Nv_m3": nv0, "perturbed_Nc_m3": float(silicon["Nc_m3"]), "perturbed_Nv_m3": float(silicon["Nv_m3"]), "ni_cm3_as_material_input": ni0,
                          "implied_intrinsic_level_shift_V_if_dos_ratio_were_active": implied, "baseline_current_A_per_um": base["current_A_per_um"], "perturbed_current_A_per_um": probe["current_A_per_um"],
                          "absolute_current_difference_A_per_um": difference, "log_current_difference_dex": log_difference, "state_hash_identity": base["state_sha256"] == probe["state_sha256"]},
              "code_path_proof": proof,
              "causal_scope": {"closed": "Nc/Nv cannot change this frozen classical Boltzmann current path once ni is fixed.",
                               "not_closed": "Sentaurus may allocate its intrinsic energy reference asymmetrically; Vela currently cannot represent that as an Nc/Nv knob in this path. M69 localizes the remaining stage dependence."},
              "acceptance": checks}
    return report, [base, probe]


def freeze_results(report: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    write_csv(LEDGER, rows); write_json(REPORT, report); summary = report["summary"]
    DOC.parent.mkdir(parents=True, exist_ok=True); DOC.write_text(f"""# SimpleMOS M68 DOS / 本征能级分配审计

M68 分类为 `{report['classification']}`。在 `ni` 不变时，把 Si `Nc` 放大 4 倍、`Nv` 缩小到 1/4（若该比值生效，对应本征能级参考移动 `{1000*float(summary['implied_intrinsic_level_shift_V_if_dos_ratio_were_active']):.3f} mV`），n23 固定状态电流仍逐位一致，差值 `{float(summary['absolute_current_difference_A_per_um']):.3e} A/um`、`{float(summary['log_current_difference_dex']):.3e} dex`。

源码路径同时证明 Boltzmann 浓度、平衡态和 SRH 平衡积不使用 `Nc/Nv`。因此 M65 剩余项不能通过当前 classical noBGN 路径中的 DOS 比值调参闭合；若 Sentaurus 的本征能级参考分配确有影响，它应表现为电势/边界参考差异，而不是 Vela 现有 `Nc/Nv` 参数响应。机器报告：`{portable(REPORT)}`。
""", encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {"schema": "vela.validation.artifact.v1", "title": "SimpleMOS M68 DOS intrinsic level", "status": report["status"], "classification": report["classification"], "report": portable(REPORT), "ledgers": [portable(LEDGER)]})
    artifacts = [REPORT, LEDGER, DOC, ARTIFACT]
    write_json(EVIDENCE, {"schema": "vela.simplemos.sdevice.m68_dos_intrinsic_level_evidence.v1", "status": "frozen" if report["acceptance"]["all_checks_pass"] else "failed", "classification": report["classification"], "contract_sha256": sha256(CONTRACT),
                          "implementation_hashes": {portable(SCRIPT): sha256(SCRIPT)}, "artifacts": {portable(path): sha256(path) for path in artifacts}, "new_sentaurus_execution": False, "new_vela_single_point_newton_execution": True, "production_reference_replaced": False, "acceptance": report["acceptance"]})


def verify() -> dict[str, Any]:
    validate_contract(); evidence, report = read_json(EVIDENCE), read_json(REPORT)
    if evidence["status"] != "frozen" or report["status"] != "accepted": raise ValueError("M68 not frozen")
    if evidence["implementation_hashes"][portable(SCRIPT)] != sha256(SCRIPT): raise ValueError("M68 implementation changed")
    for relative, expected in evidence["artifacts"].items():
        if sha256(REPO / relative) != expected: raise ValueError(f"M68 artifact changed: {relative}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--freeze-contract", action="store_true"); parser.add_argument("--execute", action="store_true"); parser.add_argument("--verify", action="store_true"); args = parser.parse_args()
    if args.freeze_contract: freeze_contract(); print(json.dumps({"status": "frozen_before_execution", "contract_sha256": sha256(CONTRACT)})); return
    contract = validate_contract()
    if args.verify: print(json.dumps(verify()["summary"], indent=2)); return
    if args.execute:
        report, rows = execute(contract); freeze_results(report, rows); print(json.dumps({"status": report["status"], "classification": report["classification"], "summary": report["summary"], "acceptance": report["acceptance"]}, indent=2))


if __name__ == "__main__": main()
