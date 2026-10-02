"""Freeze, run, and analyze SimpleMOS M69 stage residual localization."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
from typing import Any, Sequence

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m47_default_bgn_self_consistent_attribution as m47  # noqa: E402
import run_simplemos_m60_tight_convergence_port_burst as m60  # noqa: E402
import run_simplemos_m65_nobgn_intrinsic_density_attribution as m65  # noqa: E402

ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m69_stage_residual_localization_contract_v2.json"
FREEZE = ROOT / "simplemos_m69_stage_residual_localization_contract_freeze_v2.json"
M68_EVIDENCE = ROOT / "simplemos_m68_dos_intrinsic_level_evidence.json"
M69_V1_EVIDENCE = ROOT / "simplemos_m69_stage_residual_localization_evidence_v1_failed.json"
M65_EVIDENCE = ROOT / "simplemos_m65_nobgn_intrinsic_density_attribution_evidence.json"
M65_VELA = REPO / "build-release/m65_ni/vela_manifest.json"
M65_SENT = REPO / "build-release/m65_ni/sentaurus_export_manifest.json"
M65_PAIRS = ROOT / "nobgn_intrinsic_density_attribution/m65_nwell_pair_closure_ledger.csv"
M8 = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics"
OUTPUT = REPO / "build-release/m69_stage_v2"
PORTABLE = ROOT / "stage_residual_localization"
REPORT = PORTABLE / "m69_stage_residual_localization_report.json"
STATES = PORTABLE / "m69_stage_state_ledger.csv"
PAIRS = PORTABLE / "m69_pair_stage_ledger.csv"
CASES = PORTABLE / "m69_execution_case_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m69_stage_residual_localization_2026-09-03.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m69/artifact.json"
EVIDENCE = ROOT / "simplemos_m69_stage_residual_localization_evidence.json"
SCRIPT = Path(__file__).resolve(); IMPORTER = REPO / "build-release/sentaurus_import.exe"
SSH_CONFIG = Path(os.environ.get("USERPROFILE", str(Path.home()))) / ".ssh/config"
REMOTE_ROOT = "/tmp/vela_simplemos_m69_stage_20260903b"; ARCHIVE = "m69_stage_v2_results.tgz"
STAGES = ("equilibrium", "drain", "gate"); VT = 8.617333262145e-5 * 300.0


def read_json(path: Path) -> dict[str, Any]: return json.loads(path.read_text(encoding="utf-8-sig"))
def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True); path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8", newline="\n")
def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream: return list(csv.DictReader(stream))
def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
def sha256(path: Path) -> str: return hashlib.sha256(path.read_bytes()).hexdigest()
def portable(path: Path) -> str: return path.resolve().relative_to(REPO.resolve()).as_posix()
def median(values: list[float]) -> float:
    values = sorted(values); middle = len(values)//2; return values[middle] if len(values)%2 else .5*(values[middle-1]+values[middle])
def voltage_tag(value: float) -> str: return f"{value:.6f}".replace("-", "m").replace(".", "p")
def run(argv: Sequence[str], capture: bool = False) -> str:
    done = subprocess.run(list(argv), cwd=REPO, text=True, stdout=subprocess.PIPE if capture else None, stderr=subprocess.STDOUT if capture else None, encoding="utf-8", errors="replace", check=False)
    if done.returncode: raise RuntimeError(f"command failed {list(argv)}\n{done.stdout or ''}")
    return done.stdout or ""


def source_paths() -> list[Path]:
    paths = [M68_EVIDENCE, M69_V1_EVIDENCE, M65_EVIDENCE, M65_VELA, M65_SENT, M65_PAIRS, IMPORTER]
    for row in read_json(M65_VELA)["workflows"]:
        root = (REPO / row["state"]).parents[1]; paths.extend([root / "00_equilibrium/state.csv", root / "10_drain_ramp/state.csv", REPO / row["state"]])
    for row in read_json(M65_SENT)["states"]: paths.append(M8 / "sentaurus_bundle" / row["device"] / "input_fps.tdr")
    return paths


def freeze_contract() -> None:
    contract, upstream = read_json(CONTRACT), read_json(M68_EVIDENCE)
    if contract.get("schema") != "vela.simplemos.sdevice.m69_stage_residual_localization_contract.v2": raise ValueError("unexpected M69 contract")
    if upstream.get("status") != contract["upstream"]["required_m68_status"] or upstream.get("classification") != contract["upstream"]["required_m68_classification"]: raise ValueError("M68 qualification changed")
    if read_json(M69_V1_EVIDENCE).get("status") != contract["upstream"]["required_m69_v1_status"]: raise ValueError("M69 v1 failure evidence changed")
    paths = source_paths(); missing = [path for path in paths if not path.is_file()]
    if missing: raise FileNotFoundError(missing[0])
    write_json(FREEZE, {"schema": "vela.simplemos.sdevice.m69_stage_residual_localization_contract_freeze.v2", "status": "frozen_before_execution", "contract": portable(CONTRACT), "contract_sha256": sha256(CONTRACT), "upstream_hashes": {portable(path): sha256(path) for path in paths}})


def validate_contract() -> dict[str, Any]:
    contract, freeze = read_json(CONTRACT), read_json(FREEZE)
    if freeze.get("status") != "frozen_before_execution" or freeze.get("contract_sha256") != sha256(CONTRACT): raise ValueError("M69 contract not frozen")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected: raise ValueError(f"M69 input changed: {relative}")
    return contract


def deck(case: str, drain: float, gate: float) -> str:
    text = m65.sentaurus_deck(case, drain, gate)
    anchor = "  Coupled { Poisson Electron Hole }\n  Quasistationary("
    if text.count(anchor) != 1: raise RuntimeError("M69 equilibrium plot anchor changed")
    text = text.replace(anchor, f'  Coupled {{ Poisson Electron Hole }}\n  Plot(FilePrefix="{case}_equilibrium")\n  Quasistationary(', 1)
    anchor = "  ) { Coupled { Poisson Electron Hole } }\n  NewCurrentPrefix=\"IdVg_\""
    if text.count(anchor) != 1: raise RuntimeError("M69 drain plot anchor changed")
    text = text.replace(anchor, f'  ) {{ Coupled {{ Poisson Electron Hole }} }}\n  Plot(FilePrefix="{case}_drain")\n  NewCurrentPrefix="IdVg_"', 1)
    old = f'Plot(FilePrefix="{case}_state" NoOverWrite Time=(1))'
    if text.count(old) != 1: raise RuntimeError("M69 gate plot anchor changed")
    return text.replace(old, f'Plot(FilePrefix="{case}_gate" NoOverWrite Time=(1))')


def prepare(force: bool) -> dict[str, Any]:
    bundle = OUTPUT / "sentaurus_bundle"
    if force and bundle.exists(): shutil.rmtree(bundle)
    items = []
    for source in read_json(M65_SENT)["states"]:
        device, drain, gate = source["device"], float(source["drain_voltage_V"]), float(source["gate_voltage_V"])
        root = bundle / device / f"vd_{voltage_tag(drain)}"; root.mkdir(parents=True, exist_ok=True)
        tdr = M8 / "sentaurus_bundle" / device / "input_fps.tdr"; shutil.copy2(tdr, root / "input_fps.tdr")
        case = f"m69_{device}_vd_{voltage_tag(drain)}_vg_{voltage_tag(gate)}"; command = root / f"{case}_des.cmd"; command.write_text(deck(case, drain, gate), encoding="utf-8", newline="\n")
        items.append({"device": device, "drain_voltage_V": drain, "gate_voltage_V": gate, "case": case, "deck": portable(command), "deck_sha256": sha256(command)})
    manifest = {"schema": "vela.simplemos.sdevice.m69_sentaurus_manifest.v1", "status": "prepared", "contract_sha256": sha256(CONTRACT), "cases": items}; write_json(OUTPUT / "sentaurus_manifest.json", manifest); return manifest


def run_sentaurus(manifest: dict[str, Any], ssh_target: str, ssh_bin: str, scp_bin: str, ssh_config: str, jobs: int) -> str:
    ssh, scp = [ssh_bin, "-F", ssh_config], [scp_bin, "-F", ssh_config]
    banner = run([*ssh, ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'"], True).strip()
    if "T-2022.03-SP2" not in banner: raise RuntimeError(f"unexpected Sentaurus release\n{banner}")
    run([*ssh, ssh_target, f"set -eu; test ! -e {REMOTE_ROOT}; mkdir -p {REMOTE_ROOT}"])
    run([*scp, "-r", str(OUTPUT / "sentaurus_bundle"), f"{ssh_target}:{REMOTE_ROOT}/"])
    def execute(item: dict[str, Any]) -> None:
        relative = (REPO / item["deck"]).parent.relative_to(OUTPUT / "sentaurus_bundle").as_posix(); root = f"{REMOTE_ROOT}/sentaurus_bundle/{relative}"; case = item["case"]
        run([*ssh, ssh_target, f"set -eu; cd {root}; sdevice {case}_des.cmd > {case}.console.log 2>&1"])
    with ThreadPoolExecutor(max_workers=min(max(1, jobs), 8)) as pool: list(pool.map(execute, manifest["cases"]))
    run([*ssh, ssh_target, f"cd {REMOTE_ROOT}; tar -czf {ARCHIVE} sentaurus_bundle"])
    raw = OUTPUT / "sentaurus_raw"; raw.mkdir(parents=True, exist_ok=True); archive = raw / ARCHIVE; run([*scp, f"{ssh_target}:{REMOTE_ROOT}/{ARCHIVE}", str(archive)])
    with tarfile.open(archive, "r:gz") as stream: stream.extractall(raw, filter="data")
    (OUTPUT / "sentaurus_banner.txt").write_text(banner + "\n", encoding="utf-8", newline="\n"); return banner


def export(manifest: dict[str, Any]) -> dict[str, Any]:
    raw = OUTPUT / "sentaurus_raw/sentaurus_bundle"; states = []; case_rows = []
    for item in manifest["cases"]:
        relative = (REPO / item["deck"]).parent.relative_to(OUTPUT / "sentaurus_bundle"); root = raw / relative; case = item["case"]
        stage_files = {}
        for stage in STAGES:
            candidates = sorted(root.glob(f"{case}_{stage}*.tdr"))
            if len(candidates) != 1: raise RuntimeError(f"M69 {case} {stage} TDR count {len(candidates)}")
            stage_files[stage] = candidates[0]
        plots = sorted(root.glob(f"IdVg_{case}*des.plt"))
        if len(plots) != 1: raise RuntimeError(f"M69 {case} current count {len(plots)}")
        parsed = m60.parse_current(plots[0], [float(item["gate_voltage_V"])], 1e-10); current = abs(float(parsed[0]["drain_total_A_per_um"]))
        case_rows.append({**item, "current": portable(plots[0]), "final_drain_current_A_per_um": current})
        for stage, tdr in stage_files.items():
            state_id = f"{item['device']}_vd_{voltage_tag(float(item['drain_voltage_V']))}_{stage}"; target = OUTPUT / "sentaurus_exports" / state_id
            run([str(IMPORTER), "--tdr", str(tdr), "--export-dir", str(target)], True)
            states.append({**item, "stage": stage, "state_id": state_id, "tdr": portable(tdr), "export_dir": portable(target)})
    result = {"schema": "vela.simplemos.sdevice.m69_sentaurus_exports.v1", "status": "complete", "states": states, "cases": case_rows}; write_json(OUTPUT / "sentaurus_export_manifest.json", result); return result


def vela_stage_state(workflow: dict[str, Any], stage: str) -> Path:
    root = (REPO / workflow["state"]).parents[1]
    return root / ({"equilibrium": "00_equilibrium/state.csv", "drain": "10_drain_ramp/state.csv", "gate": "20_gate_sweep/state.csv"}[stage])


def analyze(contract: dict[str, Any], exports: dict[str, Any], banner: str) -> tuple[dict[str, Any], dict[str, Any]]:
    vela_manifest, m65_sent = read_json(M65_VELA), read_json(M65_SENT)
    vela = {(row["device"], float(row["drain_voltage_V"])): row for row in vela_manifest["workflows"]}; anchors = {(row["device"], float(row["drain_voltage_V"])): float(row["drain_current_A_per_um"]) for row in m65_sent["states"]}
    state_rows = []
    for item in exports["states"]:
        key = (item["device"], float(item["drain_voltage_V"])); export_dir = REPO / item["export_dir"]
        sent_psi, sent_phin, sent_n = m47.scalar_field(export_dir, "ElectrostaticPotential"), m47.scalar_field(export_dir, "eQuasiFermiPotential"), m47.scalar_field(export_dir, "eDensity")
        vela_state = m65.vela_state(vela_stage_state(vela[key], item["stage"])); coords = m47.coordinates(export_dir)
        mesh = read_json(M8 / "vela" / item["device"] / "mesh.json")
        mesh_coords = {int(row["id"]): (float(row["x"]), float(row["y"])) for row in mesh["nodes"]}
        coordinate_error = max(max(abs(coords[node][axis] - mesh_coords[node][axis]) for axis in (0, 1)) for node in sent_psi)
        _, support = m47.interface_nodes(export_dir); support = sorted(set(support) & set(vela_state["psi"]) & set(sent_psi))
        control = max(support, key=lambda node: sent_phin[node] - sent_psi[node]); dphin = vela_state["phin"][control] - sent_phin[control]; dpsi = vela_state["psi"][control] - sent_psi[control]; dbarrier = dphin - dpsi; proxy = -dbarrier/(VT*math.log(10.0))
        density = math.log10((vela_state["electrons_m3"][control]/1e6)/sent_n[control]); state_rows.append({"device": item["device"], "drain_voltage_V": item["drain_voltage_V"], "gate_voltage_V": item["gate_voltage_V"], "stage": item["stage"], "sentaurus_control_node": control, "support_node_count": len(support), "maximum_coordinate_error_um": coordinate_error, "fixed_node_phin_delta_V": dphin, "fixed_node_psi_delta_V": dpsi, "fixed_node_barrier_delta_V": dbarrier, "fixed_node_qf_barrier_proxy_dex": proxy, "fixed_node_electron_density_log_ratio_dex": density, "boltzmann_density_identity_residual_dex": density-proxy})
    state_by = {(row["device"], float(row["drain_voltage_V"]), row["stage"]): row for row in state_rows}; pair_rows = []
    for pair in read_csv(M65_PAIRS):
        low, high, drain = pair["low_device"], pair["high_device"], float(pair["drain_voltage_V"])
        for stage in STAGES:
            lo, hi = state_by[(low, drain, stage)], state_by[(high, drain, stage)]; value = float(hi["fixed_node_qf_barrier_proxy_dex"])-float(lo["fixed_node_qf_barrier_proxy_dex"])
            pair_rows.append({"low_device": low, "high_device": high, "drain_voltage_V": drain, "diagnostic_gate_voltage_V": pair["diagnostic_gate_voltage_V"], "stage": stage, "qf_barrier_pair_proxy_dex": value})
    pair_by = {(row["low_device"], row["high_device"], float(row["drain_voltage_V"]), row["stage"]): row for row in pair_rows}; increments = []
    for pair in read_csv(M65_PAIRS):
        key = (pair["low_device"], pair["high_device"], float(pair["drain_voltage_V"])); eq = float(pair_by[(*key, "equilibrium")]["qf_barrier_pair_proxy_dex"]); dr = float(pair_by[(*key, "drain")]["qf_barrier_pair_proxy_dex"]); gate = float(pair_by[(*key, "gate")]["qf_barrier_pair_proxy_dex"])
        increments.append({"equilibrium": eq, "drain": dr-eq, "gate": gate-dr, "final": gate})
        for stage, increment in (("equilibrium", eq), ("drain", dr-eq), ("gate", gate-dr)):
            pair_by[(*key, stage)]["increment_from_previous_stage_dex"] = increment
            pair_by[(*key, stage)]["final_gate_proxy_dex"] = gate
    median_abs = {stage: median([abs(row[stage]) for row in increments]) for stage in STAGES}; total = sum(median_abs.values()); shares = {stage: median_abs[stage]/total if total else 0.0 for stage in STAGES}; winner = max(shares, key=shares.get)
    classification = f"{winner}_stage_dominant" if shares[winner] >= float(contract["analysis"]["dominant_stage_minimum_median_absolute_share"]) else "mixed_stage_residual"
    replays = [math.log10(float(row["final_drain_current_A_per_um"])/anchors[(row["device"], float(row["drain_voltage_V"]))]) for row in exports["cases"]]
    checks = {"contract_frozen": read_json(FREEZE)["contract_sha256"] == sha256(CONTRACT), "sentaurus_release": contract["acceptance"]["required_sentaurus_release"] in banner, "case_count": len(exports["cases"]) == int(contract["acceptance"]["required_case_count"]), "state_counts": len(state_rows) == int(contract["acceptance"]["required_state_count_per_solver"]), "pair_stage_count": len(pair_rows) == int(contract["acceptance"]["required_pair_stage_count"]), "final_current_replay": max(abs(value) for value in replays) <= float(contract["acceptance"]["maximum_final_sentaurus_current_replay_error_dex"]), "coordinate_identity": max(float(row["maximum_coordinate_error_um"]) for row in state_rows) <= float(contract["acceptance"]["maximum_coordinate_error_um"]), "boltzmann_density_identity": max(abs(float(row["boltzmann_density_identity_residual_dex"])) for row in state_rows) <= float(contract["acceptance"]["maximum_boltzmann_density_identity_residual_dex"]), "classification_declared": classification in contract["analysis"]["classifications"], "production_reference_not_replaced": True}; checks["all_checks_pass"] = all(checks.values())
    if not checks["all_checks_pass"]: classification = "execution_or_identity_failure"
    report = {"schema": "vela.simplemos.sdevice.m69_stage_residual_localization_report.v2", "status": "accepted" if checks["all_checks_pass"] else "failed", "classification": classification, "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)}, "execution": {"sentaurus_release": banner.splitlines()[0] if banner else "", "new_sentaurus_cases": len(exports["cases"]), "new_vela_solves": 0, "independent_v2_rerun": True, "production_default_changed": False}, "summary": {"median_absolute_stage_increment_dex": median_abs, "median_absolute_stage_share": shares, "dominant_stage": winner, "dominant_stage_share": shares[winner], "median_final_qf_barrier_pair_proxy_dex": median([row["final"] for row in increments]), "maximum_final_sentaurus_current_replay_error_dex": max(abs(value) for value in replays)}, "causal_scope": {"closed": "The stage at which the residual qf-barrier observable appears is localized on the same physical path.", "not_closed": "Stage localization is not by itself a transport closure; M70 is conditional on the remaining unexplained current residual."}, "acceptance": checks}
    return report, {"states": state_rows, "pairs": list(pair_by.values()), "cases": exports["cases"]}


def freeze_results(report: dict[str, Any], rows: dict[str, Any]) -> None:
    write_csv(STATES, rows["states"]); write_csv(PAIRS, rows["pairs"]); write_csv(CASES, rows["cases"]); write_json(REPORT, report); summary=report["summary"]
    DOC.parent.mkdir(parents=True, exist_ok=True); DOC.write_text(f"""# SimpleMOS M69 阶段残差定位

M69 分类为 `{report['classification']}`。在平衡、仅漏压、诊断栅压三个状态上复用同一 noBGN 物理路径，势垒配对增量的中位绝对份额分别为：平衡 `{float(summary['median_absolute_stage_share']['equilibrium']):.2%}`、漏压 `{float(summary['median_absolute_stage_share']['drain']):.2%}`、栅压 `{float(summary['median_absolute_stage_share']['gate']):.2%}`。主导阶段为 `{summary['dominant_stage']}`。

Sentaurus 最终电流相对 M65 的最大重放误差 `{float(summary['maximum_final_sentaurus_current_replay_error_dex']):.3e} dex`。该结论定位状态差异形成阶段，不把势垒代理等同于完整漂移扩散电流闭合。机器报告：`{portable(REPORT)}`。
""", encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {"schema": "vela.validation.artifact.v1", "title": "SimpleMOS M69 stage residual localization", "status": report["status"], "classification": report["classification"], "report": portable(REPORT), "ledgers": [portable(STATES), portable(PAIRS), portable(CASES)]})
    artifacts=[REPORT,STATES,PAIRS,CASES,DOC,ARTIFACT]; write_json(EVIDENCE,{"schema":"vela.simplemos.sdevice.m69_stage_residual_localization_evidence.v2","status":"frozen" if report["acceptance"]["all_checks_pass"] else "failed","classification":report["classification"],"contract_sha256":sha256(CONTRACT),"v1_failure_evidence":portable(M69_V1_EVIDENCE),"implementation_hashes":{portable(SCRIPT):sha256(SCRIPT)},"artifacts":{portable(path):sha256(path) for path in artifacts},"new_sentaurus_execution":True,"new_vela_execution":False,"production_reference_replaced":False,"acceptance":report["acceptance"]})


def verify() -> dict[str, Any]:
    validate_contract(); evidence,report=read_json(EVIDENCE),read_json(REPORT)
    if evidence["status"]!="frozen" or report["status"]!="accepted": raise ValueError("M69 not frozen")
    if evidence["implementation_hashes"][portable(SCRIPT)]!=sha256(SCRIPT): raise ValueError("M69 implementation changed")
    for relative,expected in evidence["artifacts"].items():
        if sha256(REPO/relative)!=expected: raise ValueError(f"M69 artifact changed: {relative}")
    return report


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--freeze-contract",action="store_true"); parser.add_argument("--prepare",action="store_true"); parser.add_argument("--run-sentaurus",action="store_true"); parser.add_argument("--export",action="store_true"); parser.add_argument("--analyze",action="store_true"); parser.add_argument("--verify",action="store_true"); parser.add_argument("--force",action="store_true"); parser.add_argument("--jobs",type=int,default=8); parser.add_argument("--ssh-target",default="sentaurus"); parser.add_argument("--ssh-bin",default=r"C:\Windows\System32\OpenSSH\ssh.exe"); parser.add_argument("--scp-bin",default=r"C:\Windows\System32\OpenSSH\scp.exe"); parser.add_argument("--ssh-config",default=str(SSH_CONFIG)); args=parser.parse_args()
    if args.freeze_contract: freeze_contract(); print(json.dumps({"status":"frozen_before_execution","contract_sha256":sha256(CONTRACT)})); return
    contract=validate_contract()
    if args.verify: print(json.dumps(verify()["summary"],indent=2)); return
    manifest=prepare(args.force) if args.prepare else read_json(OUTPUT/"sentaurus_manifest.json"); banner=""
    if args.prepare and not (args.run_sentaurus or args.export or args.analyze):
        print(json.dumps({"status":"prepared","case_count":len(manifest["cases"]),"contract_sha256":sha256(CONTRACT)}));return
    if args.run_sentaurus: banner=run_sentaurus(manifest,args.ssh_target,args.ssh_bin,args.scp_bin,args.ssh_config,args.jobs)
    elif (OUTPUT/"sentaurus_banner.txt").is_file(): banner=(OUTPUT/"sentaurus_banner.txt").read_text(encoding="utf-8").strip()
    exports=export(manifest) if args.export else read_json(OUTPUT/"sentaurus_export_manifest.json")
    if args.analyze:
        report,rows=analyze(contract,exports,banner); freeze_results(report,rows); print(json.dumps({"status":report["status"],"classification":report["classification"],"summary":report["summary"],"acceptance":report["acceptance"]},indent=2))


if __name__=="__main__": main()
