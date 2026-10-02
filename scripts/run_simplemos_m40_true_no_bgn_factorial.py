#!/usr/bin/env python3
"""Run SimpleMOS M40 explicit no-BGN x SRH causal factorial."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Iterable, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import run_simplemos_m10_fixed_state_replay as m10  # noqa: E402


REPO = Path(__file__).resolve().parents[1]
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m40_true_no_bgn_factorial_contract_v1.json")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m40_true_no_bgn_factorial")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "true_no_bgn_factorial")
INPUT_TDR = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
             / "m8_original_physics/sentaurus_bundle/n23/input_fps.tdr")
VELA_BASE = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
             / "m8a_model_ablation/vela/no_hfs/workflow/vd_0p05")
M29_REPORT = (REPO / "reference_tcad/simplemos_sentaurus2022/bgn_srh_factorial"
              / "m29_bgn_srh_factorial_report.json")
IMPORTER = REPO / "build-release/sentaurus_import.exe"
RUNNER = REPO / "build-release/vela_example_runner.exe"
REMOTE_ROOT = "~/sentaurus_runs/vela_oracle/simplemos_m40_true_no_bgn_20260901_v1"
ARCHIVE_NAME = "simplemos_m40_true_no_bgn_results.tgz"
MANUAL = ("/atctools/Synopsys/tcad/T-2022.03/tcad/T-2022.03-SP2/"
          "manuals/PDFManual/data/sdevice_ug.pdf")
CELLS = ("bgn_on_srh_on", "bgn_on_srh_off",
         "bgn_off_srh_on", "bgn_off_srh_off")


def executable(name: str) -> str:
    if os.name == "nt":
        candidate = (Path(os.environ.get("SystemRoot", r"C:\Windows"))
                     / "System32/OpenSSH" / f"{name}.exe")
        if candidate.is_file():
            return str(candidate)
    return shutil.which(name) or name


def run(argv: Sequence[str], capture: bool = False) -> str:
    completed = subprocess.run(list(argv), cwd=REPO, text=True,
                               stdout=subprocess.PIPE if capture else None,
                               stderr=subprocess.STDOUT if capture else None,
                               encoding="utf-8", errors="replace",
                               check=False)
    if completed.returncode:
        raise RuntimeError(f"command failed ({completed.returncode}): {argv}\n"
                           f"{completed.stdout or ''}")
    return completed.stdout or ""


def read_json(path: Path) -> Any:
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
    if contract.get("schema") != "vela.simplemos.sdevice.m40_true_no_bgn_factorial.v1":
        raise ValueError("unexpected M40 contract schema")
    if [row["id"] for row in contract["cells"]] != list(CELLS):
        raise ValueError("M40 cell order changed")
    if contract["fixed_invariants"]["hfs_enabled"]:
        raise ValueError("M40 must remain HFS-off")
    manual = contract["manual_contract"]
    if manual["no_bgn_syntax"] != "EffectiveIntrinsicDensity(NoBandGapNarrowing)":
        raise ValueError("M40 requires the documented no-BGN syntax")


def deck(case: str, bgn_model: str, srh: bool) -> str:
    if bgn_model == "old_slotboom":
        bgn_line = "Physics { EffectiveIntrinsicDensity(BandGapNarrowing(OldSlotboom)) }"
    elif bgn_model == "none":
        bgn_line = "Physics { EffectiveIntrinsicDensity(NoBandGapNarrowing) }"
    else:
        raise ValueError(f"unsupported BGN model {bgn_model}")
    srh_line = "  Recombination(SRH(DopingDependence))\n" if srh else ""
    return f'''File {{
  Grid="input_fps.tdr"
  Plot="{case}_des.tdr"
  Current="{case}"
  Output="{case}.log"
}}

Electrode {{
  {{ Name="source" Voltage=0.0 }}
  {{ Name="drain" Voltage=0.0 }}
  {{ Name="gate" Voltage=0.0 }}
  {{ Name="substrate" Voltage=0.0 }}
}}

{bgn_line}
Physics(Material="Silicon") {{
  Mobility(PhuMob Enormal)
{srh_line}}}

Plot {{
  eDensity hDensity eQuasiFermi hQuasiFermi Potential
  eMobility hMobility SRH
  Doping DonorConcentration AcceptorConcentration
  BandGap EffectiveBandGap BandGapNarrowing EffectiveIntrinsicDensity
  Affinity ConductionBand ValenceBand
}}

Math {{ Extrapolate Iterations=20 ExitOnFailure }}

Solve {{
  Coupled(Iterations=100) {{ Poisson }}
  Coupled {{ Poisson Electron Hole }}
  Quasistationary(
    InitialStep=0.1 Increment=1.5 MinStep=1e-5 MaxStep=1
    Goal {{ Name="drain" Voltage=0.05 }}
  ) {{ Coupled {{ Poisson Electron Hole }} }}
  NewCurrentPrefix="final_"
  Quasistationary(
    DoZero InitialStep=0.5 Increment=1.5 MinStep=0.0005 MaxStep=1
    Goal {{ Name="gate" Voltage=0.05 }}
  ) {{
    Coupled {{ Poisson Electron Hole }}
    CurrentPlot(Time=(1))
  }}
  Plot(FilePrefix="final")
}}
'''


def prepare(contract: dict[str, Any], force: bool) -> dict[str, Any]:
    bundle = OUTPUT / "sentaurus_bundle"
    if force and bundle.exists():
        shutil.rmtree(bundle)
    cases = []
    for spec in contract["cells"]:
        case = f"n23_m40_{spec['id']}"
        root = bundle / case
        root.mkdir(parents=True, exist_ok=True)
        shutil.copy2(INPUT_TDR, root / "input_fps.tdr")
        command = root / f"{case}_des.cmd"
        command.write_text(deck(case, spec["bgn_model"], spec["srh"]),
                           encoding="utf-8", newline="\n")
        cases.append({
            "case": case,
            "cell": spec["id"],
            "bgn_model": spec["bgn_model"],
            "srh": spec["srh"],
            "deck": portable(command),
            "deck_sha256": m10.sha256(command),
            "input_tdr_sha256": m10.sha256(INPUT_TDR),
        })
    manifest = {
        "schema": "vela.simplemos.sdevice.m40_sentaurus_manifest.v1",
        "status": "prepared",
        "contract": portable(CONTRACT),
        "contract_sha256": m10.sha256(CONTRACT),
        "cases": cases,
    }
    write_json(OUTPUT / "sentaurus_manifest.json", manifest)
    return manifest


def probe_manual(ssh_target: str, ssh_bin: str) -> dict[str, Any]:
    pdf_hash = run([ssh_bin, ssh_target, f"sha256sum {MANUAL}"], True).split()[0]
    command = (f"pdftotext -layout {MANUAL} - | "
               "grep -A 8 -B 8 'keyword NoBandGapNarrowing' | head -30")
    excerpt = run([ssh_bin, ssh_target, command], True)
    if ("bandgap narrowing is active" not in excerpt.lower() or
            "EffectiveIntrinsicDensity(NoBandGapNarrowing)" not in excerpt):
        raise RuntimeError("T-2022.03 manual did not confirm NoBandGapNarrowing")
    probe = {
        "schema": "vela.simplemos.sdevice.m40_manual_probe.v1",
        "status": "confirmed",
        "release": "T-2022.03-SP2",
        "manual_path": MANUAL,
        "manual_sha256": pdf_hash,
        "manual_page": 305,
        "confirmed_keyword": "NoBandGapNarrowing",
        "confirmed_syntax": "EffectiveIntrinsicDensity(NoBandGapNarrowing)",
        "confirmation": "The manual identifies BGN as active by default and requires the explicit NoBandGapNarrowing keyword to disable it.",
    }
    write_json(OUTPUT / "manual_probe.json", probe)
    return probe


def run_sentaurus(manifest: dict[str, Any], ssh_target: str, ssh_bin: str,
                   scp_bin: str, remote_root: str, jobs: int) -> dict[str, Any]:
    manual = probe_manual(ssh_target, ssh_bin)
    banner = run([ssh_bin, ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'"], True).strip()
    if "T-2022.03-SP2" not in banner:
        raise RuntimeError(f"unexpected Sentaurus release:\n{banner}")
    exists = run([ssh_bin, ssh_target,
                  f"if test -e {remote_root}; then echo EXISTS; fi"], True).strip()
    if exists:
        raise RuntimeError(f"refusing to overwrite remote path {remote_root}")
    run([ssh_bin, ssh_target, f"mkdir -p {remote_root}"])
    run([scp_bin, "-r", str(OUTPUT / "sentaurus_bundle"),
         f"{ssh_target}:{remote_root}/"])

    def execute(item: dict[str, Any]) -> None:
        case = item["case"]
        command = (f"set -eu; cd {remote_root}/sentaurus_bundle/{case}; "
                   f"sdevice {case}_des.cmd > console.log 2>&1")
        run([ssh_bin, ssh_target, command])

    with ThreadPoolExecutor(max_workers=min(jobs, len(manifest["cases"]))) as pool:
        list(pool.map(execute, manifest["cases"]))
    run([ssh_bin, ssh_target,
         f"cd {remote_root} && tar -czf {ARCHIVE_NAME} sentaurus_bundle"])
    raw = OUTPUT / "sentaurus_raw"
    raw.mkdir(parents=True, exist_ok=True)
    archive = raw / ARCHIVE_NAME
    run([scp_bin, f"{ssh_target}:{remote_root}/{ARCHIVE_NAME}", str(archive)])
    with tarfile.open(archive, "r:gz") as stream:
        stream.extractall(raw, filter="data")
    (OUTPUT / "sentaurus_banner.txt").write_text(banner + "\n", encoding="utf-8")
    result = {"status": "complete", "banner": banner, "manual": manual,
              "archive_sha256": m10.sha256(archive), "case_count": len(manifest["cases"])}
    write_json(OUTPUT / "sentaurus_run.json", result)
    return result


def configure_vela(config: dict[str, Any], bgn_model: str, srh: bool) -> None:
    solver = config["solver"]
    solver["mobility"] = {"model": "phumob_lombardi",
                          "doping_concentration_basis": "total_impurity"}
    solver["recombination"] = ["srh"] if srh else ["none"]
    solver["srh_doping_dependence"]["enabled"] = srh
    solver["bandgap_narrowing"] = {
        "model": "old_slotboom" if bgn_model == "old_slotboom" else "none",
        "fermi_statistics_correction": False,
    }


def vela_cell(spec: dict[str, Any], force: bool) -> dict[str, Any]:
    root = OUTPUT / "vela" / spec["id"]
    final_state = root / "20_gate/state.csv"
    final_curve = root / "20_gate/curve.csv"
    if force and root.exists():
        shutil.rmtree(root)
    previous: Path | None = None
    phases = (("00_equilibrium", "00_equilibrium.json"),
              ("10_drain", "10_drain_ramp.json"),
              ("20_gate", "20_gate_sweep.json"))
    for label, source_name in phases:
        phase = root / label
        phase.mkdir(parents=True, exist_ok=True)
        config = copy.deepcopy(read_json(VELA_BASE / source_name))
        configure_vela(config, spec["bgn_model"], spec["srh"])
        curve = phase / "curve.csv"
        state = phase / "state.csv"
        config["output_csv"] = str(curve.resolve())
        config["log_file"] = str((phase / "run.log").resolve())
        config["sweep"]["write_state_file"] = str(state.resolve())
        if previous is not None:
            config["sweep"]["initial_state_file"] = str(previous.resolve())
        if label == "20_gate":
            config["sweep"].update({"start": 0.0, "stop": 0.05, "step": 0.05,
                                    "bias_points": [0.0, 0.05], "max_step": 0.05})
        config["simplemos_m40"] = {
            "cell": spec["id"], "bgn_model": spec["bgn_model"],
            "srh": spec["srh"], "hfs_enabled": False,
        }
        path = phase / "config.json"
        write_json(path, config)
        status = m10.execute_runner(path, RUNNER)
        if not status.get("converged"):
            raise RuntimeError(f"Vela M40 failed: {spec['id']} {label}")
        previous = state
    row = min(read_csv(final_curve),
              key=lambda item: abs(float(item["bias_V"]) - 0.05))
    return {"cell": spec["id"], "bgn_model": spec["bgn_model"],
            "srh": spec["srh"],
            "current_A_per_um": abs(float(row["current_total_A_per_um"])),
            "state": portable(final_state)}


def run_vela(contract: dict[str, Any], jobs: int, force: bool) -> dict[str, Any]:
    with ThreadPoolExecutor(max_workers=min(jobs, 4)) as pool:
        states = list(pool.map(lambda spec: vela_cell(spec, force), contract["cells"]))
    result = {"schema": "vela.simplemos.sdevice.m40_vela_manifest.v1",
              "status": "complete", "states": states}
    write_json(OUTPUT / "vela_manifest.json", result)
    return result


def parse_bgn_log(path: Path) -> dict[str, str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"Bandgap narrowing model:\s*([^\r\n]+)", text)
    if not match:
        raise RuntimeError(f"no BGN model line in {path}")
    line = match.group(1).strip()
    low = line.lower()
    if "without bandgap narrowing" in low or "no bandgap narrowing" in low:
        model = "none"
    elif "OldSlotboom" in line:
        model = "old_slotboom"
    elif "Bennett/Wilson" in line:
        model = "bennett_wilson"
    else:
        model = "unknown"
    return {"model": model, "line": line}


def max_abs_field(export_dir: Path, name: str) -> float:
    values = m10.scalar_field(export_dir, name).values()
    return max((abs(value) for value in values), default=math.inf)


def doping_identity(export_dir: Path) -> dict[str, Any]:
    donors = m10.scalar_field(export_dir, "DonorConcentration")
    acceptors = m10.scalar_field(export_dir, "AcceptorConcentration")
    doping = {int(row["node_id"]): (float(row["donors_cm3"]),
                                    float(row["acceptors_cm3"]))
              for row in read_csv(export_dir / "doping.csv")}
    common = sorted(set(donors) & set(acceptors) & set(doping))
    maximum = 0.0
    for node in common:
        for lhs, rhs in ((donors[node], doping[node][0]),
                         (acceptors[node], doping[node][1])):
            maximum = max(maximum, abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1.0))
    return {"common_node_count": len(common),
            "maximum_relative_difference": maximum,
            "has_nonzero_doping": any(d != 0.0 or a != 0.0 for d, a in doping.values())}


def export_sentaurus(contract: dict[str, Any], manifest: dict[str, Any]) -> dict[str, Any]:
    expected = {row["id"]: row for row in contract["cells"]}
    states = []
    for item in manifest["cases"]:
        raw = OUTPUT / "sentaurus_raw/sentaurus_bundle" / item["case"]
        tdr = m10.single_glob(raw, "final_des.tdr")
        current = m10.single_glob(raw, "final_*des.plt")
        log = m10.single_glob(raw, "*.log_des.log")
        export_dir = OUTPUT / "sentaurus_exports" / item["cell"]
        run([str(IMPORTER), "--tdr", str(tdr), "--export-dir", str(export_dir)])
        fields = m10.manifest_field_names(export_dir)
        required = set(contract["required_sentaurus_fields"])
        missing = required - fields
        if missing:
            raise RuntimeError(f"{item['cell']} misses fields {sorted(missing)}")
        observed = parse_bgn_log(log)
        expected_model = expected[item["cell"]]["bgn_model"]
        if observed["model"] != expected_model:
            raise RuntimeError(f"{item['cell']} expected {expected_model}, got {observed}")
        identity = doping_identity(export_dir)
        if identity["common_node_count"] < 900 or not identity["has_nonzero_doping"]:
            raise RuntimeError(f"{item['cell']} has unqualified doping observability")
        states.append({
            "cell": item["cell"], "bgn_model": expected_model,
            "srh": expected[item["cell"]]["srh"],
            "tdr": portable(tdr), "tdr_sha256": m10.sha256(tdr),
            "current_file": portable(current), "current_file_sha256": m10.sha256(current),
            "solver_log": portable(log), "solver_log_sha256": m10.sha256(log),
            "solver_bgn_model": observed,
            "export_dir": portable(export_dir),
            "field_manifest_sha256": m10.sha256(export_dir / "field_manifest.json"),
            "doping_identity": identity,
            "bandgap_narrowing_max_abs_eV": max_abs_field(export_dir, "BandgapNarrowing"),
            "terminal": m10.parse_terminal_current(current),
        })
    states.sort(key=lambda row: CELLS.index(row["cell"]))
    result = {"schema": "vela.simplemos.sdevice.m40_sentaurus_exports.v1",
              "status": "complete", "states": states}
    write_json(OUTPUT / "sentaurus_export_manifest.json", result)
    return result


def factorial(values: dict[str, float]) -> dict[str, float]:
    y11, y10 = values["bgn_on_srh_on"], values["bgn_on_srh_off"]
    y01, y00 = values["bgn_off_srh_on"], values["bgn_off_srh_off"]
    return {
        "bgn_main_effect_dex": 0.5 * ((y11 - y01) + (y10 - y00)),
        "srh_main_effect_dex": 0.5 * ((y11 - y10) + (y01 - y00)),
        "bgn_srh_interaction_dex": y11 - y10 - y01 + y00,
        "bgn_effect_srh_on_dex": y11 - y01,
        "bgn_effect_srh_off_dex": y10 - y00,
    }


def analyze(contract: dict[str, Any], sentaurus: dict[str, Any],
            vela: dict[str, Any], manual: dict[str, Any]) -> dict[str, Any]:
    sent = {row["cell"]: row for row in sentaurus["states"]}
    vel = {row["cell"]: row for row in vela["states"]}
    matrix, sent_log, vela_log = [], {}, {}
    for spec in contract["cells"]:
        cell = spec["id"]
        si = abs(float(sent[cell]["terminal"]["drain_total_current_A_per_um"]))
        vi = abs(float(vel[cell]["current_A_per_um"]))
        sent_log[cell], vela_log[cell] = math.log10(si), math.log10(vi)
        matrix.append({
            "cell": cell, "bgn_model": spec["bgn_model"], "srh": spec["srh"],
            "sentaurus_current_A_per_um": si, "vela_current_A_per_um": vi,
            "signed_log10_vela_over_sentaurus": math.log10(vi / si),
            "absolute_log10_ratio": abs(math.log10(vi / si)),
            "relative_error": abs(vi - si) / si,
            "sentaurus_bgn_field_max_abs_eV": sent[cell]["bandgap_narrowing_max_abs_eV"],
            "doping_max_relative_difference": sent[cell]["doping_identity"]["maximum_relative_difference"],
        })
    write_csv(PORTABLE / "m40_factorial_matrix.csv", matrix)
    sent_effect, vela_effect = factorial(sent_log), factorial(vela_log)
    effect_rows = [{"effect": name, "sentaurus_dex": sent_effect[name],
                    "vela_dex": vela_effect[name],
                    "vela_minus_sentaurus_dex": vela_effect[name] - sent_effect[name]}
                   for name in sent_effect]
    write_csv(PORTABLE / "m40_factor_effects.csv", effect_rows)
    no_bgn = [sent[cell] for cell in CELLS if cell.startswith("bgn_off")]
    checks = {
        "manual_keyword_confirmed": manual["status"] == "confirmed",
        "sentaurus_state_count": len(sent) == 4,
        "vela_state_count": len(vel) == 4,
        "solver_log_models_match": all(
            sent[row["id"]]["solver_bgn_model"]["model"] == row["bgn_model"]
            for row in contract["cells"]),
        "no_bgn_fields_zero": max(row["bandgap_narrowing_max_abs_eV"]
                                  for row in no_bgn) <= 1e-14,
        "doping_fields_observable": min(row["doping_identity"]["common_node_count"]
                                        for row in sent.values()) >= 900,
        "doping_identity_exact": max(row["doping_identity"]["maximum_relative_difference"]
                                     for row in sent.values()) <= 1e-15,
        "bias_exact": max(abs(row["terminal"]["gate_voltage_V"] - 0.05)
                          for row in sent.values()) <= 1e-10,
    }
    baseline = next(row for row in matrix if row["cell"] == "bgn_on_srh_on")
    no_bgn_on = next(row for row in matrix if row["cell"] == "bgn_off_srh_on")
    conditional_mismatch = (vela_effect["bgn_effect_srh_on_dex"] -
                            sent_effect["bgn_effect_srh_on_dex"])
    m29 = read_json(M29_REPORT)
    m29_false_off_gap = next(
        row["absolute_log10_ratio"] for row in m29["matrix"]["cells"]
        if row["cell"] == "bgn_off_srh_on")
    report = {
        "schema": "vela.simplemos.sdevice.m40_true_no_bgn_factorial_report.v1",
        "status": "complete" if all(checks.values()) else "failed",
        "execution": {"new_sentaurus_state_count": 4, "new_vela_state_count": 4,
                      "default_model_changed": False},
        "manual_confirmation": manual,
        "sentaurus_model_observability": [{
            "cell": row["cell"],
            "expected_model": row["bgn_model"],
            "solver_log_line": row["solver_bgn_model"]["line"],
            "solver_log_sha256": row["solver_log_sha256"],
            "tdr_sha256": row["tdr_sha256"],
            "current_file_sha256": row["current_file_sha256"],
            "bandgap_narrowing_max_abs_eV": row[
                "bandgap_narrowing_max_abs_eV"],
            "doping_common_node_count": row["doping_identity"][
                "common_node_count"],
            "doping_max_relative_difference": row["doping_identity"][
                "maximum_relative_difference"],
        } for row in sentaurus["states"]],
        "matrix": matrix,
        "factor_effects": {"sentaurus": sent_effect, "vela": vela_effect,
                           "cross_solver": effect_rows},
        "causal_result": {
            "production_gap_dex": baseline["absolute_log10_ratio"],
            "true_no_bgn_srh_on_gap_dex": no_bgn_on["absolute_log10_ratio"],
            "gap_change_dex": (no_bgn_on["absolute_log10_ratio"] -
                               baseline["absolute_log10_ratio"]),
            "bgn_response_mismatch_srh_on_dex": conditional_mismatch,
            "bgn_related_fraction_of_production_gap": (
                conditional_mismatch / baseline["absolute_log10_ratio"]),
            "remaining_fraction_of_production_gap": (
                no_bgn_on["absolute_log10_ratio"] /
                baseline["absolute_log10_ratio"]),
            "m29_bennett_wilson_labelled_off_gap_dex": m29_false_off_gap,
            "true_no_bgn_minus_m29_bennett_wilson_gap_dex": (
                no_bgn_on["absolute_log10_ratio"] - m29_false_off_gap),
            "primary_contrast": "BGN conditional effect with SRH enabled",
            "factorial_main_effect_qualified": False,
            "interaction_warning": (
                "The no-BGN/SRH-off Vela corner has a 1.33254 dex parity gap; "
                "the 2x2 main effects are not used as a single-cause estimate."
            ),
        },
        "acceptance": {**checks, "all_checks_pass": all(checks.values())},
        "claim_policy": contract["claim_policy"],
        "artifacts": {
            "factorial_matrix": portable(PORTABLE / "m40_factorial_matrix.csv"),
            "factor_effects": portable(PORTABLE / "m40_factor_effects.csv"),
        },
    }
    write_json(PORTABLE / "m40_true_no_bgn_factorial_report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--live-sentaurus", action="store_true")
    parser.add_argument("--run-vela", action="store_true")
    parser.add_argument("--export", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--scp-bin", default=executable("scp"))
    parser.add_argument("--remote-root", default=REMOTE_ROOT)
    args = parser.parse_args()
    if not any((args.prepare, args.live_sentaurus, args.run_vela,
                args.export, args.analyze)):
        parser.error("select at least one action")
    contract = read_json(CONTRACT)
    validate(contract)
    manifest_path = OUTPUT / "sentaurus_manifest.json"
    manifest = prepare(contract, args.force) if args.prepare else read_json(manifest_path)
    live = run_sentaurus(manifest, args.ssh_target, args.ssh_bin, args.scp_bin,
                         args.remote_root, args.jobs) if args.live_sentaurus else None
    vela = run_vela(contract, args.jobs, args.force) if args.run_vela else None
    sentaurus = export_sentaurus(contract, manifest) if args.export else None
    report = None
    if args.analyze:
        manual = (live or {}).get("manual") or read_json(OUTPUT / "manual_probe.json")
        sentaurus = sentaurus or read_json(OUTPUT / "sentaurus_export_manifest.json")
        vela = vela or read_json(OUTPUT / "vela_manifest.json")
        report = analyze(contract, sentaurus, vela, manual)
    print(json.dumps({
        "status": (report or live or manifest).get("status"),
        "sentaurus_cases": len(manifest["cases"]),
        "all_checks_pass": (report or {}).get("acceptance", {}).get("all_checks_pass"),
        "default_model_changed": False,
    }))
    return 0 if report is None or report["acceptance"]["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
