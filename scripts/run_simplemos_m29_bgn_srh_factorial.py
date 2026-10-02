#!/usr/bin/env python3
"""Run the SimpleMOS M29 HFS-off BGN x SRH causal factorial."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Iterable


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import run_simplemos_m10_fixed_state_replay as m10  # noqa: E402


REPO = Path(__file__).resolve().parents[1]
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m29_bgn_srh_factorial_contract_v1.json")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m29_bgn_srh_factorial")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "bgn_srh_factorial")
INPUT_TDR = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
             / "m8_original_physics/sentaurus_bundle/n23/input_fps.tdr")
M28 = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
       / "m28_exact_bias_path_audit")
VELA_BASE = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
             / "m8a_model_ablation/vela/no_hfs/workflow/vd_0p05")
IMPORTER = REPO / "build-release/sentaurus_import.exe"
RUNNER = REPO / "build-release/vela_example_runner.exe"
CELLS = ("bgn_on_srh_on", "bgn_on_srh_off",
         "bgn_off_srh_on", "bgn_off_srh_off")


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


def run(command: list[str]) -> str:
    completed = subprocess.run(command, cwd=REPO, text=True,
                               capture_output=True, check=False)
    if completed.returncode:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {command}\n"
            f"{completed.stderr or completed.stdout}")
    return completed.stdout


def validate(contract: dict[str, Any]) -> None:
    if contract.get("schema") != "vela.simplemos.sdevice.m29_bgn_srh_factorial.v1":
        raise ValueError("unexpected M29 contract schema")
    if [row["id"] for row in contract["cells"]] != list(CELLS):
        raise ValueError("M29 cell order changed")
    if contract["fixed_invariants"]["hfs_enabled"]:
        raise ValueError("M29 must remain HFS-off")


def deck(case: str, bgn: bool, srh: bool) -> str:
    bgn_line = "Physics { EffectiveIntrinsicDensity(OldSlotboom) }\n" if bgn else ""
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

{bgn_line}Physics(Material="Silicon") {{
  Mobility(PhuMob Enormal)
{srh_line}}}

Plot {{
  eDensity hDensity eQuasiFermi hQuasiFermi Potential
  eMobility hMobility SRH BandGap BandGapNarrowing
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
        if spec["reuse_m28"]:
            continue
        case = f"n23_m29_{spec['id']}"
        root = bundle / case
        root.mkdir(parents=True, exist_ok=True)
        shutil.copy2(INPUT_TDR, root / "input_fps.tdr")
        command = root / f"{case}_des.cmd"
        command.write_text(deck(case, spec["bgn"], spec["srh"]),
                           encoding="utf-8", newline="\n")
        cases.append({
            "case": case,
            "cell": spec["id"],
            "bgn": spec["bgn"],
            "srh": spec["srh"],
            "deck": portable(command),
            "deck_sha256": m10.sha256(command),
            "input_tdr_sha256": m10.sha256(INPUT_TDR),
        })
    manifest = {
        "schema": "vela.simplemos.sdevice.m29_sentaurus_manifest.v1",
        "status": "prepared",
        "contract": portable(CONTRACT),
        "contract_sha256": m10.sha256(CONTRACT),
        "cases": cases,
    }
    write_json(OUTPUT / "sentaurus_manifest.json", manifest)
    return manifest


def configure_physics(config: dict[str, Any], bgn: bool, srh: bool) -> None:
    solver = config["solver"]
    solver["mobility"] = {
        "model": "phumob_lombardi",
        "doping_concentration_basis": "total_impurity",
    }
    solver["recombination"] = ["srh"] if srh else ["none"]
    solver["srh_doping_dependence"]["enabled"] = srh
    solver["bandgap_narrowing"] = {
        "model": "old_slotboom" if bgn else "none",
        "fermi_statistics_correction": False,
    }


def vela_cell(spec: dict[str, Any], force: bool) -> dict[str, Any]:
    root = OUTPUT / "vela" / spec["id"]
    final_state = root / "20_gate/state.csv"
    final_curve = root / "20_gate/curve.csv"
    if final_state.is_file() and final_curve.is_file() and not force:
        row = min(read_csv(final_curve),
                  key=lambda item: abs(float(item["bias_V"]) - 0.05))
        return {"cell": spec["id"], "bgn": spec["bgn"], "srh": spec["srh"],
                "current_A_per_um": abs(float(row["current_total_A_per_um"])),
                "state": portable(final_state)}
    previous: Path | None = None
    phases = (("00_equilibrium", "00_equilibrium.json"),
              ("10_drain", "10_drain_ramp.json"),
              ("20_gate", "20_gate_sweep.json"))
    for label, source_name in phases:
        phase = root / label
        phase.mkdir(parents=True, exist_ok=True)
        config = read_json(VELA_BASE / source_name)
        configure_physics(config, spec["bgn"], spec["srh"])
        curve = phase / "curve.csv"
        state = phase / "state.csv"
        config["output_csv"] = str(curve.resolve())
        config["log_file"] = str((phase / "run.log").resolve())
        config["sweep"]["write_state_file"] = str(state.resolve())
        if previous is not None:
            config["sweep"]["initial_state_file"] = str(previous.resolve())
        if label == "20_gate":
            config["sweep"].update({
                "start": 0.0, "stop": 0.05, "step": 0.05,
                "bias_points": [0.0, 0.05], "max_step": 0.05,
            })
        config["simplemos_m29"] = {
            "cell": spec["id"], "bgn": spec["bgn"], "srh": spec["srh"],
            "hfs_enabled": False, "fixed_route": True,
        }
        path = phase / "config.json"
        write_json(path, config)
        status = m10.execute_runner(path, RUNNER)
        if not status.get("converged"):
            raise RuntimeError(f"Vela M29 failed: {spec['id']} {label}")
        previous = state
    rows = read_csv(final_curve)
    row = min(rows, key=lambda item: abs(float(item["bias_V"]) - 0.05))
    if not math.isclose(float(row["bias_V"]), 0.05, abs_tol=1e-12):
        raise ValueError(f"missing exact Vela bias for {spec['id']}")
    return {"cell": spec["id"], "bgn": spec["bgn"], "srh": spec["srh"],
            "current_A_per_um": abs(float(row["current_total_A_per_um"])),
            "state": portable(final_state)}


def run_vela(contract: dict[str, Any], jobs: int, force: bool) -> dict[str, Any]:
    with ThreadPoolExecutor(max_workers=min(jobs, 4)) as pool:
        states = list(pool.map(lambda spec: vela_cell(spec, force), contract["cells"]))
    result = {"schema": "vela.simplemos.sdevice.m29_vela_manifest.v1",
              "status": "complete", "states": states}
    write_json(OUTPUT / "vela_manifest.json", result)
    return result


def export_sentaurus(contract: dict[str, Any], manifest: dict[str, Any]) -> dict[str, Any]:
    states = []
    baseline_root = M28 / "sentaurus_raw/sentaurus_bundle/n23_original_initial_dozero_no_hfs"
    baseline_export = M28 / "sentaurus_exports/original_initial_dozero/no_hfs"
    baseline_current = m10.single_glob(baseline_root, "final_*des.plt")
    states.append({
        "cell": "bgn_on_srh_on", "bgn": True, "srh": True, "reused_m28": True,
        "tdr": portable(m10.single_glob(baseline_root, "final_des.tdr")),
        "current_file": portable(baseline_current),
        "export_dir": portable(baseline_export),
        "terminal": m10.parse_terminal_current(baseline_current),
    })
    for item in manifest["cases"]:
        raw = OUTPUT / "sentaurus_raw/sentaurus_bundle" / item["case"]
        tdr = m10.single_glob(raw, "final_des.tdr")
        current = m10.single_glob(raw, "final_*des.plt")
        export_dir = OUTPUT / "sentaurus_exports" / item["cell"]
        run([str(IMPORTER), "--tdr", str(tdr), "--export-dir", str(export_dir)])
        required = {"ElectrostaticPotential", "eQuasiFermiPotential", "eDensity"}
        missing = required - m10.manifest_field_names(export_dir)
        if missing:
            raise RuntimeError(f"{item['cell']} misses {sorted(missing)}")
        states.append({
            "cell": item["cell"], "bgn": item["bgn"], "srh": item["srh"],
            "reused_m28": False, "tdr": portable(tdr),
            "tdr_sha256": m10.sha256(tdr), "current_file": portable(current),
            "current_file_sha256": m10.sha256(current),
            "export_dir": portable(export_dir),
            "field_manifest_sha256": m10.sha256(export_dir / "field_manifest.json"),
            "terminal": m10.parse_terminal_current(current),
        })
    states.sort(key=lambda row: CELLS.index(row["cell"]))
    result = {"schema": "vela.simplemos.sdevice.m29_sentaurus_exports.v1",
              "status": "complete", "states": states}
    write_json(OUTPUT / "sentaurus_export_manifest.json", result)
    return result


def percentile(values: Iterable[float], fraction: float) -> float:
    return float(m10.percentile(list(values), fraction))


def pearson(left: list[float], right: list[float]) -> float:
    lm, rm = sum(left) / len(left), sum(right) / len(right)
    numerator = sum((a - lm) * (b - rm) for a, b in zip(left, right))
    denominator = math.sqrt(sum((a - lm) ** 2 for a in left)
                            * sum((b - rm) ** 2 for b in right))
    return numerator / denominator if denominator else 0.0


def state_field(path: Path, name: str) -> dict[int, float]:
    return {int(row["node_id"]): float(row[name]) for row in read_csv(path)}


def factorial(values: dict[str, float]) -> dict[str, float]:
    y11 = values["bgn_on_srh_on"]
    y10 = values["bgn_on_srh_off"]
    y01 = values["bgn_off_srh_on"]
    y00 = values["bgn_off_srh_off"]
    return {
        "bgn_main_effect_dex": 0.5 * ((y11 - y01) + (y10 - y00)),
        "srh_main_effect_dex": 0.5 * ((y11 - y10) + (y01 - y00)),
        "bgn_srh_interaction_dex": y11 - y10 - y01 + y00,
        "bgn_effect_srh_on_dex": y11 - y01,
        "bgn_effect_srh_off_dex": y10 - y00,
        "srh_effect_bgn_on_dex": y11 - y10,
        "srh_effect_bgn_off_dex": y01 - y00,
    }


def analyze(contract: dict[str, Any], sentaurus: dict[str, Any],
            vela: dict[str, Any]) -> dict[str, Any]:
    sent = {row["cell"]: row for row in sentaurus["states"]}
    vel = {row["cell"]: row for row in vela["states"]}
    matrix_rows = []
    sent_log, vela_log = {}, {}
    for spec in contract["cells"]:
        cell = spec["id"]
        si = abs(float(sent[cell]["terminal"]["drain_total_current_A_per_um"]))
        vi = abs(float(vel[cell]["current_A_per_um"]))
        sent_log[cell], vela_log[cell] = math.log10(si), math.log10(vi)
        matrix_rows.append({
            "cell": cell, "bgn": spec["bgn"], "srh": spec["srh"],
            "sentaurus_current_A_per_um": si, "vela_current_A_per_um": vi,
            "signed_log10_vela_over_sentaurus": math.log10(vi / si),
            "absolute_log10_ratio": abs(math.log10(vi / si)),
            "relative_error": abs(vi - si) / si,
        })
    write_csv(PORTABLE / "m29_factorial_matrix.csv", matrix_rows)
    sent_effect = factorial(sent_log)
    vela_effect = factorial(vela_log)
    effect_rows = []
    for name in sent_effect:
        effect_rows.append({
            "effect": name,
            "sentaurus_dex": sent_effect[name],
            "vela_dex": vela_effect[name],
            "vela_minus_sentaurus_dex": vela_effect[name] - sent_effect[name],
        })
    write_csv(PORTABLE / "m29_factor_effects.csv", effect_rows)

    sent_fields, vela_fields = {}, {}
    for cell in CELLS:
        sent_root = REPO / sent[cell]["export_dir"]
        sent_fields[cell] = {
            "phin": m10.scalar_field(sent_root, "eQuasiFermiPotential"),
            "logn": {node: math.log10(max(value, 1e-300))
                     for node, value in m10.scalar_field(sent_root, "eDensity").items()},
        }
        state = REPO / vel[cell]["state"]
        vela_fields[cell] = {
            "phin": state_field(state, "phin"),
            "logn": {node: math.log10(max(value, 1e-300))
                     for node, value in state_field(state, "electrons_m3").items()},
        }
    node_rows = []
    for factor in ("bgn", "srh"):
        for field in ("phin", "logn"):
            common = sorted(set.intersection(*(
                set(sent_fields[cell][field]) for cell in CELLS)))
            for node in common:
                def effect(source: dict[str, dict[str, dict[int, float]]]) -> float:
                    v = {cell: source[cell][field][node] for cell in CELLS}
                    if factor == "bgn":
                        return 0.5 * ((v["bgn_on_srh_on"] - v["bgn_off_srh_on"])
                                      + (v["bgn_on_srh_off"] - v["bgn_off_srh_off"]))
                    return 0.5 * ((v["bgn_on_srh_on"] - v["bgn_on_srh_off"])
                                  + (v["bgn_off_srh_on"] - v["bgn_off_srh_off"]))
                node_rows.append({"factor": factor, "field": field, "node_id": node,
                                  "sentaurus_effect": effect(sent_fields),
                                  "vela_effect": effect(vela_fields)})
    write_csv(PORTABLE / "m29_node_factor_effects.csv", node_rows)
    spatial = []
    for factor in ("bgn", "srh"):
        for field in ("phin", "logn"):
            rows = [row for row in node_rows
                    if row["factor"] == factor and row["field"] == field]
            left = [float(row["sentaurus_effect"]) for row in rows]
            right = [float(row["vela_effect"]) for row in rows]
            spatial.append({
                "factor": factor, "field": field, "node_count": len(rows),
                "pearson": pearson(left, right),
                "median_absolute_difference": percentile(
                    [abs(a - b) for a, b in zip(left, right)], 0.5),
                "p95_absolute_difference": percentile(
                    [abs(a - b) for a, b in zip(left, right)], 0.95),
            })
    write_csv(PORTABLE / "m29_spatial_effect_summary.csv", spatial)

    baseline_gap = next(row for row in matrix_rows
                        if row["cell"] == "bgn_on_srh_on")["absolute_log10_ratio"]
    best = min(matrix_rows, key=lambda row: row["absolute_log10_ratio"])
    values = [float(value) for row in matrix_rows for key, value in row.items()
              if key not in ("cell", "bgn", "srh")]
    report = {
        "schema": "vela.simplemos.sdevice.m29_bgn_srh_factorial_report.v1",
        "status": "complete",
        "execution": {"new_sentaurus_state_count": 3, "vela_state_count": 4,
                      "reused_m28_state_count": 1, "default_model_changed": False},
        "acceptance": {
            "maximum_bias_error_V": max(abs(float(row["terminal"]["gate_voltage_V"]) - 0.05)
                                         for row in sentaurus["states"]),
            "nonfinite_fraction": sum(not math.isfinite(value) for value in values) / len(values),
            "common_silicon_node_count": min(row["node_count"] for row in spatial),
        },
        "matrix": {"cells": matrix_rows, "baseline_absolute_gap_dex": baseline_gap,
                   "best_gap_cell": best["cell"],
                   "best_absolute_gap_dex": best["absolute_log10_ratio"]},
        "factor_effects": {"sentaurus": sent_effect, "vela": vela_effect,
                           "cross_solver": effect_rows},
        "spatial_effects": spatial,
        "artifacts": {
            "factorial_matrix": portable(PORTABLE / "m29_factorial_matrix.csv"),
            "factor_effects": portable(PORTABLE / "m29_factor_effects.csv"),
            "node_factor_effects": portable(PORTABLE / "m29_node_factor_effects.csv"),
            "spatial_effect_summary": portable(PORTABLE / "m29_spatial_effect_summary.csv"),
        },
        "claim_policy": contract["claim_policy"],
    }
    checks = {
        "sentaurus_state_count": len(sentaurus["states"]) == 4,
        "vela_state_count": len(vela["states"]) == 4,
        "bias_error": report["acceptance"]["maximum_bias_error_V"] <= 1e-10,
        "finite_values": report["acceptance"]["nonfinite_fraction"] == 0.0,
        "common_nodes": report["acceptance"]["common_silicon_node_count"] >= 900,
    }
    report["acceptance"]["checks"] = checks
    report["acceptance"]["all_checks_pass"] = all(checks.values())
    write_json(PORTABLE / "m29_bgn_srh_factorial_report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--run-vela", action="store_true")
    parser.add_argument("--export", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    contract = read_json(CONTRACT)
    validate(contract)
    manifest = prepare(contract, args.force) if args.prepare else read_json(
        OUTPUT / "sentaurus_manifest.json")
    vela = run_vela(contract, args.jobs, args.force) if args.run_vela else None
    sentaurus = export_sentaurus(contract, manifest) if args.export else None
    if args.analyze:
        vela = vela or read_json(OUTPUT / "vela_manifest.json")
        sentaurus = sentaurus or read_json(OUTPUT / "sentaurus_export_manifest.json")
        report = analyze(contract, sentaurus, vela)
        print(json.dumps({"status": report["status"],
                          "all_checks_pass": report["acceptance"]["all_checks_pass"],
                          "baseline_gap_dex": report["matrix"]["baseline_absolute_gap_dex"],
                          "best_gap_cell": report["matrix"]["best_gap_cell"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
