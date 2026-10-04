#!/usr/bin/env python3
"""Run paired Sentaurus/Vela model ablations for the SimpleMOS M8 residual."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
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
import compare_simplemos_m4_controlled_matrix as compare_m4  # noqa: E402
import run_simplemos_m3_workflow as m3  # noqa: E402
import run_simplemos_m4_controlled_matrix as m4  # noqa: E402
import run_simplemos_m8_original_matrix as m8  # noqa: E402


DEFAULT_CONTRACT = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_m8a_model_ablation_contract_v1.json"
)
DEFAULT_OUTPUT = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m8a_model_ablation"
)
DEFAULT_TDR = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m8_upstream_tdrs/n23_fps.tdr"
)
DEFAULT_MATERIALS = (
    REPO / "reference_tcad/transportmodels_sentaurus2022/vela"
    / "materials_sentaurus2022.json"
)
DEFAULT_VALIDATION_CONTRACT = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_sdevice_validation_contract_v1.json"
)
DEFAULT_REMOTE_ROOT = (
    "~/sentaurus_runs/vela_oracle/simplemos_m8a_model_ablation_20260828"
)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO).as_posix()
    except ValueError:
        return str(resolved)


def executable(name: str) -> str:
    if os.name == "nt":
        candidate = (Path(os.environ.get("SystemRoot", r"C:\Windows"))
                     / "System32" / "OpenSSH" / f"{name}.exe")
        if candidate.is_file():
            return str(candidate)
    return shutil.which(name) or name


def run(argv: Sequence[str], *, capture: bool = False) -> str:
    completed = subprocess.run(
        list(argv), cwd=REPO, check=True, text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
    )
    return completed.stdout or ""


def voltage_tag(value: float) -> str:
    return format(value, ".12g").replace("-", "m").replace(".", "p")


def validate_contract(contract: dict[str, Any]) -> None:
    expected_ids = [
        "full", "no_srh", "plain_srh", "no_bgn", "no_enormal",
        "no_hfs", "phumob_only", "constant_mu",
    ]
    variants = contract["variants"]
    if [item["id"] for item in variants] != expected_ids:
        raise ValueError("M8-A variants must use the frozen ablation order")
    if contract["device"]["id"] != "n23":
        raise ValueError("M8-A first round must use n23")
    if [float(value) for value in contract["bias_matrix"][
            "drain_voltages_V"]] != [0.05]:
        raise ValueError("M8-A first round must use Vd=0.05 V")
    lattice = [float(value) for value in
               contract["bias_matrix"]["gate_lattice"]["values_V"]]
    if len(lattice) != 51:
        raise ValueError("M8-A gate lattice must contain 51 points")
    for index, value in enumerate(lattice):
        if not math.isclose(value, index * 0.05, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError("M8-A gate lattice must be exactly 0:0.05:2.5 V")
    if contract["comparison"]["interpolation"] != "forbidden":
        raise ValueError("M8-A forbids interpolation")


def sentaurus_physics(variant: dict[str, Any]) -> str:
    sections: list[str] = []
    if variant["old_slotboom"]:
        sections.append("Physics { EffectiveIntrinsicDensity(OldSlotboom) }")
    silicon: list[str] = []
    models = list(variant["sentaurus_mobility"])
    if models:
        silicon.append("   Mobility(" + " ".join(models) + ")")
    recombination = variant["recombination"]
    if recombination == "srh_doping_dependence":
        silicon.append("   Recombination(SRH(DopingDependence))")
    elif recombination == "srh_plain":
        silicon.append("   Recombination(SRH)")
    elif recombination != "none":
        raise ValueError(f"unsupported recombination {recombination!r}")
    if silicon:
        sections.append(
            'Physics(Material="Silicon") {\n' + "\n".join(silicon) + "\n}")
    return "\n".join(sections)


def sentaurus_deck(case: str, vd: float, variant: dict[str, Any],
                    intervals: int) -> str:
    physics = sentaurus_physics(variant)
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

{physics}

Math {{ Extrapolate Iterations=20 ExitOnFailure }}

Solve {{
   Coupled(Iterations=100) {{ Poisson }}
   Coupled {{ Poisson Electron Hole }}
   Quasistationary(
      InitialStep=0.1 Increment=1.5 MinStep=1e-5 MaxStep=1
      Goal {{ Name="drain" Voltage={vd:.17g} }}
   ) {{ Coupled {{ Poisson Electron Hole }} }}
   NewCurrentPrefix="IdVg_"
   Quasistationary(
      DoZero InitialStep=0.01 Increment=1.5 MinStep=1e-5 MaxStep=0.05
      Goal {{ Name="gate" Voltage=2.5 }}
   ) {{
      Coupled {{ Poisson Electron Hole }}
      CurrentPlot(Time=(Range=(0 1) Intervals={intervals}))
   }}
}}
'''


def prepare_sentaurus_bundle(contract: dict[str, Any], contract_path: Path,
                              tdr: Path, output_dir: Path) -> dict[str, Any]:
    bundle = output_dir / "sentaurus_bundle"
    bundle.mkdir(parents=True, exist_ok=True)
    intervals = int(contract["bias_matrix"]["gate_lattice"]["point_count"]) - 1
    cases: list[dict[str, Any]] = []
    for variant in contract["variants"]:
        variant_id = str(variant["id"])
        case_dir = bundle / variant_id
        case_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(tdr, case_dir / "input_fps.tdr")
        for vd in contract["bias_matrix"]["drain_voltages_V"]:
            case = f"{variant_id}_vd_{voltage_tag(float(vd))}"
            deck = case_dir / f"{case}_des.cmd"
            deck.write_text(
                sentaurus_deck(case, float(vd), variant, intervals),
                encoding="utf-8", newline="\n")
            cases.append({
                "case": case,
                "variant": variant_id,
                "drain_voltage_V": float(vd),
                "deck": portable_path(deck),
                "deck_sha256": sha256(deck),
                "expected_plot": f"IdVg_{case}_des.plt",
                "parameter_snapshot": f"{case}_models.par",
            })
    manifest = {
        "schema": "vela.simplemos.sdevice.m8a_sentaurus_matrix.v1",
        "status": "prepared",
        "contract": portable_path(contract_path),
        "contract_sha256": sha256(contract_path),
        "tdr": portable_path(tdr),
        "tdr_sha256": sha256(tdr),
        "cases": cases,
    }
    write_json(output_dir / "sentaurus_matrix_manifest.json", manifest)
    return manifest


def run_sentaurus_vm(manifest: dict[str, Any], output_dir: Path,
                     ssh_target: str, ssh_bin: str, scp_bin: str,
                     remote_root: str, jobs: int) -> str:
    if jobs < 1:
        raise ValueError("Sentaurus jobs must be at least one")
    banner = run(
        [ssh_bin, ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'"],
        capture=True).strip()
    if "T-2022.03-SP2" not in banner:
        raise RuntimeError(f"unexpected Sentaurus release:\n{banner}")
    run([ssh_bin, ssh_target,
         f"test ! -e {remote_root} && mkdir -p {remote_root}"])
    run([scp_bin, "-r", str(output_dir / "sentaurus_bundle"),
         f"{ssh_target}:{remote_root}/"])

    def execute_case(index_case: tuple[int, dict[str, Any]]) -> None:
        index, case = index_case
        name = str(case["case"])
        variant = str(case["variant"])
        print(f"[{index}/{len(manifest['cases'])}] Sentaurus {name}", flush=True)
        command = (
            f"set -eu; cd {remote_root}/sentaurus_bundle/{variant}; "
            f"test ! -e models.par; "
            f"sdevice -P {name}_des.cmd > {name}.parameters.log 2>&1; "
            f"mv models.par {name}_models.par; "
            f"sdevice {name}_des.cmd > {name}.console.log 2>&1"
        )
        run([ssh_bin, ssh_target, command])

    indexed = list(enumerate(manifest["cases"], start=1))
    if jobs == 1:
        for item in indexed:
            execute_case(item)
    else:
        with ThreadPoolExecutor(max_workers=min(jobs, len(indexed))) as executor:
            list(executor.map(execute_case, indexed))

    archive_name = "simplemos_m8a_results.tgz"
    run([ssh_bin, ssh_target,
         f"cd {remote_root} && tar -czf {archive_name} sentaurus_bundle"])
    raw = output_dir / "sentaurus_raw"
    raw.mkdir(parents=True, exist_ok=True)
    archive = raw / archive_name
    run([scp_bin, f"{ssh_target}:{remote_root}/{archive_name}", str(archive)])
    with tarfile.open(archive, "r:gz") as stream:
        stream.extractall(raw, filter="data")
    (output_dir / "sentaurus_banner.txt").write_text(
        banner + "\n", encoding="utf-8", newline="\n")
    return banner


def extract_references(contract: dict[str, Any], manifest: dict[str, Any],
                       output_dir: Path, banner: str) -> dict[str, Any]:
    lattice = [float(value) for value in
               contract["bias_matrix"]["gate_lattice"]["values_V"]]
    tolerance = float(contract["comparison"]["exact_bias_tolerance_V"])
    raw = output_dir / "sentaurus_raw/sentaurus_bundle"
    reference_dir = output_dir / "sentaurus_reference"
    reference_dir.mkdir(parents=True, exist_ok=True)
    artifacts: list[dict[str, Any]] = []
    for case in manifest["cases"]:
        source_dir = raw / case["variant"]
        rows = m4.exact_plt_curve(
            source_dir / case["expected_plot"], lattice, tolerance)
        output = reference_dir / f"{case['case']}_reference.csv"
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        parameter_snapshot = source_dir / case["parameter_snapshot"]
        artifacts.append({
            "case": case["case"],
            "variant": case["variant"],
            "path": portable_path(output),
            "sha256": sha256(output),
            "point_count": len(rows),
            "parameter_snapshot": portable_path(parameter_snapshot),
            "parameter_snapshot_sha256": sha256(parameter_snapshot),
        })
    result = {
        "schema": "vela.simplemos.sdevice.m8a_reference_curves.v1",
        "status": "qualified",
        "sentaurus_banner": banner,
        "contract_sha256": manifest["contract_sha256"],
        "tdr_sha256": manifest["tdr_sha256"],
        "interpolation": "forbidden",
        "artifacts": artifacts,
    }
    write_json(reference_dir / "reference_manifest.json", result)
    return result


def mobility_config(variant: dict[str, Any]) -> dict[str, Any]:
    config: dict[str, Any] = {
        "model": variant["vela_mobility"],
        "doping_concentration_basis": "total_impurity",
    }
    if "field" in str(variant["vela_mobility"]):
        config.update({
            "high_field_driving_force": "quasi_fermi_gradient",
            "high_field_gradient_discretization": "transport_cell_vector",
        })
    return config


def apply_variant_physics(config: dict[str, Any], variant: dict[str, Any],
                          contract: dict[str, Any]) -> None:
    solver = config["solver"]
    solver["mobility"] = mobility_config(variant)
    solver["bandgap_narrowing"] = {
        "model": "old_slotboom" if variant["old_slotboom"] else "none",
        "fermi_statistics_correction": False,
    }
    defaults = contract["srh_defaults"]
    recombination = variant["recombination"]
    if recombination == "none":
        solver["recombination"] = ["none"]
        solver["srh_doping_dependence"]["enabled"] = False
    elif recombination == "srh_plain":
        solver["recombination"] = ["srh"]
        solver["srh_doping_dependence"]["enabled"] = False
        solver["taun"] = float(defaults["electron_lifetime_s"])
        solver["taup"] = float(defaults["hole_lifetime_s"])
    elif recombination == "srh_doping_dependence":
        solver["recombination"] = ["srh"]
        srh = solver["srh_doping_dependence"]
        srh["enabled"] = True
        for carrier, lifetime_key in (
                ("electron", "electron_lifetime_s"),
                ("hole", "hole_lifetime_s")):
            srh[carrier].update({
                "tau_min_s": 0.0,
                "tau_max_s": float(defaults[lifetime_key]),
                "reference_doping_m3": float(
                    defaults["reference_doping_cm3"]),
                "gamma": float(defaults["gamma"]),
            })
    else:
        raise ValueError(f"unsupported recombination {recombination!r}")


def prepare_vela_workflows(contract: dict[str, Any], tdr: Path,
                           output_dir: Path, materials: Path,
                           importer: Path) -> dict[str, Any]:
    neutral = output_dir / "neutral/n23"
    m8.import_tdr(tdr, neutral, importer)
    vela = output_dir / "vela"
    vela.mkdir(parents=True, exist_ok=True)
    mesh = vela / "mesh.json"
    write_json(mesh, m4.neutral_mesh(neutral))
    lattice = [float(value) for value in
               contract["bias_matrix"]["gate_lattice"]["values_V"]]
    stages: list[dict[str, Any]] = []
    for variant in contract["variants"]:
        variant_id = str(variant["id"])
        base = vela / variant_id / "base.json"
        config = m4.base_config(
            mesh, neutral / "doping.csv", materials, mobility_config(variant))
        apply_variant_physics(config, variant, contract)
        write_json(base, config)
        workflow_dir = vela / variant_id / "workflow"
        workflow = m3.materialize(
            base, DEFAULT_VALIDATION_CONTRACT, workflow_dir,
            [float(value) for value in contract["bias_matrix"][
                "drain_voltages_V"]])
        for stage in workflow["stages"]:
            config_path = Path(stage["config"])
            stage_config = read_json(config_path)
            if stage["phase"] == "equilibrium":
                stage_config["solver"]["quasi_fermi_update_limit_V"] = 0.1
            elif stage["phase"] == "drain_ramp":
                stage_config["solver"].update({
                    "line_search_mode": "block_filter",
                    "residual_filter_gamma": 1e-4,
                    "residual_filter_envelope_factor": 2.0,
                })
                stage_config["sweep"]["step"] = 0.005
                stage_config["sweep"]["max_step"] = 0.01
            else:
                stage_config["sweep"]["bias_points"] = lattice
                stage_config["sweep"]["step"] = 0.05
                stage_config["sweep"]["min_step"] = min(
                    float(stage_config["sweep"]["min_step"]), 0.05)
                stage_config["sweep"]["max_step"] = max(
                    float(stage_config["sweep"]["max_step"]), 0.05)
                stage["exact_gate_lattice_V"] = lattice
            stage_config["simplemos_m8a"] = {
                "variant": variant_id,
                "paired_sentaurus_configuration": True,
            }
            write_json(config_path, stage_config)
            stage["config_sha256"] = sha256(config_path)
        workflow["simplemos_m8a"] = {"variant": variant_id, "device": "n23"}
        write_json(workflow_dir / "workflow_manifest.json", workflow)
        stages.extend(workflow["stages"])
    result = {
        "schema": "vela.simplemos.sdevice.m8a_vela_matrix.v1",
        "status": "materialized",
        "tdr_sha256": sha256(tdr),
        "mesh": portable_path(mesh),
        "mesh_sha256": sha256(mesh),
        "stages": stages,
    }
    write_json(vela / "vela_matrix_manifest.json", result)
    return result


def execute_vela(contract: dict[str, Any], output_dir: Path,
                 runner: Path, jobs: int) -> dict[str, Any]:
    if jobs < 1:
        raise ValueError("Vela jobs must be at least one")

    def execute_variant(index_variant: tuple[int, dict[str, Any]]) -> dict[str, Any]:
        index, variant = index_variant
        variant_id = str(variant["id"])
        print(f"[{index}/{len(contract['variants'])}] Vela {variant_id}", flush=True)
        workflow_dir = output_dir / "vela" / variant_id / "workflow"
        manifest = workflow_dir / "workflow_manifest.json"
        workflow = m3.execute(read_json(manifest), runner.resolve(), workflow_dir)
        return {
            "variant": variant_id,
            "status": workflow["status"],
            "manifest": portable_path(manifest),
            "manifest_sha256": sha256(manifest),
        }

    indexed = list(enumerate(contract["variants"], start=1))
    if jobs == 1:
        results = [execute_variant(item) for item in indexed]
    else:
        with ThreadPoolExecutor(max_workers=min(jobs, len(indexed))) as executor:
            results = list(executor.map(execute_variant, indexed))
    report = {
        "schema": "vela.simplemos.sdevice.m8a_vela_execution.v1",
        "status": "accepted" if all(
            item["status"] == "accepted" for item in results) else "fail",
        "runner": portable_path(runner),
        "runner_sha256": sha256(runner),
        "variants": results,
    }
    write_json(output_dir / "vela/vela_execution_report.json", report)
    return report


def diagnostic_metrics(rows: list[dict[str, Any]],
                       contract: dict[str, Any]) -> dict[str, Any]:
    windows = contract["diagnostic_windows"]

    def at_bias(target: float) -> dict[str, Any]:
        matches = [row for row in rows if math.isclose(
            float(row["gate_voltage_V"]), target, rel_tol=0.0, abs_tol=1e-12)]
        if len(matches) != 1:
            raise ValueError(f"expected one comparison row at Vg={target:g} V")
        row = matches[0]
        return {
            "gate_voltage_V": target,
            "sentaurus_current_A_per_um": row["sentaurus_current_A_per_um"],
            "vela_current_A_per_um": row["vela_current_A_per_um"],
            "absolute_log10_ratio": row["absolute_log10_ratio"],
            "relative_error": row["relative_error"],
        }

    weak_rows = [row for row in rows if
                 float(windows["weak_inversion_min_V"]) <=
                 float(row["gate_voltage_V"]) <=
                 float(windows["weak_inversion_max_V"])]
    weak_max = max(weak_rows, key=lambda row: float(row["absolute_log10_ratio"]))
    return {
        "deep_off": [at_bias(float(value)) for value in
                     windows["deep_off_gate_voltages_V"]],
        "weak_inversion_max": {
            "gate_voltage_V": weak_max["gate_voltage_V"],
            "absolute_log10_ratio": weak_max["absolute_log10_ratio"],
            "relative_error": weak_max["relative_error"],
        },
        "strong_inversion": at_bias(float(
            windows["strong_inversion_gate_voltage_V"])),
    }


def compare_results(contract: dict[str, Any], output_dir: Path) -> dict[str, Any]:
    comparisons = output_dir / "comparisons"
    comparisons.mkdir(parents=True, exist_ok=True)
    cases: list[dict[str, Any]] = []
    for variant in contract["variants"]:
        variant_id = str(variant["id"])
        for vd in contract["bias_matrix"]["drain_voltages_V"]:
            tag = voltage_tag(float(vd))
            case = f"{variant_id}_vd_{tag}"
            reference = output_dir / "sentaurus_reference" / f"{case}_reference.csv"
            candidate = (output_dir / "vela" / variant_id / "workflow"
                         / f"vd_{tag}" / "20_gate_sweep.csv")
            result = compare_m4.compare_case(reference, candidate, contract)
            rows = result.pop("rows")
            result["diagnostic_metrics"] = diagnostic_metrics(rows, contract)
            comparison_csv = comparisons / f"{case}_comparison.csv"
            compare_m4.write_case_csv(comparison_csv, rows)
            result.update({
                "case": case,
                "variant": variant_id,
                "drain_voltage_V": float(vd),
                "comparison_csv": portable_path(comparison_csv),
                "comparison_csv_sha256": sha256(comparison_csv),
            })
            cases.append(result)
    baseline = next(item for item in cases if item["variant"] == "full")
    for item in cases:
        item["effect_vs_full"] = {
            "maximum_log10_ratio_change": (
                item["maximum_absolute_log10_ratio_above_floor"]
                - baseline["maximum_absolute_log10_ratio_above_floor"]),
            "maximum_relative_error_change": (
                item["maximum_relative_error_above_floor"]
                - baseline["maximum_relative_error_above_floor"]),
            "deep_off_vg_0p05_log10_ratio_change": (
                item["diagnostic_metrics"]["deep_off"][1]["absolute_log10_ratio"]
                - baseline["diagnostic_metrics"]["deep_off"][1][
                    "absolute_log10_ratio"]),
            "weak_inversion_max_log10_ratio_change": (
                item["diagnostic_metrics"]["weak_inversion_max"][
                    "absolute_log10_ratio"]
                - baseline["diagnostic_metrics"]["weak_inversion_max"][
                    "absolute_log10_ratio"]),
        }
    report = {
        "schema": "vela.simplemos.sdevice.m8a_comparison.v1",
        "status": "complete" if all(
            item["status"] in {"pass", "fail"} for item in cases)
        else "qualification_failed",
        "comparison_basis": "paired identical physics configuration per variant",
        "interpolation": "forbidden",
        "cases": cases,
    }
    write_json(comparisons / "comparison_report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--tdr", type=Path, default=DEFAULT_TDR)
    parser.add_argument("--materials", type=Path, default=DEFAULT_MATERIALS)
    parser.add_argument("--importer", type=Path,
                        default=REPO / "build-release/sentaurus_import.exe")
    parser.add_argument("--runner", type=Path,
                        default=REPO / "build-release/vela_example_runner.exe")
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("--reuse-prepared", action="store_true")
    parser.add_argument("--live-sentaurus", action="store_true")
    parser.add_argument("--extract-existing", action="store_true")
    parser.add_argument("--execute-vela", action="store_true")
    parser.add_argument("--compare", action="store_true")
    parser.add_argument("--vela-jobs", type=int, default=1)
    parser.add_argument("--sentaurus-jobs", type=int, default=1)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--scp-bin", default=executable("scp"))
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    args = parser.parse_args()

    contract_path = args.contract.resolve()
    contract = read_json(contract_path)
    validate_contract(contract)
    output_dir = args.output_dir.resolve()
    if args.clean and output_dir.exists():
        shutil.rmtree(output_dir)
    tdr = args.tdr.resolve()
    if sha256(tdr) != contract["device"]["tdr_sha256"]:
        raise ValueError("n23 TDR hash does not match the frozen M8-A contract")

    if args.reuse_prepared:
        sentaurus_manifest = read_json(
            output_dir / "sentaurus_matrix_manifest.json")
    else:
        sentaurus_manifest = prepare_sentaurus_bundle(
            contract, contract_path, tdr, output_dir)
        prepare_vela_workflows(
            contract, tdr, output_dir, args.materials.resolve(),
            args.importer.resolve())

    vela_report = None
    if args.execute_vela:
        vela_report = execute_vela(
            contract, output_dir, args.runner.resolve(), args.vela_jobs)

    if args.live_sentaurus:
        banner = run_sentaurus_vm(
            sentaurus_manifest, output_dir, args.ssh_target, args.ssh_bin,
            args.scp_bin, args.remote_root, args.sentaurus_jobs)
        extract_references(contract, sentaurus_manifest, output_dir, banner)
    elif args.extract_existing:
        banner = (output_dir / "sentaurus_banner.txt").read_text(
            encoding="utf-8").strip()
        extract_references(contract, sentaurus_manifest, output_dir, banner)

    comparison = compare_results(contract, output_dir) if args.compare else None
    summary = {
        "variants": len(contract["variants"]),
        "sentaurus_cases": len(sentaurus_manifest["cases"]),
        "vela_status": (vela_report or {}).get("status"),
        "comparison_status": (comparison or {}).get("status"),
    }
    print(json.dumps(summary))
    return 1 if vela_report is not None and vela_report["status"] != "accepted" else 0


if __name__ == "__main__":
    raise SystemExit(main())
