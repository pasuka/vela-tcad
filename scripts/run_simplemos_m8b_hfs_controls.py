#!/usr/bin/env python3
"""Run Sentaurus-only HFS control diagnostics against the frozen M8-A baseline."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import json
import math
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import tarfile
from typing import Any, Sequence


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import compare_simplemos_m4_controlled_matrix as compare_m4  # noqa: E402
import run_simplemos_m4_controlled_matrix as m4  # noqa: E402
import run_simplemos_m8_original_matrix as m8  # noqa: E402
import run_simplemos_m8a_model_ablation as m8a  # noqa: E402
import run_simplemos_m8a_confirmation as confirmation  # noqa: E402

# Re-export the shared digest helper so regression checks can verify the
# complete frozen evidence chain through this orchestration module.
sha256 = m8a.sha256


DEFAULT_CONTRACT = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_m8b_hfs_controls_contract_v1.json"
)
DEFAULT_OUTPUT = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m8b_hfs_controls"
)
DEFAULT_BASELINE = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m8a_confirmation"
)
DEFAULT_RAW = REPO / "build-release/m8b_hfs_raw"
DEFAULT_TDRS = confirmation.DEFAULT_TDRS
DEFAULT_REMOTE_ROOT = (
    "~/sentaurus_runs/vela_oracle/simplemos_m8b_hfs_controls_20260828_v2"
)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8",
                    newline="\n")


def run(argv: Sequence[str], *, capture: bool = False) -> str:
    completed = subprocess.run(
        list(argv), cwd=REPO, check=True, text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None)
    return completed.stdout or ""


def extract_archive(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:gz") as stream:
        stream.extractall(destination, filter="data")


def validate_contract(contract: dict[str, Any]) -> None:
    if [item["id"] for item in contract["devices"]] != ["n17", "n21"]:
        raise ValueError("M8-B devices must be n17 and n21")
    expected = [
        "explicit_gradqf", "eparallel", "qf_at_contacts",
        "no_parallel_boundary", "refdens_efield_1e8",
    ]
    if [item["id"] for item in contract["variants"]] != expected:
        raise ValueError("M8-B HFS controls must use the frozen order")
    if [float(value) for value in contract["bias_matrix"][
            "drain_voltages_V"]] != [0.05, 1.0]:
        raise ValueError("M8-B drain biases must be 0.05 and 1 V")
    lattice = [float(value) for value in contract["bias_matrix"][
        "gate_lattice"]["values_V"]]
    if len(lattice) != 51 or any(not math.isclose(
            value, index * 0.05, rel_tol=0.0, abs_tol=1e-12)
            for index, value in enumerate(lattice)):
        raise ValueError("M8-B gate lattice must be 0:0.05:2.5 V")
    if contract["comparison"]["interpolation"] != "forbidden":
        raise ValueError("M8-B forbids bias interpolation")


def validate_inputs(contract: dict[str, Any], tdrs: dict[str, Path],
                    baseline: Path) -> None:
    expected = {item["id"]: item["tdr_sha256"]
                for item in contract["devices"]}
    for device, digest in expected.items():
        path = tdrs[device]
        if m8a.sha256(path) != digest:
            raise ValueError(f"{device} TDR hash does not match M8-B contract")
        for vd in contract["bias_matrix"]["drain_voltages_V"]:
            tag = m8a.voltage_tag(float(vd))
            reference = (baseline / "sentaurus_reference"
                         / f"{device}_full_vd_{tag}_reference.csv")
            candidate = (baseline / "vela" / device / "full" / "workflow"
                         / f"vd_{tag}" / "20_gate_sweep.csv")
            if not reference.is_file() or not candidate.is_file():
                raise FileNotFoundError(
                    f"missing frozen M8-A baseline for {device}, Vd={vd:g}")


def math_block(options: list[str]) -> str:
    lines = ["Math {", "   Extrapolate Iterations=20 ExitOnFailure"]
    lines.extend(f"   {option}" for option in options)
    lines.append("}")
    return "\n".join(lines)


def sentaurus_deck(case: str, vd: float, variant: dict[str, Any],
                   intervals: int,
                   numerical_retry: dict[str, Any] | None = None) -> str:
    full = {
        "old_slotboom": True,
        "sentaurus_mobility": ["PhuMob", "HighFieldSaturation", "Enormal"],
        "recombination": "srh_doping_dependence",
    }
    text = m8a.sentaurus_deck(case, vd, full, intervals)
    text = text.replace(
        "Mobility(PhuMob HighFieldSaturation Enormal)",
        f"Mobility(PhuMob {variant['hfs_option']} Enormal)")
    text = text.replace(
        "Math { Extrapolate Iterations=20 ExitOnFailure }",
        math_block(list(variant["math_options"])))
    if numerical_retry:
        iterations = int(numerical_retry["coupled_iterations"])
        min_step = float(numerical_retry["gate_min_step"])
        text = text.replace(
            "   Extrapolate Iterations=20 ExitOnFailure",
            f"   Extrapolate Iterations={iterations} ExitOnFailure")
        text = text.replace(
            "      DoZero InitialStep=0.01 Increment=1.5 MinStep=1e-5 MaxStep=0.05",
            "      DoZero InitialStep=0.01 Increment=1.5 "
            f"MinStep={min_step:g} MaxStep=0.05")
        if "line_search_damping" in numerical_retry:
            damping = float(numerical_retry["line_search_damping"])
            target = "      Coupled { Poisson Electron Hole }\n      CurrentPlot"
            replacement = (
                f"      Coupled(Iterations={iterations} "
                f"LineSearchDamping={damping:g}) "
                "{ Poisson Electron Hole }\n      CurrentPlot")
            if target not in text:
                raise ValueError("gate-sweep Coupled block was not found")
            text = text.replace(target, replacement, 1)
    return text


def prepare_sentaurus(contract: dict[str, Any], contract_path: Path,
                       tdrs: dict[str, Path], output: Path) -> dict[str, Any]:
    bundle = output / "sentaurus_bundle"
    bundle.mkdir(parents=True, exist_ok=True)
    intervals = len(contract["bias_matrix"]["gate_lattice"]["values_V"]) - 1
    retries = {item["case"]: item for item in contract.get(
        "numerical_retry_policy", {}).get("cases", [])}
    cases: list[dict[str, Any]] = []
    for device in contract["devices"]:
        device_id = str(device["id"])
        device_dir = bundle / device_id
        device_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(tdrs[device_id], device_dir / "input_fps.tdr")
        for variant in contract["variants"]:
            variant_id = str(variant["id"])
            variant_dir = device_dir / variant_id
            variant_dir.mkdir(parents=True, exist_ok=True)
            for vd in contract["bias_matrix"]["drain_voltages_V"]:
                tag = m8a.voltage_tag(float(vd))
                case = f"{device_id}_{variant_id}_vd_{tag}"
                # sdevice -P always writes models.par in the current working
                # directory.  Keep every bias in a separate directory so
                # parallel parameter snapshots cannot race each other.
                case_dir = variant_dir / case
                case_dir.mkdir(parents=True, exist_ok=True)
                deck = case_dir / f"{case}_des.cmd"
                text = sentaurus_deck(
                    case, float(vd), variant, intervals,
                    retries.get(case)).replace(
                        'Grid="input_fps.tdr"', 'Grid="../../input_fps.tdr"')
                deck.write_text(text, encoding="utf-8", newline="\n")
                cases.append({
                    "case": case,
                    "device": device_id,
                    "variant": variant_id,
                    "drain_voltage_V": float(vd),
                    "deck": m8a.portable_path(deck),
                    "deck_sha256": m8a.sha256(deck),
                    "expected_plot": f"IdVg_{case}_des.plt",
                    "parameter_snapshot": f"{case}_models.par",
                })
    manifest = {
        "schema": "vela.simplemos.sdevice.m8b_hfs_controls_sentaurus.v1",
        "status": "prepared",
        "contract": m8a.portable_path(contract_path),
        "contract_sha256": m8a.sha256(contract_path),
        "tdr_sha256": {key: m8a.sha256(value) for key, value in tdrs.items()},
        "cases": cases,
    }
    write_json(output / "sentaurus_matrix_manifest.json", manifest)
    return manifest


def run_sentaurus(manifest: dict[str, Any], output: Path, ssh_target: str,
                  ssh_bin: str, scp_bin: str, remote_root: str,
                  jobs: int, raw_extract: Path) -> str:
    if jobs < 1:
        raise ValueError("Sentaurus jobs must be at least one")
    banner = run(
        [ssh_bin, "-n", ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'"],
        capture=True).strip()
    if "T-2022.03-SP2" not in banner:
        raise RuntimeError(f"unexpected Sentaurus release:\n{banner}")
    run([ssh_bin, "-n", ssh_target, f"mkdir -p {remote_root}"])
    run([scp_bin, "-r", str(output / "sentaurus_bundle"),
         f"{ssh_target}:{remote_root}/"])

    def execute(index_case: tuple[int, dict[str, Any]]) -> None:
        index, case = index_case
        name = str(case["case"])
        remote_dir = (f"{remote_root}/sentaurus_bundle/{case['device']}/"
                      f"{case['variant']}/{name}")
        print(f"[{index}/{len(manifest['cases'])}] Sentaurus {name}",
              flush=True)
        command = (
            f"set -eu; cd {remote_dir}; "
            f"if test -s {name}.console.log && "
            f"grep -q 'Good Bye' {name}.console.log && "
            f"! grep -q 'Exit due to failure' {name}.console.log; "
            f"then exit 0; fi; "
            f"sdevice -P {name}_des.cmd > {name}.parameters.log 2>&1; "
            f"mv -f models.par {name}_models.par; "
            f"sdevice {name}_des.cmd > {name}.console.log 2>&1"
        )
        run([ssh_bin, "-n", ssh_target, command])

    indexed = list(enumerate(manifest["cases"], start=1))
    with ThreadPoolExecutor(max_workers=min(jobs, len(indexed))) as executor:
        list(executor.map(execute, indexed))

    archive_name = "simplemos_m8b_hfs_controls_results.tgz"
    run([ssh_bin, "-n", ssh_target,
         f"cd {remote_root} && tar -czf {archive_name} sentaurus_bundle"])
    raw = output / "sentaurus_raw"
    raw.mkdir(parents=True, exist_ok=True)
    archive = raw / archive_name
    run([scp_bin, f"{ssh_target}:{remote_root}/{archive_name}", str(archive)])
    (output / "sentaurus_banner.txt").write_text(
        banner + "\n", encoding="utf-8", newline="\n")
    extract_archive(archive, raw_extract)
    return banner


def extract_references(contract: dict[str, Any], manifest: dict[str, Any],
                       output: Path, raw_extract: Path,
                       banner: str) -> dict[str, Any]:
    lattice = [float(value) for value in contract["bias_matrix"][
        "gate_lattice"]["values_V"]]
    tolerance = float(contract["comparison"]["exact_bias_tolerance_V"])
    destination = output / "sentaurus_reference"
    destination.mkdir(parents=True, exist_ok=True)
    artifacts = []
    for case in manifest["cases"]:
        source_dir = (raw_extract / "sentaurus_bundle"
                      / case["device"] / case["variant"] / case["case"])
        rows = m4.exact_plt_curve(
            source_dir / case["expected_plot"], lattice, tolerance)
        path = destination / f"{case['case']}_reference.csv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]),
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        parameters = source_dir / case["parameter_snapshot"]
        artifacts.append({
            "case": case["case"],
            "device": case["device"],
            "variant": case["variant"],
            "point_count": len(rows),
            "path": m8a.portable_path(path),
            "sha256": m8a.sha256(path),
            "parameter_snapshot": m8a.portable_path(parameters),
            "parameter_snapshot_sha256": m8a.sha256(parameters),
        })
    report = {
        "schema": "vela.simplemos.sdevice.m8b_hfs_controls_reference.v1",
        "status": "qualified",
        "sentaurus_banner": banner,
        "contract_sha256": manifest["contract_sha256"],
        "tdr_sha256": manifest["tdr_sha256"],
        "interpolation": "forbidden",
        "artifacts": artifacts,
    }
    write_json(destination / "reference_manifest.json", report)
    return report


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _response_rows(baseline_reference: Path, controlled_reference: Path,
                   comparison_csv: Path) -> list[dict[str, Any]]:
    baseline = _read_csv(baseline_reference)
    controlled = _read_csv(controlled_reference)
    comparison = _read_csv(comparison_csv)
    if not (len(baseline) == len(controlled) == len(comparison)):
        raise ValueError("response inputs must use the same exact bias lattice")
    rows = []
    for base, control, result in zip(baseline, controlled, comparison):
        vb = float(base["gate_voltage_V"])
        vc = float(control["gate_voltage_V"])
        vr = float(result["gate_voltage_V"])
        if not (math.isclose(vb, vc, abs_tol=1e-12) and
                math.isclose(vb, vr, abs_tol=1e-12)):
            raise ValueError("response inputs are not bias aligned")
        base_current = abs(float(base["drain_total_current_A_per_um"]))
        controlled_current = abs(float(
            control["drain_total_current_A_per_um"]))
        response = math.log10(controlled_current / base_current)
        rows.append({
            "gate_voltage_V": vb,
            "baseline_sentaurus_current_A_per_um": base_current,
            "controlled_sentaurus_current_A_per_um": controlled_current,
            "sentaurus_log10_response": response,
            "absolute_sentaurus_log10_response": abs(response),
            "vela_current_A_per_um": float(result["vela_current_A_per_um"]),
            "controlled_absolute_log10_error": float(
                result["absolute_log10_ratio"]),
        })
    return rows


def _write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]),
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def compare(contract: dict[str, Any], output: Path,
            baseline: Path) -> dict[str, Any]:
    destination = output / "comparisons"
    destination.mkdir(parents=True, exist_ok=True)
    baseline_report = read_json(baseline / "comparisons/comparison_report.json")
    baseline_cases = {
        (item["device"], float(item["drain_voltage_V"])): item
        for item in baseline_report["cases"] if item["variant"] == "full"
    }
    cases = []
    for device in contract["devices"]:
        device_id = str(device["id"])
        for variant in contract["variants"]:
            variant_id = str(variant["id"])
            for vd in contract["bias_matrix"]["drain_voltages_V"]:
                vd_value = float(vd)
                tag = m8a.voltage_tag(vd_value)
                case = f"{device_id}_{variant_id}_vd_{tag}"
                controlled = (output / "sentaurus_reference"
                              / f"{case}_reference.csv")
                vela = (baseline / "vela" / device_id / "full" / "workflow"
                        / f"vd_{tag}" / "20_gate_sweep.csv")
                baseline_reference = (baseline / "sentaurus_reference"
                                      / f"{device_id}_full_vd_{tag}_reference.csv")
                result = compare_m4.compare_case(controlled, vela, contract)
                comparison_rows = result.pop("rows")
                comparison_csv = destination / f"{case}_comparison.csv"
                compare_m4.write_case_csv(comparison_csv, comparison_rows)
                responses = _response_rows(
                    baseline_reference, controlled, comparison_csv)
                response_csv = destination / f"{case}_response.csv"
                _write_rows(response_csv, responses)
                baseline_case = baseline_cases[(device_id, vd_value)]
                response_values = [float(row["sentaurus_log10_response"])
                                   for row in responses]
                abs_errors = [float(row["controlled_absolute_log10_error"])
                              for row in responses]
                result.update({
                    "case": case,
                    "device": device_id,
                    "variant": variant_id,
                    "dimension": variant["dimension"],
                    "drain_voltage_V": vd_value,
                    "comparison_csv": m8a.portable_path(comparison_csv),
                    "comparison_csv_sha256": m8a.sha256(comparison_csv),
                    "response_csv": m8a.portable_path(response_csv),
                    "response_csv_sha256": m8a.sha256(response_csv),
                    "maximum_sentaurus_response_dex": max(
                        abs(value) for value in response_values),
                    "median_sentaurus_response_dex": statistics.median(
                        abs(value) for value in response_values),
                    "mean_absolute_log10_error": statistics.fmean(abs_errors),
                    "effect_vs_default": {
                        "maximum_log10_error_change": (
                            result["maximum_absolute_log10_ratio_above_floor"]
                            - baseline_case[
                                "maximum_absolute_log10_ratio_above_floor"]),
                        "median_log10_error_change": (
                            result["median_absolute_log10_ratio_above_floor"]
                            - baseline_case[
                                "median_absolute_log10_ratio_above_floor"]),
                    },
                })
                cases.append(result)

    policy = contract["attribution_policy"]
    summaries = []
    for variant in contract["variants"]:
        variant_cases = [item for item in cases
                         if item["variant"] == variant["id"]]
        max_changes = [item["effect_vs_default"][
            "maximum_log10_error_change"] for item in variant_cases]
        responses = [item["maximum_sentaurus_response_dex"]
                     for item in variant_cases]
        summaries.append({
            "variant": variant["id"],
            "dimension": variant["dimension"],
            "condition_count": len(variant_cases),
            "maximum_response_dex": max(responses),
            "minimum_response_dex": min(responses),
            "maximum_error_changes_dex": max_changes,
            "improves_all_conditions": all(value < 0.0 for value in max_changes),
            "worsens_all_conditions": all(value > 0.0 for value in max_changes),
            "material_in_all_conditions": all(
                value >= float(policy["material_response_min_dex"])
                for value in responses),
        })
    null = next(item for item in summaries
                if item["variant"] == "explicit_gradqf")
    report = {
        "schema": "vela.simplemos.sdevice.m8b_hfs_controls_comparison.v1",
        "status": "complete",
        "comparison_basis": (
            "fixed M8-A full-physics Vela candidate; one Sentaurus HFS "
            "control changed from the reused M8-A default reference"),
        "interpolation": "forbidden",
        "null_control_pass": (
            null["maximum_response_dex"] <=
            float(policy["null_control_max_response_dex"])),
        "variant_summaries": summaries,
        "cases": cases,
    }
    write_json(destination / "comparison_report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--reuse-prepared", action="store_true")
    parser.add_argument("--live-sentaurus", action="store_true")
    parser.add_argument("--extract-existing", action="store_true")
    parser.add_argument("--compare", action="store_true")
    parser.add_argument("--sentaurus-jobs", type=int, default=4)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=m8a.executable("ssh"))
    parser.add_argument("--scp-bin", default=m8a.executable("scp"))
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    args = parser.parse_args()

    contract_path = args.contract.resolve()
    contract = read_json(contract_path)
    validate_contract(contract)
    tdrs = {key: value.resolve() for key, value in DEFAULT_TDRS.items()}
    baseline = args.baseline_dir.resolve()
    validate_inputs(contract, tdrs, baseline)
    output = args.output_dir.resolve()
    raw_extract = args.raw_dir.resolve()
    if args.reuse_prepared:
        manifest = read_json(output / "sentaurus_matrix_manifest.json")
    else:
        manifest = prepare_sentaurus(
            contract, contract_path, tdrs, output)

    if args.live_sentaurus:
        banner = run_sentaurus(
            manifest, output, args.ssh_target, args.ssh_bin, args.scp_bin,
            args.remote_root, args.sentaurus_jobs, raw_extract)
        extract_references(contract, manifest, output, raw_extract, banner)
    elif args.extract_existing:
        banner_path = output / "sentaurus_banner.txt"
        if banner_path.is_file():
            banner = banner_path.read_text(encoding="utf-8").strip()
        else:
            banner = run([
                args.ssh_bin, "-n", args.ssh_target,
                "sdevice -h 2>&1 | sed -n '1,5p'"], capture=True).strip()
            if "T-2022.03-SP2" not in banner:
                raise RuntimeError(f"unexpected Sentaurus release:\n{banner}")
            banner_path.write_text(
                banner + "\n", encoding="utf-8", newline="\n")
        archive = (output / "sentaurus_raw"
                   / "simplemos_m8b_hfs_controls_results.tgz")
        extract_archive(archive, raw_extract)
        extract_references(contract, manifest, output, raw_extract, banner)

    comparison = compare(contract, output, baseline) if args.compare else None
    print(json.dumps({
        "sentaurus_cases": len(manifest["cases"]),
        "comparison_status": (comparison or {}).get("status"),
        "null_control_pass": (comparison or {}).get("null_control_pass"),
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
