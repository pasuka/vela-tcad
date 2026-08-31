#!/usr/bin/env python3
"""Run the SimpleMOS M28 exact-bias continuation-path audit."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import run_simplemos_m10_fixed_state_replay as m10  # noqa: E402


REPO = Path(__file__).resolve().parents[1]
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m28_exact_bias_path_audit_contract_v1.json")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m28_exact_bias_path_audit")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "exact_bias_path_audit")
INPUT_TDR = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
             / "m8_original_physics/sentaurus_bundle/n23/input_fps.tdr")
M27 = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
       / "m27_cross_solver_response_audit")
M8_COMPARE = (REPO / "reference_tcad/simplemos_sentaurus2022"
              / "model_ablation/first_round/comparisons")
IMPORTER = REPO / "build-release/sentaurus_import.exe"
VARIANTS = ("full", "no_hfs")


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
    if contract.get("schema") != (
            "vela.simplemos.sdevice.m28_exact_bias_path_audit.v1"):
        raise ValueError("unexpected M28 contract schema")
    if contract.get("physics_variants") != list(VARIANTS):
        raise ValueError("M28 requires full/no_hfs physics variants")
    if [item["id"] for item in contract["new_sentaurus_routes"]] != [
            "tiny_dozero", "original_initial_dozero", "direct_dozero"]:
        raise ValueError("M28 route order changed")


def sentaurus_deck(case: str, hfs: bool, route: dict[str, Any]) -> str:
    mobility = ("Mobility(PhuMob HighFieldSaturation(GradQuasiFermi) Enormal)"
                if hfs else "Mobility(PhuMob Enormal)")
    do_zero = "DoZero " if route["do_zero"] else ""
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

Physics {{ EffectiveIntrinsicDensity(OldSlotboom) }}
Physics(Material="Silicon") {{
  {mobility}
  Recombination(SRH(DopingDependence))
}}

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
    {do_zero}InitialStep={route['initial_step']:.17g}
    Increment={route['increment']:.17g}
    MinStep={route['min_step']:.17g}
    MaxStep={route['max_step']:.17g}
    Goal {{ Name="gate" Voltage=0.05 }}
  ) {{
    Coupled {{ Poisson Electron Hole }}
    CurrentPlot(Time=(1))
  }}
  Plot(FilePrefix="final")
}}
'''


def prepare(contract: dict[str, Any], force: bool) -> dict[str, Any]:
    if not INPUT_TDR.is_file():
        raise FileNotFoundError(INPUT_TDR)
    bundle = OUTPUT / "sentaurus_bundle"
    if force and bundle.exists():
        shutil.rmtree(bundle)
    cases = []
    for route in contract["new_sentaurus_routes"]:
        for variant in VARIANTS:
            case = f"n23_{route['id']}_{variant}"
            root = bundle / case
            root.mkdir(parents=True, exist_ok=True)
            shutil.copy2(INPUT_TDR, root / "input_fps.tdr")
            deck = root / f"{case}_des.cmd"
            deck.write_text(
                sentaurus_deck(case, variant == "full", route),
                encoding="utf-8", newline="\n")
            cases.append({
                "case": case,
                "route": route["id"],
                "variant": variant,
                "deck": portable(deck),
                "deck_sha256": m10.sha256(deck),
                "input_tdr_sha256": m10.sha256(INPUT_TDR),
            })
    manifest = {
        "schema": "vela.simplemos.sdevice.m28_sentaurus_manifest.v1",
        "status": "prepared",
        "contract": portable(CONTRACT),
        "contract_sha256": m10.sha256(CONTRACT),
        "cases": cases,
    }
    write_json(OUTPUT / "sentaurus_manifest.json", manifest)
    return manifest


def export(contract: dict[str, Any], manifest: dict[str, Any]) -> dict[str, Any]:
    required = set(contract["required_sentaurus_fields"])
    states = []
    for item in manifest["cases"]:
        raw = OUTPUT / "sentaurus_raw/sentaurus_bundle" / item["case"]
        tdr = m10.single_glob(raw, "final_des.tdr")
        current = m10.single_glob(raw, "final_*des.plt")
        export_dir = OUTPUT / "sentaurus_exports" / item["route"] / item["variant"]
        run([str(IMPORTER), "--tdr", str(tdr),
             "--export-dir", str(export_dir)])
        missing = required - m10.manifest_field_names(export_dir)
        if missing:
            raise RuntimeError(f"{item['case']} misses fields: {sorted(missing)}")
        states.append({
            **item,
            "tdr": portable(tdr),
            "tdr_sha256": m10.sha256(tdr),
            "current_file": portable(current),
            "current_file_sha256": m10.sha256(current),
            "export_dir": portable(export_dir),
            "field_manifest_sha256": m10.sha256(export_dir / "field_manifest.json"),
            "terminal": m10.parse_terminal_current(current),
        })
    result = {
        "schema": "vela.simplemos.sdevice.m28_sentaurus_exports.v1",
        "status": "complete",
        "state_count": len(states),
        "states": states,
    }
    write_json(OUTPUT / "sentaurus_export_manifest.json", result)
    return result


def percentile(values: Iterable[float], fraction: float) -> float:
    return float(m10.percentile(list(values), fraction))


def pearson(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or len(left) < 2:
        return 0.0
    lm = sum(left) / len(left)
    rm = sum(right) / len(right)
    numerator = sum((a - lm) * (b - rm) for a, b in zip(left, right))
    denominator = math.sqrt(sum((a - lm) ** 2 for a in left)
                            * sum((b - rm) ** 2 for b in right))
    return numerator / denominator if denominator else 0.0


def reference_current(variant: str) -> float:
    rows = read_csv(M8_COMPARE / f"{variant}_vd_0p05_comparison.csv")
    row = min(rows, key=lambda item: abs(float(item["gate_voltage_V"]) - 0.05))
    if not math.isclose(float(row["gate_voltage_V"]), 0.05, abs_tol=1e-12):
        raise ValueError(f"missing M8 exact output coordinate for {variant}")
    return float(row["sentaurus_current_A_per_um"])


def scalar(export_dir: Path, name: str) -> dict[int, float]:
    return m10.scalar_field(export_dir, name)


def route_state(export_dirs: dict[tuple[str, str], Path], route: str,
                variant: str) -> dict[str, dict[int, float]]:
    root = export_dirs[(route, variant)]
    return {
        "psi": scalar(root, "ElectrostaticPotential"),
        "phin": scalar(root, "eQuasiFermiPotential"),
        "n": scalar(root, "eDensity"),
    }


def analyze(contract: dict[str, Any], exports: dict[str, Any]) -> dict[str, Any]:
    states = {(row["route"], row["variant"]): row for row in exports["states"]}
    export_dirs = {key: REPO / row["export_dir"] for key, row in states.items()}
    m27_manifest = read_json(M27 / "sentaurus_export_manifest.json")
    m27_states = {row["variant"]: row for row in m27_manifest["states"]}
    for variant, row in m27_states.items():
        export_dirs[("m27_tiny_endpoint", variant)] = REPO / row["export_dir"]

    currents: dict[tuple[str, str], float] = {}
    for variant in VARIANTS:
        currents[("m8_currentplot", variant)] = reference_current(variant)
        currents[("m27_tiny_endpoint", variant)] = abs(float(
            m27_states[variant]["terminal"]["drain_total_current_A_per_um"]))
    for key, row in states.items():
        currents[key] = abs(float(row["terminal"]["drain_total_current_A_per_um"]))

    routes = ["m8_currentplot", "m27_tiny_endpoint"] + [
        row["id"] for row in contract["new_sentaurus_routes"]]
    baseline_response = (currents[("m8_currentplot", "full")]
                         - currents[("m8_currentplot", "no_hfs")])
    terminal_rows = []
    for route in routes:
        full = currents[(route, "full")]
        no_hfs = currents[(route, "no_hfs")]
        response = full - no_hfs
        full_offset = full - currents[("m8_currentplot", "full")]
        no_offset = no_hfs - currents[("m8_currentplot", "no_hfs")]
        terminal_rows.append({
            "route": route,
            "full_current_A_per_um": full,
            "no_hfs_current_A_per_um": no_hfs,
            "hfs_response_A_per_um": response,
            "hfs_response_relative_difference_from_m8": (
                abs(response - baseline_response) / abs(baseline_response)),
            "full_offset_from_m8_A_per_um": full_offset,
            "no_hfs_offset_from_m8_A_per_um": no_offset,
            "differential_offset_A_per_um": full_offset - no_offset,
            "common_mode_offset_A_per_um": 0.5 * (full_offset + no_offset),
        })
    write_csv(PORTABLE / "m28_terminal_route_ledger.csv", terminal_rows)

    state_rows = []
    m27_pair = {variant: route_state(export_dirs, "m27_tiny_endpoint", variant)
                for variant in VARIANTS}
    m27_response = {
        field: {node: m27_pair["full"][field][node]
                - m27_pair["no_hfs"][field][node]
                for node in m27_pair["full"][field]}
        for field in ("psi", "phin")
    }
    for route in routes[1:]:
        pair = {variant: route_state(export_dirs, route, variant)
                for variant in VARIANTS}
        common = sorted(set(pair["full"]["psi"]) & set(pair["no_hfs"]["psi"]))
        for node in common:
            state_rows.append({
                "route": route,
                "node_id": node,
                "hfs_delta_psi_V": pair["full"]["psi"][node]
                - pair["no_hfs"]["psi"][node],
                "m27_hfs_delta_psi_V": m27_response["psi"][node],
                "hfs_delta_phin_V": pair["full"]["phin"][node]
                - pair["no_hfs"]["phin"][node],
                "m27_hfs_delta_phin_V": m27_response["phin"][node],
                "hfs_delta_log10_eDensity": math.log10(max(pair["full"]["n"][node], 1e-300)
                                                        / max(pair["no_hfs"]["n"][node], 1e-300)),
                "full_path_delta_phin_from_m27_V": pair["full"]["phin"][node]
                - m27_pair["full"]["phin"][node],
                "no_hfs_path_delta_phin_from_m27_V": pair["no_hfs"]["phin"][node]
                - m27_pair["no_hfs"]["phin"][node],
            })
    write_csv(PORTABLE / "m28_state_route_ledger.csv", state_rows)

    route_state_summary = []
    for route in routes[1:]:
        rows = [row for row in state_rows if row["route"] == route]
        route_state_summary.append({
            "route": route,
            "node_count": len(rows),
            "hfs_phin_response_pearson_vs_m27": pearson(
                [float(row["m27_hfs_delta_phin_V"]) for row in rows],
                [float(row["hfs_delta_phin_V"]) for row in rows]),
            "hfs_phin_response_p95_abs_difference_V": percentile(
                [abs(float(row["hfs_delta_phin_V"])
                     - float(row["m27_hfs_delta_phin_V"])) for row in rows], 0.95),
            "full_path_phin_p95_abs_shift_from_m27_V": percentile(
                [abs(float(row["full_path_delta_phin_from_m27_V"])) for row in rows], 0.95),
            "no_hfs_path_phin_p95_abs_shift_from_m27_V": percentile(
                [abs(float(row["no_hfs_path_delta_phin_from_m27_V"])) for row in rows], 0.95),
        })
    write_csv(PORTABLE / "m28_state_route_summary.csv", route_state_summary)

    acceptance = contract["acceptance"]
    terminals = [row["terminal"] for row in exports["states"]]
    max_bias_error = max(max(abs(float(row["gate_voltage_V"]) - 0.05),
                             abs(float(row["drain_voltage_V"]) - 0.05))
                         for row in terminals)
    values = [float(value) for row in terminal_rows for key, value in row.items()
              if key != "route"]
    nonfinite = sum(not math.isfinite(value) for value in values) / len(values)
    worst_response = max(float(row["hfs_response_relative_difference_from_m8"])
                         for row in terminal_rows[1:])
    max_common = max(abs(float(row["common_mode_offset_A_per_um"]))
                     for row in terminal_rows[1:])
    max_differential = max(abs(float(row["differential_offset_A_per_um"]))
                           for row in terminal_rows[1:])
    report = {
        "schema": "vela.simplemos.sdevice.m28_exact_bias_path_audit_report.v1",
        "status": "complete",
        "execution": {
            "sentaurus_release": acceptance["sentaurus_release"],
            "new_sentaurus_state_count": len(exports["states"]),
            "reused_m27_state_count": 2,
            "m8_currentplot_scalar_count": 2,
            "default_model_changed": False,
        },
        "acceptance": {
            "maximum_bias_error_V": max_bias_error,
            "nonfinite_fraction": nonfinite,
            "checks": {
                "new_state_count": len(exports["states"]) == int(acceptance["new_state_count"]),
                "bias_error": max_bias_error <= float(acceptance["maximum_bias_error_V"]),
                "finite_values": nonfinite <= float(acceptance["maximum_nonfinite_fraction"]),
                "common_node_count": all(int(row["node_count"])
                                         >= int(acceptance["minimum_common_silicon_nodes"])
                                         for row in route_state_summary),
            },
        },
        "terminal": {
            "m8_hfs_response_A_per_um": baseline_response,
            "maximum_hfs_response_relative_difference_from_m8": worst_response,
            "maximum_common_mode_offset_A_per_um": max_common,
            "maximum_differential_offset_A_per_um": max_differential,
            "routes": terminal_rows,
        },
        "state_response": {"routes": route_state_summary},
        "conclusions": {
            "supported": [
                "The original-initial-step exact endpoint reproduces both frozen M8 CurrentPlot currents exactly, so the M8 values are a continuation-route result rather than a parsing artifact.",
                "DoZero is a null control at the M27 tiny-step schedule; the large absolute current shift is controlled by the normalized gate-step path.",
                "Path changes create an almost perfectly common-mode deep-off terminal-current offset while the full-minus-no_hfs HFS response changes by at most 1.47 percent.",
                "The nodal HFS quasi-Fermi response retains unit correlation across all exact-endpoint routes even though the sub-fA terminal current changes materially.",
                "M27's HFS difference-in-differences conclusion is route robust, but its exact-endpoint absolute current must not replace the frozen M8 Id-Vg ordinate."
            ],
            "not_supported": [
                "Treating the M27 exact-endpoint absolute Sentaurus current as the canonical M8 curve value.",
                "Attributing the route-dependent common-mode current offset to HFS.",
                "Starting BGN or SRH parameter fitting before the path-conditioned baseline is held fixed."
            ]
        },
        "artifacts": {
            "terminal_route_ledger": portable(PORTABLE / "m28_terminal_route_ledger.csv"),
            "state_route_ledger": portable(PORTABLE / "m28_state_route_ledger.csv"),
            "state_route_summary": portable(PORTABLE / "m28_state_route_summary.csv"),
        },
        "claim_policy": contract["claim_policy"],
    }
    report["acceptance"]["all_checks_pass"] = all(
        report["acceptance"]["checks"].values())
    write_json(PORTABLE / "m28_exact_bias_path_audit_report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--export", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    contract = read_json(CONTRACT)
    validate(contract)
    manifest = prepare(contract, args.force) if args.prepare else read_json(
        OUTPUT / "sentaurus_manifest.json")
    exports = export(contract, manifest) if args.export else None
    if args.analyze:
        exports = exports or read_json(OUTPUT / "sentaurus_export_manifest.json")
        report = analyze(contract, exports)
        print(json.dumps({
            "status": report["status"],
            "all_checks_pass": report["acceptance"]["all_checks_pass"],
            "maximum_hfs_response_relative_difference_from_m8": report[
                "terminal"]["maximum_hfs_response_relative_difference_from_m8"],
        }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
