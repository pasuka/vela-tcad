#!/usr/bin/env python3
"""Run the 16-state SimpleMOS Sentaurus-to-Vela fixed-state replay."""

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
import statistics
import subprocess
import sys
import tarfile
from typing import Any, Iterable, Sequence


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m9_refdens_scan as m9  # noqa: E402
import sentaurus_import  # noqa: E402


CONTRACT = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_m10_fixed_state_replay_contract_v1.json"
)
OUTPUT = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m10_fixed_state_replay"
)
M9_OUTPUT = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m9_hfs_diagnostics"
)
BASELINE = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m8a_confirmation"
)
RUNNER = REPO / "build-release/vela_example_runner.exe"
IMPORTER = REPO / "build-release/sentaurus_import.exe"
REMOTE_ROOT = "~/sentaurus_runs/vela_oracle/simplemos_m10_fixed_state_20260828_v1"
ARCHIVE_NAME = "simplemos_m10_fixed_state_results.tgz"
DEVICES = ("n17", "n21")
DRAIN_CASES = (("0p05", 0.05), ("1", 1.0))
GATE_CASES = (("0", 0.0), ("0p05", 0.05), ("0p8", 0.8), ("2p5", 2.5))
Q = 1.602176634e-19


def executable(name: str) -> str:
    if os.name == "nt":
        candidate = (
            Path(os.environ.get("SystemRoot", r"C:\Windows"))
            / "System32" / "OpenSSH" / f"{name}.exe"
        )
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


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8", newline="\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"refusing to write an empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO).as_posix()
    except ValueError:
        return str(resolved)


def validate_contract(contract: dict[str, Any]) -> None:
    if contract.get("schema") != "vela.simplemos.sdevice.m10_fixed_state_replay.v1":
        raise ValueError("unexpected M10 contract schema")
    policy = contract["default_model_policy"]
    if policy["modify_vela_default_model"] or not policy["diagnostics_are_read_only"]:
        raise ValueError("M10 must remain a read-only diagnostic")
    if [item["id"] for item in contract["devices"]] != list(DEVICES):
        raise ValueError("M10 devices must be n17 and n21")
    if [float(value) for value in contract["bias_matrix"]["drain_voltages_V"]] != [0.05, 1.0]:
        raise ValueError("unexpected M10 drain matrix")
    if [float(value) for value in contract["bias_matrix"]["gate_voltages_V"]] != [0.0, 0.05, 0.8, 2.5]:
        raise ValueError("unexpected M10 gate matrix")
    if int(contract["bias_matrix"]["state_count"]) != 16:
        raise ValueError("M10 must contain sixteen states")
    m9.validate_baseline_guard({"default_model_policy": policy})


def sentaurus_deck(case: str, drain_voltage: float) -> str:
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
  Mobility(PhuMob HighFieldSaturation(GradQuasiFermi) Enormal)
  Recombination(SRH(DopingDependence))
}}

Plot {{
  eDensity hDensity
  TotalCurrent/Vector eCurrent/Vector hCurrent/Vector
  eMobility hMobility eVelocity hVelocity
  eQuasiFermi hQuasiFermi
  ElectricField/Vector Potential SpaceCharge
  Doping DonorConcentration AcceptorConcentration
  SRH
  eGradQuasiFermi/Vector hGradQuasiFermi/Vector
  eEparallel hEparallel eENormal hENormal
  BandGap BandGapNarrowing Affinity
  ConductionBand ValenceBand
}}

Math {{ Extrapolate Iterations=20 ExitOnFailure }}

Solve {{
  Coupled(Iterations=100) {{ Poisson }}
  Coupled {{ Poisson Electron Hole }}
  NewCurrentPrefix="vg_0_"
  Quasistationary(
    InitialStep=0.1 Increment=1.5 MinStep=1e-5 MaxStep=1
    Goal {{ Name="drain" Voltage={drain_voltage:.17g} }}
  ) {{
    Coupled {{ Poisson Electron Hole }}
    CurrentPlot(Time=(1))
  }}
  Plot(FilePrefix="vg_0")

  NewCurrentPrefix="vg_0p05_"
  Quasistationary(
    InitialStep=0.01 Increment=1.5 MinStep=1e-5 MaxStep=0.05
    Goal {{ Name="gate" Voltage=0.05 }}
  ) {{
    Coupled {{ Poisson Electron Hole }}
    CurrentPlot(Time=(1))
  }}
  Plot(FilePrefix="vg_0p05")

  NewCurrentPrefix="vg_0p8_"
  Quasistationary(
    InitialStep=0.01 Increment=1.5 MinStep=1e-5 MaxStep=0.05
    Goal {{ Name="gate" Voltage=0.8 }}
  ) {{
    Coupled {{ Poisson Electron Hole }}
    CurrentPlot(Time=(1))
  }}
  Plot(FilePrefix="vg_0p8")

  NewCurrentPrefix="vg_2p5_"
  Quasistationary(
    InitialStep=0.01 Increment=1.5 MinStep=1e-5 MaxStep=0.05
    Goal {{ Name="gate" Voltage=2.5 }}
  ) {{
    Coupled {{ Poisson Electron Hole }}
    CurrentPlot(Time=(1))
  }}
  Plot(FilePrefix="vg_2p5")
}}
'''


def prepare(contract_path: Path, output: Path) -> dict[str, Any]:
    contract = read_json(contract_path)
    validate_contract(contract)
    bundle = output / "sentaurus_bundle"
    cases = []
    for device in DEVICES:
        source_tdr = (
            M9_OUTPUT / "sentaurus_bundle" / device / "input_fps.tdr")
        if not source_tdr.is_file():
            raise FileNotFoundError(source_tdr)
        for drain_tag, drain_voltage in DRAIN_CASES:
            case = f"{device}_vd_{drain_tag}"
            root = bundle / case
            root.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source_tdr, root / "input_fps.tdr")
            deck = root / f"{case}_des.cmd"
            deck.write_text(
                sentaurus_deck(case, drain_voltage),
                encoding="utf-8", newline="\n")
            cases.append({
                "case": case,
                "device": device,
                "drain_tag": drain_tag,
                "drain_voltage_V": drain_voltage,
                "deck": portable(deck),
                "deck_sha256": sha256(deck),
                "tdr_sha256": sha256(source_tdr),
            })
    manifest = {
        "schema": "vela.simplemos.sdevice.m10_sentaurus_manifest.v1",
        "status": "prepared",
        "contract": portable(contract_path),
        "contract_sha256": sha256(contract_path),
        "cases": cases,
    }
    write_json(output / "sentaurus_manifest.json", manifest)
    return manifest


def run_sentaurus(manifest: dict[str, Any], output: Path, ssh_target: str,
                   ssh_bin: str, scp_bin: str, remote_root: str,
                   jobs: int) -> str:
    banner = run(
        [ssh_bin, ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'"],
        capture=True).strip()
    if "T-2022.03-SP2" not in banner:
        raise RuntimeError(f"Unexpected Sentaurus release:\n{banner}")
    run([ssh_bin, ssh_target, f"mkdir -p {remote_root}"])
    run([scp_bin, "-r", str(output / "sentaurus_bundle"),
         f"{ssh_target}:{remote_root}/"])

    def execute(item: dict[str, Any]) -> None:
        case = item["case"]
        command = (
            f"set -eu; cd {remote_root}/sentaurus_bundle/{case}; "
            f"sdevice {case}_des.cmd > console.log 2>&1")
        run([ssh_bin, ssh_target, command])

    with ThreadPoolExecutor(max_workers=min(jobs, len(manifest["cases"]))) as pool:
        list(pool.map(execute, manifest["cases"]))
    run([ssh_bin, ssh_target,
         f"cd {remote_root} && tar -czf {ARCHIVE_NAME} sentaurus_bundle"])
    raw = output / "sentaurus_raw"
    raw.mkdir(parents=True, exist_ok=True)
    archive = raw / ARCHIVE_NAME
    run([scp_bin, f"{ssh_target}:{remote_root}/{ARCHIVE_NAME}", str(archive)])
    with tarfile.open(archive, "r:gz") as stream:
        stream.extractall(raw, filter="data")
    (output / "sentaurus_banner.txt").write_text(
        banner + "\n", encoding="utf-8")
    return banner


def parse_terminal_current(path: Path) -> dict[str, float]:
    text = path.read_text(errors="ignore")
    datasets = sentaurus_import.parse_quoted_list(text, "datasets")
    rows = sentaurus_import.parse_values_block(text, len(datasets))
    if not rows:
        raise RuntimeError(f"no terminal-current rows in {path}")
    values = dict(zip(datasets, rows[-1], strict=True))
    return {
        "gate_voltage_V": float(values["gate OuterVoltage"]),
        "drain_voltage_V": float(values["drain OuterVoltage"]),
        "drain_total_current_A_per_um": float(values["drain TotalCurrent"]),
        "drain_electron_current_A_per_um": float(values["drain eCurrent"]),
        "drain_hole_current_A_per_um": float(values["drain hCurrent"]),
    }


def single_glob(root: Path, pattern: str) -> Path:
    matches = sorted(root.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(
            f"expected one {pattern} in {root}, found {len(matches)}")
    return matches[0]


def manifest_field_names(export_dir: Path, region: int = 0) -> set[str]:
    manifest = read_json(export_dir / "field_manifest.json")
    return {
        str(item["name"]) for item in manifest["fields"]
        if int(item["region"]) == region
        and item.get("mapping_status") == "complete"
    }


def export_sentaurus(contract: dict[str, Any], manifest: dict[str, Any],
                     output: Path) -> dict[str, Any]:
    raw = output / "sentaurus_raw" / "sentaurus_bundle"
    states = []
    required = set(contract["required_sentaurus_fields"])
    for item in manifest["cases"]:
        root = raw / item["case"]
        for gate_tag, gate_voltage in GATE_CASES:
            tdr = single_glob(root, f"vg_{gate_tag}_des.tdr")
            current = single_glob(root, f"vg_{gate_tag}_*des.plt")
            state_id = f"{item['case']}_vg_{gate_tag}"
            export_dir = output / "sentaurus_exports" / state_id
            run([str(IMPORTER), "--tdr", str(tdr),
                 "--export-dir", str(export_dir)])
            missing = required - manifest_field_names(export_dir)
            if missing:
                raise RuntimeError(
                    f"{state_id} misses Sentaurus fields: {sorted(missing)}")
            terminal = parse_terminal_current(current)
            if not math.isclose(
                    terminal["gate_voltage_V"], gate_voltage, abs_tol=1e-10):
                raise RuntimeError(
                    f"{state_id} current bias mismatch: {terminal}")
            states.append({
                "state": state_id,
                "case": item["case"],
                "device": item["device"],
                "drain_tag": item["drain_tag"],
                "drain_voltage_V": item["drain_voltage_V"],
                "gate_tag": gate_tag,
                "gate_voltage_V": gate_voltage,
                "tdr": portable(tdr),
                "tdr_sha256": sha256(tdr),
                "terminal_current_file": portable(current),
                "terminal_current_file_sha256": sha256(current),
                "export_dir": portable(export_dir),
                "field_manifest_sha256": sha256(
                    export_dir / "field_manifest.json"),
                "sentaurus_terminal": terminal,
            })
    if len(states) != 16:
        raise RuntimeError(f"expected sixteen exported states, got {len(states)}")
    result = {
        "schema": "vela.simplemos.sdevice.m10_sentaurus_exports.v1",
        "status": "complete",
        "state_count": len(states),
        "interpolation": "forbidden",
        "states": states,
    }
    write_json(output / "sentaurus_export_manifest.json", result)
    return result


def scalar_field(export_dir: Path, name: str, region: int = 0) -> dict[int, float]:
    rows = read_csv(export_dir / "fields" / f"{name}_region{region}.csv")
    return {int(row["node_id"]): float(row["component0"]) for row in rows}


def vector_field(export_dir: Path, name: str, region: int = 0) -> dict[int, tuple[float, float]]:
    rows = read_csv(export_dir / "fields" / f"{name}_region{region}.csv")
    return {
        int(row["node_id"]): (float(row["component0"]), float(row["component1"]))
        for row in rows
    }


def merge_all_region_scalar(export_dir: Path, name: str) -> dict[int, float]:
    manifest = read_json(export_dir / "field_manifest.json")
    result: dict[int, float] = {}
    for item in manifest["fields"]:
        if item.get("name") != name or item.get("mapping_status") != "complete":
            continue
        for node, value in scalar_field(
                export_dir, name, int(item["region"])).items():
            if node in result and not math.isclose(
                    result[node], value, rel_tol=1e-8, abs_tol=1e-10):
                raise RuntimeError(
                    f"{name} has inconsistent region traces at node {node}")
            result[node] = value
    return result


def baseline_case_map() -> dict[str, dict[str, Any]]:
    manifest = read_json(
        M9_OUTPUT / "vela_edge_probes/edge_probe_manifest.json")
    return {item["case"]: item for item in manifest["cases"]}


def make_sentaurus_state(export_dir: Path, baseline_state: Path,
                         output: Path) -> dict[str, Any]:
    rows = read_csv(baseline_state)
    if not rows:
        raise ValueError(f"empty baseline state: {baseline_state}")
    fieldnames = [
        "node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"]
    required_columns = {
        "node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"}
    if not required_columns.issubset(fieldnames):
        raise ValueError(f"baseline state columns are incomplete: {baseline_state}")
    by_node = {int(row["node_id"]): row for row in rows}
    psi = merge_all_region_scalar(export_dir, "ElectrostaticPotential")
    phin = scalar_field(export_dir, "eQuasiFermiPotential")
    phip = scalar_field(export_dir, "hQuasiFermiPotential")
    electrons = scalar_field(export_dir, "eDensity")
    holes = scalar_field(export_dir, "hDensity")
    silicon = set(phin)
    if not (silicon == set(phip) == set(electrons) == set(holes)):
        raise RuntimeError("Sentaurus carrier fields do not share one node set")
    if set(psi) != set(by_node):
        missing = set(by_node) - set(psi)
        extra = set(psi) - set(by_node)
        raise RuntimeError(
            f"Sentaurus/Vela node map mismatch: missing={len(missing)}, extra={len(extra)}")
    for node, row in by_node.items():
        row["psi"] = format(psi[node], ".17g")
        if node in silicon:
            row["phin"] = format(phin[node], ".17g")
            row["phip"] = format(phip[node], ".17g")
            row["electrons_m3"] = format(electrons[node] * 1.0e6, ".17g")
            row["holes_m3"] = format(holes[node] * 1.0e6, ".17g")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for node in sorted(by_node):
            writer.writerow({name: by_node[node][name] for name in fieldnames})
    return {
        "node_count": len(by_node),
        "silicon_node_count": len(silicon),
        "state_sha256": sha256(output),
    }


def make_drive_overrides(export_dir: Path, edge_template: Path,
                         output_dir: Path) -> dict[str, Path]:
    electron = vector_field(export_dir, "eGradQuasiFermi")
    hole = vector_field(export_dir, "hGradQuasiFermi")
    edges = read_csv(edge_template)
    methods: dict[str, list[dict[str, Any]]] = {
        "mean_endpoint_magnitude": [],
        "magnitude_of_mean_vector": [],
    }
    for row in edges:
        edge_id = int(row["edge_id"])
        n0, n1 = int(row["node0"]), int(row["node1"])
        for method in methods:
            e_drive = 0.0
            h_drive = 0.0
            if n0 in electron and n1 in electron:
                if method == "mean_endpoint_magnitude":
                    e_drive = 0.5 * (
                        math.hypot(*electron[n0]) + math.hypot(*electron[n1]))
                    h_drive = 0.5 * (
                        math.hypot(*hole[n0]) + math.hypot(*hole[n1]))
                else:
                    e_drive = math.hypot(
                        0.5 * (electron[n0][0] + electron[n1][0]),
                        0.5 * (electron[n0][1] + electron[n1][1]))
                    h_drive = math.hypot(
                        0.5 * (hole[n0][0] + hole[n1][0]),
                        0.5 * (hole[n0][1] + hole[n1][1]))
            methods[method].append({
                "edge_id": edge_id,
                "electron_drive_V_m": e_drive * 100.0,
                "hole_drive_V_m": h_drive * 100.0,
            })
    result = {}
    for method, rows_out in methods.items():
        path = output_dir / f"sentaurus_drive_{method}.csv"
        write_csv(path, rows_out)
        result[method] = path
    return result


def probe_config(source: Path, state: Path, output: Path,
                 drain_voltage: float, gate_voltage: float,
                 simulation_type: str,
                 drive_override: Path | None = None,
                 provenance: str | None = None) -> dict[str, Any]:
    config = read_json(source)
    config["simulation_type"] = simulation_type
    config["state_file"] = str(state.resolve())
    config["output_csv"] = str(output.resolve())
    config.pop("sweep", None)
    config.pop("log_file", None)
    for contact in config["contacts"]:
        if contact["name"] == "drain":
            contact["bias"] = drain_voltage
        elif contact["name"] == "gate":
            contact["bias"] = gate_voltage
    if drive_override is not None:
        config["mobility_drive_override_csv"] = str(drive_override.resolve())
        config["mobility_drive_override_provenance"] = provenance
    return config


def execute_runner(config: Path, runner: Path) -> dict[str, Any]:
    environment = os.environ.copy()
    environment["PATH"] = (
        r"D:\msys64\ucrt64\bin" + os.pathsep + environment.get("PATH", ""))
    completed = subprocess.run(
        [str(runner), "--config", str(config)], cwd=REPO,
        text=True, capture_output=True, env=environment, check=False)
    (config.parent / f"{config.stem}.stdout.txt").write_text(
        completed.stdout, encoding="utf-8")
    (config.parent / f"{config.stem}.stderr.txt").write_text(
        completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(
            f"Vela probe failed for {config}: "
            f"{completed.stderr or completed.stdout}")
    return json.loads(completed.stdout.strip().splitlines()[-1])


def percentile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not ordered:
        return math.nan
    position = (len(ordered) - 1) * fraction
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def error_stats(values: Iterable[float]) -> dict[str, float | int]:
    data = [float(value) for value in values if math.isfinite(float(value))]
    return {
        "count": len(data),
        "median": percentile(data, 0.5),
        "p90": percentile(data, 0.9),
        "p95": percentile(data, 0.95),
        "maximum": max(data, default=math.nan),
    }


def log_error(reference: float, candidate: float) -> float:
    if reference <= 0.0 or candidate <= 0.0:
        return math.nan
    return abs(math.log10(candidate) - math.log10(reference))


def contact_nodes(export_dir: Path, name: str) -> set[int]:
    for row in read_csv(export_dir / "contacts.csv"):
        if row["name"] == name:
            return {int(value) for value in row["node_ids"].split(";") if value}
    raise KeyError(f"contact {name} not found in {export_dir}")


def drain_cut_current(rows: list[dict[str, str]], nodes: set[int],
                      mobility_scale: dict[int, tuple[float, float]] | None = None
                      ) -> dict[str, float | int]:
    electron_particles = 0.0
    hole_particles = 0.0
    count = 0
    for row in rows:
        n0, n1 = int(row["node0"]), int(row["node1"])
        at0, at1 = n0 in nodes, n1 in nodes
        if at0 == at1:
            continue
        sign = 1.0 if at0 else -1.0
        edge = int(row["edge_id"])
        e_scale, h_scale = ((1.0, 1.0) if mobility_scale is None
                            else mobility_scale.get(edge, (1.0, 1.0)))
        electron_particles += sign * e_scale * float(
            row["electron_particle_line_flux_per_m_s"])
        hole_particles += sign * h_scale * float(
            row["hole_particle_line_flux_per_m_s"])
        count += 1
    electron = -Q * electron_particles * 1.0e-6
    hole = Q * hole_particles * 1.0e-6
    return {
        "crossing_edge_count": count,
        "electron_A_per_um": electron,
        "hole_A_per_um": hole,
        "total_A_per_um": electron + hole,
    }


def analyze_case(state: dict[str, Any], run_dir: Path) -> dict[str, Any]:
    export_dir = REPO / state["export_dir"]
    baseline_mobility = read_csv(run_dir / "vela_drive_edge_mobility.csv")
    sent_drive_mobility = read_csv(
        run_dir / "sentaurus_mean_endpoint_magnitude_edge_mobility.csv")
    sensitivity_mobility = read_csv(
        run_dir / "sentaurus_magnitude_of_mean_vector_edge_mobility.csv")
    sg_rows = read_csv(run_dir / "vela_drive_sg_edges.csv")
    sent_mu_e = scalar_field(export_dir, "eMobility")
    sent_mu_h = scalar_field(export_dir, "hMobility")
    sent_grad_e = vector_field(export_dir, "eGradQuasiFermi")
    sent_current_e = vector_field(export_dir, "eCurrentDensity")
    sent_current_h = vector_field(export_dir, "hCurrentDensity")
    baseline_by_edge = {int(row["edge_id"]): row for row in baseline_mobility}
    primary_by_edge = {int(row["edge_id"]): row for row in sent_drive_mobility}
    sensitivity_by_edge = {int(row["edge_id"]): row for row in sensitivity_mobility}

    spatial_rows = []
    sent_projected_e: dict[int, float] = {}
    final_mobility_scale: dict[int, tuple[float, float]] = {}
    primary_drive_scale: dict[int, tuple[float, float]] = {}
    sensitivity_drive_scale: dict[int, tuple[float, float]] = {}
    for row in sg_rows:
        edge = int(row["edge_id"])
        n0, n1 = int(row["node0"]), int(row["node1"])
        if n0 not in sent_mu_e or n1 not in sent_mu_e:
            continue
        dx = float(row["x1"]) - float(row["x0"])
        dy = float(row["y1"]) - float(row["y0"])
        length = math.hypot(dx, dy)
        if length <= 0.0:
            continue
        tx, ty = dx / length, dy / length
        sent_e_mu = 0.5 * (sent_mu_e[n0] + sent_mu_e[n1]) * 1.0e-4
        sent_h_mu = 0.5 * (sent_mu_h[n0] + sent_mu_h[n1]) * 1.0e-4
        base_e_mu = float(row["electron_mobility_m2_V_s"])
        base_h_mu = float(row["hole_mobility_m2_V_s"])
        primary_e_mu = float(primary_by_edge[edge][
            "electron_final_mobility_m2_V_s"])
        primary_h_mu = float(primary_by_edge[edge][
            "hole_final_mobility_m2_V_s"])
        sensitivity_e_mu = float(sensitivity_by_edge[edge][
            "electron_final_mobility_m2_V_s"])
        sensitivity_h_mu = float(sensitivity_by_edge[edge][
            "hole_final_mobility_m2_V_s"])
        final_mobility_scale[edge] = (
            sent_e_mu / base_e_mu if base_e_mu > 0.0 else 1.0,
            sent_h_mu / base_h_mu if base_h_mu > 0.0 else 1.0)
        primary_drive_scale[edge] = (
            primary_e_mu / base_e_mu if base_e_mu > 0.0 else 1.0,
            primary_h_mu / base_h_mu if base_h_mu > 0.0 else 1.0)
        sensitivity_drive_scale[edge] = (
            sensitivity_e_mu / base_e_mu if base_e_mu > 0.0 else 1.0,
            sensitivity_h_mu / base_h_mu if base_h_mu > 0.0 else 1.0)
        jex = 0.5 * (sent_current_e[n0][0] + sent_current_e[n1][0])
        jey = 0.5 * (sent_current_e[n0][1] + sent_current_e[n1][1])
        jhx = 0.5 * (sent_current_h[n0][0] + sent_current_h[n1][0])
        jhy = 0.5 * (sent_current_h[n0][1] + sent_current_h[n1][1])
        sent_e_line = abs((jex * tx + jey * ty) * 1.0e4 * float(row["couple_m"]))
        sent_h_line = abs((jhx * tx + jhy * ty) * 1.0e4 * float(row["couple_m"]))
        sent_projected_e[edge] = sent_e_line
        grad_mean_mag = 0.5 * (
            math.hypot(*sent_grad_e[n0]) + math.hypot(*sent_grad_e[n1]))
        spatial_rows.append({
            "edge_id": edge,
            "node0": n0,
            "node1": n1,
            "x_mid_um": 0.5e6 * (float(row["x0"]) + float(row["x1"])),
            "y_mid_um": 0.5e6 * (float(row["y0"]) + float(row["y1"])),
            "sentaurus_eGradQF_mean_magnitude_V_cm": grad_mean_mag,
            "vela_eGradQF_drive_V_cm": float(baseline_by_edge[edge][
                "electron_mobility_field_V_m"]) / 100.0,
            "sentaurus_eMobility_cm2_V_s": sent_e_mu * 1.0e4,
            "vela_eMobility_vela_drive_cm2_V_s": base_e_mu * 1.0e4,
            "vela_eMobility_sentaurus_drive_cm2_V_s": primary_e_mu * 1.0e4,
            "vela_eMobility_sentaurus_drive_sensitivity_cm2_V_s": sensitivity_e_mu * 1.0e4,
            "sentaurus_projected_eLineCurrent_A_m": sent_e_line,
            "sentaurus_projected_hLineCurrent_A_m": sent_h_line,
            "vela_sg_eLineCurrent_A_m": abs(
                Q * float(row["electron_particle_line_flux_per_m_s"])),
            "vela_sg_hLineCurrent_A_m": abs(
                Q * float(row["hole_particle_line_flux_per_m_s"])),
        })
    maximum_current = max(sent_projected_e.values(), default=0.0)
    active = {
        edge for edge, value in sent_projected_e.items()
        if value > maximum_current * 1.0e-3
    }
    for row in spatial_rows:
        edge = int(row["edge_id"])
        row["active_current_edge"] = edge in active
        row["drive_abs_error_dex"] = log_error(
            float(row["sentaurus_eGradQF_mean_magnitude_V_cm"]),
            float(row["vela_eGradQF_drive_V_cm"]))
        row["mobility_vela_drive_abs_error_dex"] = log_error(
            float(row["sentaurus_eMobility_cm2_V_s"]),
            float(row["vela_eMobility_vela_drive_cm2_V_s"]))
        row["mobility_sentaurus_drive_abs_error_dex"] = log_error(
            float(row["sentaurus_eMobility_cm2_V_s"]),
            float(row["vela_eMobility_sentaurus_drive_cm2_V_s"]))
        row["mobility_sentaurus_drive_sensitivity_abs_error_dex"] = log_error(
            float(row["sentaurus_eMobility_cm2_V_s"]),
            float(row["vela_eMobility_sentaurus_drive_sensitivity_cm2_V_s"]))
        row["sg_line_current_abs_error_dex"] = log_error(
            float(row["sentaurus_projected_eLineCurrent_A_m"]),
            float(row["vela_sg_eLineCurrent_A_m"]))
    write_csv(run_dir / "edge_replay_comparison.csv", spatial_rows)

    def selected_stats(column: str, active_only: bool) -> dict[str, float | int]:
        return error_stats(
            float(row[column]) for row in spatial_rows
            if (not active_only or row["active_current_edge"])
            and math.isfinite(float(row[column])))

    drain_nodes = contact_nodes(export_dir, "drain")
    terminal_native = drain_cut_current(sg_rows, drain_nodes)
    terminal_sent_drive = drain_cut_current(
        sg_rows, drain_nodes, primary_drive_scale)
    terminal_sent_drive_sensitivity = drain_cut_current(
        sg_rows, drain_nodes, sensitivity_drive_scale)
    terminal_sent_mobility = drain_cut_current(
        sg_rows, drain_nodes, final_mobility_scale)
    sent_terminal = state["sentaurus_terminal"]
    sent_current = abs(float(sent_terminal["drain_total_current_A_per_um"]))
    terminal_variants = {
        "vela_drive_vela_hfs": terminal_native,
        "sentaurus_drive_vela_hfs": terminal_sent_drive,
        "sentaurus_drive_vela_hfs_sensitivity": terminal_sent_drive_sensitivity,
        "sentaurus_final_mobility": terminal_sent_mobility,
    }
    for value in terminal_variants.values():
        value["absolute_log10_error_dex"] = log_error(
            sent_current, abs(float(value["total_A_per_um"])))
        value["relative_error"] = (
            abs(abs(float(value["total_A_per_um"])) - sent_current)
            / max(sent_current, 1.0e-300))
    return {
        "state": state["state"],
        "device": state["device"],
        "drain_voltage_V": state["drain_voltage_V"],
        "gate_voltage_V": state["gate_voltage_V"],
        "silicon_edge_count": len(spatial_rows),
        "active_current_edge_count": len(active),
        "drive": {
            "all_edges_abs_error_dex": selected_stats(
                "drive_abs_error_dex", False),
            "active_edges_abs_error_dex": selected_stats(
                "drive_abs_error_dex", True),
        },
        "electron_mobility": {
            "vela_drive": {
                "all_edges_abs_error_dex": selected_stats(
                    "mobility_vela_drive_abs_error_dex", False),
                "active_edges_abs_error_dex": selected_stats(
                    "mobility_vela_drive_abs_error_dex", True),
            },
            "sentaurus_drive": {
                "all_edges_abs_error_dex": selected_stats(
                    "mobility_sentaurus_drive_abs_error_dex", False),
                "active_edges_abs_error_dex": selected_stats(
                    "mobility_sentaurus_drive_abs_error_dex", True),
            },
            "sentaurus_drive_sensitivity": {
                "all_edges_abs_error_dex": selected_stats(
                    "mobility_sentaurus_drive_sensitivity_abs_error_dex", False),
                "active_edges_abs_error_dex": selected_stats(
                    "mobility_sentaurus_drive_sensitivity_abs_error_dex", True),
            },
        },
        "projected_electron_line_current": {
            "all_edges_abs_error_dex": selected_stats(
                "sg_line_current_abs_error_dex", False),
            "active_edges_abs_error_dex": selected_stats(
                "sg_line_current_abs_error_dex", True),
            "interpretation": (
                "Sentaurus nodal current projected onto Vela primal edges; "
                "localization diagnostic, not an exact discretization identity"),
        },
        "sentaurus_terminal_current_A_per_um": float(
            sent_terminal["drain_total_current_A_per_um"]),
        "terminal_current_replay": terminal_variants,
        "artifacts": {
            "edge_comparison": portable(
                run_dir / "edge_replay_comparison.csv"),
        },
    }


def execute_replay(export_manifest: dict[str, Any], output: Path,
                   runner: Path, jobs: int) -> dict[str, Any]:
    if not runner.is_file():
        raise FileNotFoundError(runner)
    baseline = baseline_case_map()

    def execute(state: dict[str, Any]) -> dict[str, Any]:
        state_id = state["state"]
        root = output / "replay" / state_id
        root.mkdir(parents=True, exist_ok=True)
        baseline_id = (
            f"{state['device']}_vd_{state['drain_tag']}_vg_{state['gate_tag']}")
        source = baseline[baseline_id]
        baseline_state = Path(source["state_file"])
        state_file = root / "sentaurus_state_for_vela.csv"
        state_info = make_sentaurus_state(
            REPO / state["export_dir"], baseline_state, state_file)
        drive_files = make_drive_overrides(
            REPO / state["export_dir"], Path(source["edge_csv"]), root)
        probes = [
            ("vela_drive_edge_mobility", "edge_mobility_probe", None, None),
            ("vela_drive_sg_edges", "sg_edge_flux_probe", None, None),
            ("sentaurus_mean_endpoint_magnitude_edge_mobility",
             "edge_mobility_probe", drive_files["mean_endpoint_magnitude"],
             "sentaurus_mean_endpoint_grad_qf_magnitude"),
            ("sentaurus_magnitude_of_mean_vector_edge_mobility",
             "edge_mobility_probe", drive_files["magnitude_of_mean_vector"],
             "sentaurus_magnitude_of_mean_endpoint_grad_qf_vector"),
        ]
        statuses = {}
        for label, simulation_type, override, provenance in probes:
            result_csv = root / f"{label}.csv"
            config_path = root / f"{label}.json"
            write_json(config_path, probe_config(
                Path(source["probe_config"]), state_file, result_csv,
                float(state["drain_voltage_V"]),
                float(state["gate_voltage_V"]), simulation_type,
                override, provenance))
            statuses[label] = execute_runner(config_path, runner)
        result = analyze_case(state, root)
        result["state_import"] = state_info
        result["probe_statuses"] = statuses
        write_json(root / "replay_summary.json", result)
        return result

    with ThreadPoolExecutor(max_workers=min(jobs, len(export_manifest["states"]))) as pool:
        cases = list(pool.map(execute, export_manifest["states"]))
    report = {
        "schema": "vela.simplemos.sdevice.m10_fixed_state_replay_evidence.v1",
        "status": "complete",
        "case_count": len(cases),
        "default_model_changed": False,
        "sentaurus_release": "T-2022.03-SP2",
        "cases": cases,
    }
    write_json(output / "fixed_state_replay_report.json", report)
    write_case_summary(output / "fixed_state_replay_summary.csv", report)
    return report


def write_case_summary(path: Path, report: dict[str, Any]) -> None:
    rows = []
    for item in report["cases"]:
        replay = item["terminal_current_replay"]
        rows.append({
            "state": item["state"],
            "device": item["device"],
            "drain_voltage_V": item["drain_voltage_V"],
            "gate_voltage_V": item["gate_voltage_V"],
            "active_current_edge_count": item["active_current_edge_count"],
            "drive_active_p95_error_dex": item["drive"][
                "active_edges_abs_error_dex"]["p95"],
            "mobility_vela_drive_active_p95_error_dex": item[
                "electron_mobility"]["vela_drive"][
                    "active_edges_abs_error_dex"]["p95"],
            "mobility_sentaurus_drive_active_p95_error_dex": item[
                "electron_mobility"]["sentaurus_drive"][
                    "active_edges_abs_error_dex"]["p95"],
            "sentaurus_terminal_current_A_per_um": item[
                "sentaurus_terminal_current_A_per_um"],
            "terminal_vela_drive_error_dex": replay[
                "vela_drive_vela_hfs"]["absolute_log10_error_dex"],
            "terminal_sentaurus_drive_error_dex": replay[
                "sentaurus_drive_vela_hfs"]["absolute_log10_error_dex"],
            "terminal_sentaurus_mobility_error_dex": replay[
                "sentaurus_final_mobility"]["absolute_log10_error_dex"],
        })
    write_csv(path, rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    parser.add_argument("--runner", type=Path, default=RUNNER)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--live-sentaurus", action="store_true")
    parser.add_argument("--extract-sentaurus", action="store_true")
    parser.add_argument("--replay", action="store_true")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--scp-bin", default=executable("scp"))
    parser.add_argument("--remote-root", default=REMOTE_ROOT)
    args = parser.parse_args()
    if args.jobs < 1:
        raise ValueError("jobs must be positive")
    if not any((args.prepare, args.live_sentaurus,
                args.extract_sentaurus, args.replay)):
        parser.error("select at least one action")
    contract_path = args.contract.resolve()
    contract = read_json(contract_path)
    validate_contract(contract)
    output = args.output_dir.resolve()
    manifest_path = output / "sentaurus_manifest.json"
    manifest = (prepare(contract_path, output)
                if args.prepare or not manifest_path.is_file()
                else read_json(manifest_path))
    if args.live_sentaurus:
        run_sentaurus(
            manifest, output, args.ssh_target, args.ssh_bin,
            args.scp_bin, args.remote_root, args.jobs)
    export_path = output / "sentaurus_export_manifest.json"
    exports = (export_sentaurus(contract, manifest, output)
               if args.extract_sentaurus
               else read_json(export_path) if export_path.is_file()
               else None)
    report = None
    if args.replay:
        if exports is None:
            raise FileNotFoundError(export_path)
        report = execute_replay(
            exports, output, args.runner.resolve(), args.jobs)
    print(json.dumps({
        "status": (report or {}).get("status", "prepared"),
        "sentaurus_cases": len(manifest["cases"]),
        "exported_states": len((exports or {}).get("states", [])),
        "replayed_states": len((report or {}).get("cases", [])),
        "default_model_changed": False,
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
