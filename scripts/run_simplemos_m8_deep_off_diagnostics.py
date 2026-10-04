#!/usr/bin/env python3
"""Collect matched n17/n21 SimpleMOS M8 deep-off diagnostic states."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import tarfile
from typing import Any, Sequence


REPO = Path(__file__).resolve().parents[1]
M8_ROOT = (
    REPO / "build-release" / "reference_tcad" / "simplemos_sentaurus2022"
    / "m8_original_physics"
)
OUTPUT = (
    REPO / "build-release" / "reference_tcad" / "simplemos_sentaurus2022"
    / "m8_deep_off_diagnostics"
)
IMPORTER = REPO / "build-release" / "sentaurus_import.exe"
RUNNER = REPO / "build-release" / "vela_example_runner.exe"
REMOTE_ROOT = "~/sentaurus_runs/vela_oracle/simplemos_m8_deep_off_20260827"
DEVICES = ("n17", "n21")
DRAIN_CASES = (("0p05", 0.05), ("1", 1.0))


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


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sentaurus_deck(device: str, drain_voltage: float) -> str:
    return f'''File {{
  Grid="input_fps.tdr"
  Plot="{device}_diagnostic_des.tdr"
  Current="{device}_diagnostic"
  Output="{device}_diagnostic.log"
}}

Electrode {{
  {{ Name="source" Voltage=0.0 }}
  {{ Name="drain" Voltage=0.0 }}
  {{ Name="gate" Voltage=0.0 }}
  {{ Name="substrate" Voltage=0.0 }}
}}

Physics {{ EffectiveIntrinsicDensity(OldSlotboom) }}
Physics(Material="Silicon") {{
  Mobility(PhuMob HighFieldSaturation Enormal)
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
  Plot(FilePrefix="equilibrium")
  Quasistationary(
    InitialStep=0.1 Increment=1.5 MinStep=1e-5 MaxStep=1
    Goal {{ Name="drain" Voltage={drain_voltage:.17g} }}
  ) {{ Coupled {{ Poisson Electron Hole }} }}
  Plot(FilePrefix="drain_state")
  NewCurrentPrefix="IdVg_"
  Quasistationary(
    DoZero InitialStep=0.01 Increment=1.5 MinStep=1e-5 MaxStep=0.05
    Goal {{ Name="gate" Voltage=0.05 }}
  ) {{
    Coupled {{ Poisson Electron Hole }}
    CurrentPlot(Time=(0; 1))
    Plot(FilePrefix="deepoff" NoOverWrite Time=(0; 1))
  }}
}}
'''


def prepare_sentaurus() -> list[dict[str, Any]]:
    bundle = OUTPUT / "sentaurus_bundle"
    cases: list[dict[str, Any]] = []
    for device in DEVICES:
        source_tdr = M8_ROOT / "sentaurus_bundle" / device / "input_fps.tdr"
        if not source_tdr.is_file():
            raise FileNotFoundError(source_tdr)
        for drain_tag, drain_voltage in DRAIN_CASES:
            case = f"{device}_vd_{drain_tag}"
            case_dir = bundle / case
            case_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source_tdr, case_dir / "input_fps.tdr")
            deck = case_dir / f"{case}_des.cmd"
            deck.write_text(
                sentaurus_deck(case, drain_voltage), encoding="utf-8", newline="\n")
            cases.append({
                "case": case,
                "device": device,
                "drain_voltage_V": drain_voltage,
                "deck": str(deck.resolve()),
                "deck_sha256": sha256(deck),
                "tdr_sha256": sha256(source_tdr),
            })
    write_json(OUTPUT / "sentaurus_manifest.json", {
        "schema": "vela.simplemos.sdevice.m8_deep_off_sentaurus.v1",
        "status": "prepared",
        "cases": cases,
    })
    return cases


def run_sentaurus(cases: list[dict[str, Any]], ssh_target: str,
                   ssh_bin: str, scp_bin: str) -> None:
    banner = run(
        [ssh_bin, ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'"],
        capture=True).strip()
    if "T-2022.03-SP2" not in banner:
        raise RuntimeError(f"Unexpected Sentaurus release:\n{banner}")
    run([ssh_bin, ssh_target, f"mkdir -p {REMOTE_ROOT}"])
    run([scp_bin, "-r", str(OUTPUT / "sentaurus_bundle"),
         f"{ssh_target}:{REMOTE_ROOT}/"])
    for index, case in enumerate(cases, start=1):
        print(f"[{index}/{len(cases)}] Sentaurus {case['case']}", flush=True)
        run([ssh_bin, ssh_target,
             f"set -eu; cd {REMOTE_ROOT}/sentaurus_bundle/{case['case']}; "
             f"sdevice {case['case']}_des.cmd > console.log 2>&1"])
    archive_name = "simplemos_m8_deep_off_results.tgz"
    run([ssh_bin, ssh_target,
         f"cd {REMOTE_ROOT} && tar -czf {archive_name} sentaurus_bundle"])
    raw = OUTPUT / "sentaurus_raw"
    raw.mkdir(parents=True, exist_ok=True)
    archive = raw / archive_name
    run([scp_bin, f"{ssh_target}:{REMOTE_ROOT}/{archive_name}", str(archive)])
    with tarfile.open(archive, "r:gz") as stream:
        stream.extractall(raw, filter="data")
    (OUTPUT / "sentaurus_banner.txt").write_text(
        banner + "\n", encoding="utf-8")


def export_sentaurus(cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    raw = OUTPUT / "sentaurus_raw" / "sentaurus_bundle"
    exports: list[dict[str, Any]] = []
    for case in cases:
        case_raw = raw / case["case"]
        equilibrium = sorted(case_raw.glob("equilibrium*.tdr"))
        drain = sorted(case_raw.glob("drain_state*.tdr"))
        deepoff = sorted(case_raw.glob("deepoff*.tdr"))
        if len(equilibrium) != 1 or len(drain) != 1 or len(deepoff) != 2:
            raise RuntimeError(
                f"Unexpected snapshots for {case['case']}: "
                f"equilibrium={len(equilibrium)}, drain={len(drain)}, "
                f"deepoff={len(deepoff)}")
        snapshots = [
            ("equilibrium", None, equilibrium[0]),
            ("drain_vg_0", 0.0, deepoff[0]),
            ("drain_vg_0p05", 0.05, deepoff[1]),
        ]
        for state, gate_voltage, tdr in snapshots:
            export_dir = OUTPUT / "sentaurus_exports" / case["case"] / state
            run([str(IMPORTER), "--tdr", str(tdr), "--export-dir", str(export_dir)])
            exports.append({
                "case": case["case"],
                "state": state,
                "gate_voltage_V": gate_voltage,
                "tdr": str(tdr.resolve()),
                "tdr_sha256": sha256(tdr),
                "export_dir": str(export_dir.resolve()),
                "field_manifest": str((export_dir / "field_manifest.json").resolve()),
            })
    write_json(OUTPUT / "sentaurus_export_manifest.json", {
        "schema": "vela.simplemos.sdevice.m8_deep_off_exports.v1",
        "status": "complete",
        "states": exports,
    })
    return exports


def diagnostics_block(root: Path) -> dict[str, Any]:
    return {
        "contact_edge": {
            "enabled": True,
            "contacts": ["source", "drain", "substrate"],
            "csv_file": str(root / "contact_edges.csv"),
        },
        "terminal_balance": {
            "enabled": True,
            "contacts": ["source", "drain", "gate", "substrate"],
            "csv_file": str(root / "terminal_balance.csv"),
        },
        "terminal_current_method_compare": {
            "enabled": True,
            "contacts": ["source", "drain", "substrate"],
            "csv_file": str(root / "terminal_current_methods.csv"),
        },
        "continuity_balance": {
            "enabled": True,
            "contacts": ["source", "drain", "substrate"],
            "csv_file": str(root / "continuity_balance.csv"),
        },
        "srh_balance": {
            "enabled": True,
            "material": "Si",
            "drain_contact": "drain",
            "substrate_contact": "substrate",
            "kcl_contacts": ["source", "drain", "gate", "substrate"],
            "resolution_margin_ratio": 10.0,
            "csv_file": str(root / "srh_balance.csv"),
        },
    }


def prepare_vela() -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for device in DEVICES:
        for drain_tag, drain_voltage in DRAIN_CASES:
            source_dir = M8_ROOT / "vela" / device / "workflow" / f"vd_{drain_tag}"
            source_config = read_json(source_dir / "20_gate_sweep.json")
            drain_state = source_dir / "10_drain_ramp_accepted_state.csv"
            for state, points in (("vg_0", [0.0]), ("vg_0p05", [0.0, 0.05])):
                root = OUTPUT / "vela" / f"{device}_vd_{drain_tag}" / state
                root.mkdir(parents=True, exist_ok=True)
                config = json.loads(json.dumps(source_config))
                config["output_csv"] = str(root / "iv.csv")
                config["log_file"] = str(root / "run.log")
                config["sweep"].update({
                    "start": points[0],
                    "stop": points[-1],
                    "step": 0.05,
                    "bias_points": points,
                    "initial_state_file": str(drain_state.resolve()),
                    "write_state_file": str((root / "accepted_state.csv").resolve()),
                    "write_vtk": True,
                    "vtk_prefix": str((root / "vtk" / "state").resolve()),
                    "diagnostics": diagnostics_block(root.resolve()),
                })
                config["simplemos_m8_repair"] = {
                    "diagnostic_only": True,
                    "source_drain_voltage_V": drain_voltage,
                    "target_gate_voltage_V": points[-1],
                }
                config_path = root / "simulation.json"
                write_json(config_path, config)
                cases.append({
                    "case": f"{device}_vd_{drain_tag}_{state}",
                    "device": device,
                    "drain_voltage_V": drain_voltage,
                    "gate_voltage_V": points[-1],
                    "config": str(config_path.resolve()),
                    "config_sha256": sha256(config_path),
                })
    write_json(OUTPUT / "vela_diagnostic_manifest.json", {
        "schema": "vela.simplemos.sdevice.m8_deep_off_vela.v1",
        "status": "prepared",
        "cases": cases,
    })
    return cases


def run_vela(cases: list[dict[str, Any]]) -> None:
    for index, case in enumerate(cases, start=1):
        print(f"[{index}/{len(cases)}] Vela {case['case']}", flush=True)
        run([str(RUNNER), "--config", case["config"]])


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def csv_field(path: Path) -> dict[int, tuple[float, ...]]:
    values: dict[int, tuple[float, ...]] = {}
    for row in csv_rows(path):
        values[int(row["node_id"])] = tuple(
            float(row[key]) for key in sorted(row) if key.startswith("component"))
    return values


def percentile(values: Sequence[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    position = fraction * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def voltage_error(reference: dict[int, float], candidate: dict[int, float],
                  node_ids: set[int]) -> dict[str, float | int]:
    errors = [candidate[node] - reference[node] for node in node_ids]
    absolute = [abs(value) for value in errors]
    return {
        "nodes": len(errors),
        "signed_median_V": statistics.median(errors),
        "absolute_median_V": statistics.median(absolute),
        "absolute_p90_V": percentile(absolute, 0.90),
        "absolute_max_V": max(absolute),
    }


def log_ratio_error(reference: dict[int, float], candidate: dict[int, float],
                    node_ids: set[int], *, floor: float = 1e-300
                    ) -> dict[str, float | int]:
    ratios = [
        math.log10(max(abs(candidate[node]), floor) /
                   max(abs(reference[node]), floor))
        for node in node_ids
    ]
    absolute = [abs(value) for value in ratios]
    return {
        "nodes": len(ratios),
        "median_log10_vela_over_sentaurus": statistics.median(ratios),
        "p10_log10_vela_over_sentaurus": percentile(ratios, 0.10),
        "p90_log10_vela_over_sentaurus": percentile(ratios, 0.90),
        "absolute_median_dex": statistics.median(absolute),
        "absolute_p90_dex": percentile(absolute, 0.90),
        "absolute_max_dex": max(absolute),
    }


def vector_magnitudes(values: dict[int, tuple[float, ...]]) -> dict[int, float]:
    return {
        node: math.sqrt(sum(component * component for component in vector))
        for node, vector in values.items()
    }


def read_vtk_point_fields(path: Path, requested: set[str]
                          ) -> dict[str, dict[int, tuple[float, ...]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    try:
        start = next(index for index, line in enumerate(lines)
                     if line.startswith("POINT_DATA "))
    except StopIteration as error:
        raise RuntimeError(f"POINT_DATA missing from {path}") from error
    count = int(lines[start].split()[1])
    fields: dict[str, dict[int, tuple[float, ...]]] = {}
    index = start + 1
    while index < len(lines):
        tokens = lines[index].split()
        if not tokens:
            index += 1
            continue
        if tokens[0] == "SCALARS":
            name = tokens[1]
            components = int(tokens[3]) if len(tokens) >= 4 else 1
            index += 2
        elif tokens[0] == "VECTORS":
            name = tokens[1]
            components = 3
            index += 1
        else:
            index += 1
            continue
        if index + count > len(lines):
            raise RuntimeError(f"Truncated VTK field {name} in {path}")
        if name in requested:
            field: dict[int, tuple[float, ...]] = {}
            for node in range(count):
                values = tuple(float(value) for value in lines[index + node].split())
                if len(values) != components:
                    raise RuntimeError(
                        f"Unexpected component count for {name} node {node}")
                field[node] = values
            fields[name] = field
        index += count
    missing = requested - fields.keys()
    if missing:
        raise RuntimeError(f"Missing VTK fields in {path}: {sorted(missing)}")
    return fields


def selector_nodes(export_dir: Path) -> tuple[dict[str, set[int]],
                                              dict[int, tuple[float, float]]]:
    coordinates = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in csv_rows(export_dir / "nodes.csv")
    }
    silicon: set[int] = set()
    oxide: set[int] = set()
    for row in csv_rows(export_dir / "elements.csv"):
        nodes = {int(row["node0"]), int(row["node1"]), int(row["node2"])}
        if row["material"] == "Si":
            silicon.update(nodes)
        if row["region"] == "Oxide_1":
            oxide.update(nodes)
    interface = silicon & oxide
    channel = {node for node in interface if abs(coordinates[node][1]) <= 0.125}
    if not channel:
        raise RuntimeError(f"No central channel interface nodes in {export_dir}")
    surface_x = statistics.median(coordinates[node][0] for node in channel)
    near_channel = {
        node for node in silicon
        if abs(coordinates[node][1]) <= 0.25
        and abs(coordinates[node][0] - surface_x) <= 0.10
    }
    return {
        "all_silicon": silicon,
        "si_oxide_interface": interface,
        "central_channel_interface": channel,
        "near_channel_silicon": near_channel,
    }, coordinates


def scalar_component(values: dict[int, tuple[float, ...]],
                     scale: float = 1.0) -> dict[int, float]:
    return {node: vector[0] * scale for node, vector in values.items()}


def integrate_silicon_nodal_rate(export_dir: Path,
                                 values_cm3_s: dict[int, float]) -> float:
    coordinates = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in csv_rows(export_dir / "nodes.csv")
    }
    integral_cm3_s_per_um = 0.0
    for row in csv_rows(export_dir / "elements.csv"):
        if row["material"] != "Si":
            continue
        nodes = (int(row["node0"]), int(row["node1"]), int(row["node2"]))
        (x0, y0), (x1, y1), (x2, y2) = (
            coordinates[node] for node in nodes)
        area_um2 = 0.5 * abs(
            (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))
        mean_rate = sum(values_cm3_s[node] for node in nodes) / 3.0
        # 2-D TCAD current is reported per micrometre of out-of-plane width.
        integral_cm3_s_per_um += mean_rate * area_um2 * 1.0e-12
    return 1.602176634e-19 * integral_cm3_s_per_um


def independent_srh_rate(n: dict[int, float], p: dict[int, float],
                         ni: dict[int, float], donors: dict[int, float],
                         acceptors: dict[int, float], *, basis: str,
                         electron_tau_max_s: float
                         ) -> dict[int, float]:
    rate: dict[int, float] = {}
    for node in n:
        if basis == "total_impurity":
            doping = abs(donors[node]) + abs(acceptors[node])
        elif basis == "net_doping":
            doping = abs(donors[node] - acceptors[node])
        else:
            raise ValueError(basis)
        tau_n = electron_tau_max_s / (1.0 + doping / 1.0e16)
        tau_p = 3.0e-6 / (1.0 + doping / 1.0e16)
        denominator = tau_p * (n[node] + ni[node]) + tau_n * (p[node] + ni[node])
        rate[node] = (n[node] * p[node] - ni[node] * ni[node]) / denominator
    return rate


def analyze_fields() -> dict[str, Any]:
    requested_vtk = {
        "SRHRecombinationCm3PerS",
        "BandgapNarrowing",
        "ElectronMobilityCm2PerVs",
        "HoleMobilityCm2PerVs",
        "SentaurusElectronCurrentDensityVector",
        "SentaurusHoleCurrentDensityVector",
    }
    cases: list[dict[str, Any]] = []
    detail_rows: list[dict[str, Any]] = []
    for device in DEVICES:
        for drain_tag, drain_voltage in DRAIN_CASES:
            case_name = f"{device}_vd_{drain_tag}"
            for state, vela_state_name in (
                    ("equilibrium", None),
                    ("drain_vg_0", "vg_0"),
                    ("drain_vg_0p05", "vg_0p05")):
                export_dir = OUTPUT / "sentaurus_exports" / case_name / state
                selectors, coordinates = selector_nodes(export_dir)
                sentaurus = {
                    "psi": scalar_component(csv_field(
                        export_dir / "fields" / "ElectrostaticPotential_region0.csv")),
                    "electrons": scalar_component(csv_field(
                        export_dir / "fields" / "eDensity_region0.csv")),
                    "holes": scalar_component(csv_field(
                        export_dir / "fields" / "hDensity_region0.csv")),
                    "phin": scalar_component(csv_field(
                        export_dir / "fields" / "eQuasiFermiPotential_region0.csv")),
                    "phip": scalar_component(csv_field(
                        export_dir / "fields" / "hQuasiFermiPotential_region0.csv")),
                }
                if vela_state_name is None:
                    state_path = (M8_ROOT / "vela" / device / "workflow" /
                                  f"vd_{drain_tag}" /
                                  "00_equilibrium_accepted_state.csv")
                    vtk_fields: dict[str, dict[int, tuple[float, ...]]] = {}
                else:
                    vela_root = OUTPUT / "vela" / case_name / vela_state_name
                    state_path = vela_root / "accepted_state.csv"
                    vtk_paths = sorted((vela_root / "vtk").glob("*.vtk"))
                    if not vtk_paths:
                        raise FileNotFoundError(f"No VTK states in {vela_root / 'vtk'}")
                    vtk_path = vtk_paths[-1]
                    vtk_fields = read_vtk_point_fields(vtk_path, requested_vtk)
                state_rows = csv_rows(state_path)
                vela = {
                    "psi": {int(row["node_id"]): float(row["psi"])
                            for row in state_rows},
                    "electrons": {int(row["node_id"]):
                                  float(row["electrons_m3"]) / 1e6
                                  for row in state_rows},
                    "holes": {int(row["node_id"]):
                              float(row["holes_m3"]) / 1e6
                              for row in state_rows},
                    "phin": {int(row["node_id"]): float(row["phin"])
                             for row in state_rows},
                    "phip": {int(row["node_id"]): float(row["phip"])
                             for row in state_rows},
                }
                selector_metrics: dict[str, Any] = {}
                integrated_srh: dict[str, float] | None = None
                if vtk_fields:
                    sentaurus_srh = scalar_component(csv_field(
                        export_dir / "fields" / "srhRecombination_region0.csv"))
                    vela_srh = scalar_component(
                        vtk_fields["SRHRecombinationCm3PerS"])
                    sentaurus_integral = integrate_silicon_nodal_rate(
                        export_dir, sentaurus_srh)
                    vela_integral = integrate_silicon_nodal_rate(
                        export_dir, vela_srh)
                    integrated_srh = {
                        "sentaurus_A_per_um": sentaurus_integral,
                        "vela_A_per_um": vela_integral,
                        "signed_ratio_vela_over_sentaurus": (
                            vela_integral / sentaurus_integral),
                    }
                    equilibrium_dir = (OUTPUT / "sentaurus_exports" / case_name /
                                       "equilibrium" / "fields")
                    equilibrium_n = scalar_component(csv_field(
                        equilibrium_dir / "eDensity_region0.csv"))
                    equilibrium_p = scalar_component(csv_field(
                        equilibrium_dir / "hDensity_region0.csv"))
                    ni = {node: math.sqrt(equilibrium_n[node] * equilibrium_p[node])
                          for node in equilibrium_n}
                    donors = scalar_component(csv_field(
                        export_dir / "fields" / "DonorConcentration_region0.csv"))
                    acceptors = scalar_component(csv_field(
                        export_dir / "fields" / "AcceptorConcentration_region0.csv"))
                    integrated_srh[
                        "independent_transportmodels_override_A_per_um"] = (
                        integrate_silicon_nodal_rate(
                            export_dir, independent_srh_rate(
                                sentaurus["electrons"], sentaurus["holes"], ni,
                                donors, acceptors, basis="total_impurity",
                                electron_tau_max_s=3.0e-8)))
                    integrated_srh[
                        "independent_simplemos_default_A_per_um"] = (
                        integrate_silicon_nodal_rate(
                            export_dir, independent_srh_rate(
                                sentaurus["electrons"], sentaurus["holes"], ni,
                                donors, acceptors, basis="total_impurity",
                                electron_tau_max_s=1.0e-5)))
                for selector, nodes in selectors.items():
                    metrics: dict[str, Any] = {
                        "potential": voltage_error(
                            sentaurus["psi"], vela["psi"], nodes),
                        "electron_density": log_ratio_error(
                            sentaurus["electrons"], vela["electrons"], nodes),
                        "hole_density": log_ratio_error(
                            sentaurus["holes"], vela["holes"], nodes),
                        "electron_quasi_fermi": voltage_error(
                            sentaurus["phin"], vela["phin"], nodes),
                        "hole_quasi_fermi": voltage_error(
                            sentaurus["phip"], vela["phip"], nodes),
                    }
                    if vtk_fields:
                        scalar_pairs = (
                            ("srh", "srhRecombination", "SRHRecombinationCm3PerS"),
                            ("bandgap_narrowing", "BandgapNarrowing", "BandgapNarrowing"),
                            ("electron_mobility", "eMobility", "ElectronMobilityCm2PerVs"),
                            ("hole_mobility", "hMobility", "HoleMobilityCm2PerVs"),
                        )
                        for label, sentaurus_name, vela_name in scalar_pairs:
                            reference = scalar_component(csv_field(
                                export_dir / "fields" /
                                f"{sentaurus_name}_region0.csv"))
                            candidate = scalar_component(vtk_fields[vela_name])
                            metrics[label] = log_ratio_error(
                                reference, candidate, nodes, floor=1e-100)
                        for carrier, sentaurus_name, vela_name in (
                                ("electron_current_density", "eCurrentDensity",
                                 "SentaurusElectronCurrentDensityVector"),
                                ("hole_current_density", "hCurrentDensity",
                                 "SentaurusHoleCurrentDensityVector")):
                            reference = vector_magnitudes(csv_field(
                                export_dir / "fields" /
                                f"{sentaurus_name}_region0.csv"))
                            candidate = vector_magnitudes(vtk_fields[vela_name])
                            metrics[carrier] = log_ratio_error(
                                reference, candidate, nodes, floor=1e-100)
                    selector_metrics[selector] = metrics
                channel = selectors["central_channel_interface"]
                for node in sorted(channel, key=lambda item: coordinates[item][1]):
                    detail_rows.append({
                        "case": case_name,
                        "state": state,
                        "node_id": node,
                        "x_um": coordinates[node][0],
                        "y_um": coordinates[node][1],
                        "sentaurus_psi_V": sentaurus["psi"][node],
                        "vela_psi_V": vela["psi"][node],
                        "sentaurus_eDensity_cm3": sentaurus["electrons"][node],
                        "vela_eDensity_cm3": vela["electrons"][node],
                        "sentaurus_phin_V": sentaurus["phin"][node],
                        "vela_phin_V": vela["phin"][node],
                    })
                cases.append({
                    "case": case_name,
                    "device": device,
                    "drain_voltage_V": drain_voltage,
                    "state": state,
                    "selector_node_counts": {
                        name: len(nodes) for name, nodes in selectors.items()},
                    "integrated_srh": integrated_srh,
                    "metrics": selector_metrics,
                })
    detail_path = OUTPUT / "central_channel_field_comparison.csv"
    detail_path.parent.mkdir(parents=True, exist_ok=True)
    with detail_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(detail_rows[0]))
        writer.writeheader()
        writer.writerows(detail_rows)
    report = {
        "schema": "vela.simplemos.sdevice.m8_deep_off_field_compare.v1",
        "status": "complete",
        "convention": "ratios are Vela/Sentaurus; densities are compared in cm^-3",
        "cases": cases,
        "central_channel_detail_csv": str(detail_path.resolve()),
    }
    write_json(OUTPUT / "field_comparison_report.json", report)
    return report


def summarize(cases: list[dict[str, Any]]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for case in cases:
        root = Path(case["config"]).parent
        iv = csv_rows(root / "iv.csv")[-1]
        methods = [row for row in csv_rows(root / "terminal_current_methods.csv")
                   if row["contact"] == "drain"][-1]
        srh = csv_rows(root / "srh_balance.csv")[-1]
        rows.append({
            "case": case["case"],
            "device": case["device"],
            "drain_voltage_V": case["drain_voltage_V"],
            "gate_voltage_V": case["gate_voltage_V"],
            "current_A_per_um": float(iv["current_total_A_per_um"]),
            "electron_current_A_per_um": float(iv["current_electron_A_per_um"]),
            "hole_current_A_per_um": float(iv["current_hole_A_per_um"]),
            "sg_current_A_per_um": float(methods["I_sgflux_A_per_um"]),
            "residual_current_A_per_um": float(methods["I_residual_A_per_um"]),
            "qf_floor_current_A_per_um": float(
                methods["I_sgflux_with_qf_floor_A_per_um"]),
            "srh_net_current_A_per_um": float(srh["srh_net_current_A_per_um"]),
            "kcl_residual_A_per_um": float(
                srh["four_terminal_kcl_residual_A_per_um"]),
            "srh_numerical_status": srh["numerical_status"],
        })
    summary = {
        "schema": "vela.simplemos.sdevice.m8_deep_off_diagnostic_summary.v1",
        "status": "complete",
        "cases": rows,
    }
    write_json(OUTPUT / "vela_diagnostic_summary.json", summary)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--run-vela", action="store_true")
    parser.add_argument("--live-sentaurus", action="store_true")
    parser.add_argument("--export-sentaurus", action="store_true")
    parser.add_argument("--summarize", action="store_true")
    parser.add_argument("--analyze-fields", action="store_true")
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--scp-bin", default=executable("scp"))
    args = parser.parse_args()

    if not any((args.prepare, args.run_vela, args.live_sentaurus,
                args.export_sentaurus, args.summarize, args.analyze_fields)):
        parser.error("select at least one action")
    sentaurus_cases = prepare_sentaurus()
    vela_cases = prepare_vela()
    if args.run_vela:
        run_vela(vela_cases)
    if args.live_sentaurus:
        run_sentaurus(
            sentaurus_cases, args.ssh_target, args.ssh_bin, args.scp_bin)
    if args.export_sentaurus:
        export_sentaurus(sentaurus_cases)
    if args.summarize:
        print(json.dumps(summarize(vela_cases), indent=2))
    if args.analyze_fields:
        report = analyze_fields()
        print(json.dumps({"status": report["status"],
                          "cases": len(report["cases"])}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
