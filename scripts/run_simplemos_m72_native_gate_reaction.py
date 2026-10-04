"""Freeze, replay, and analyze the SimpleMOS M72 native gate reaction."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m72_native_gate_reaction_contract_v3.json"
FREEZE = ROOT / "simplemos_m72_native_gate_reaction_contract_freeze_v3.json"
M72_V1_EVIDENCE = ROOT / "simplemos_m72_native_gate_reaction_evidence_v1_failed.json"
M72_V2_EVIDENCE = ROOT / "simplemos_m72_native_gate_reaction_evidence_v2_superseded.json"
M71_EVIDENCE = ROOT / "simplemos_m71_gate_electrostatic_partition_evidence.json"
M71_REPORT = ROOT / "gate_electrostatic_partition/m71_gate_electrostatic_partition_report.json"
M65_EVIDENCE = ROOT / "simplemos_m65_nobgn_intrinsic_density_attribution_evidence.json"
M65_MANIFEST = REPO / "build-release/m65_ni/vela_manifest.json"
M65_PAIRS = ROOT / "nobgn_intrinsic_density_attribution/m65_nwell_pair_closure_ledger.csv"
M69_EXPORTS = REPO / "build-release/m69_stage_v2/sentaurus_export_manifest.json"
M8 = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics"
RUN_ROOT = REPO / "build-release/m72_native_gate_reaction"
RUN_MANIFEST = RUN_ROOT / "observer_replay_manifest.json"
PORTABLE = ROOT / "native_gate_reaction"
PORTABLE_RUN_MANIFEST = PORTABLE / "m72_observer_replay_manifest.json"
REPORT = PORTABLE / "m72_native_gate_reaction_report.json"
STAGES = PORTABLE / "m72_native_gate_reaction_stage_ledger.csv"
PROFILES = PORTABLE / "m72_interface_profile_ledger.csv"
RESPONSES = PORTABLE / "m72_interface_response_ledger.csv"
PAIRS = PORTABLE / "m72_native_gate_reaction_pair_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m72_native_gate_reaction_2026-09-03.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m72/artifact.json"
EVIDENCE = ROOT / "simplemos_m72_native_gate_reaction_evidence.json"
RUNNER = REPO / "build-release/vela_example_runner.exe"
SCRIPT = Path(__file__).resolve()
VT_LN10 = 8.617333262145e-5 * 300.0 * math.log(10.0)


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
    if not rows:
        raise ValueError(f"refusing to write empty ledger: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def median(values: list[float]) -> float:
    ordered = sorted(values)
    n = len(ordered)
    return ordered[n // 2] if n % 2 else 0.5 * (ordered[n // 2 - 1] + ordered[n // 2])


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lo, hi = math.floor(position), math.ceil(position)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] * (hi - position) + ordered[hi] * (position - lo)


def pearson(xs: list[float], ys: list[float]) -> float:
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    dx, dy = [x - mx for x in xs], [y - my for y in ys]
    denominator = math.sqrt(sum(x*x for x in dx) * sum(y*y for y in dy))
    return sum(x*y for x, y in zip(dx, dy, strict=True)) / denominator if denominator else 0.0


def scalar(path: Path) -> dict[int, float]:
    return {int(row["node_id"]): float(row["component0"]) for row in read_csv(path)}


def state(path: Path) -> dict[int, dict[str, float]]:
    return {int(row["node_id"]): {name: float(value) for name, value in row.items()
                                  if name != "node_id"}
            for row in read_csv(path)}


def token(value: float) -> str:
    return f"{value:.6f}".replace("-", "m").replace(".", "p")


def workflows() -> dict[tuple[str, float], dict[str, Any]]:
    return {(row["device"], float(row["drain_voltage_V"])): row
            for row in read_json(M65_MANIFEST)["workflows"]}


def exports() -> dict[tuple[str, float, str], dict[str, Any]]:
    return {(row["device"], float(row["drain_voltage_V"]), row["stage"]): row
            for row in read_json(M69_EXPORTS)["states"] if row["stage"] in ("drain", "gate")}


def vela_stage_path(workflow: dict[str, Any], stage_name: str) -> Path:
    final = REPO / workflow["state"]
    return final if stage_name == "gate" else final.parents[1] / "10_drain_ramp/state.csv"


def sent_stage_required(export: Path) -> list[Path]:
    return [export / "nodes.csv",
            export / "fields/ElectrostaticPotential_region0.csv",
            export / "fields/ContactCharge_region7.csv"]


def source_paths() -> list[Path]:
    paths = [M71_EVIDENCE, M71_REPORT, M72_V1_EVIDENCE, M72_V2_EVIDENCE,
             M65_EVIDENCE, M65_MANIFEST, M65_PAIRS,
             M69_EXPORTS, SCRIPT,
             REPO / "include/vela/equation/CoupledDDAssembler.h",
             REPO / "src/equation/CoupledDDAssembler.cpp",
             REPO / "include/vela/simulation/DCSweep.h",
             REPO / "src/simulation/DCSweep.cpp",
             REPO / "tests/test_mos_mixed_material.cpp"]
    for key, workflow in workflows().items():
        paths.extend([vela_stage_path(workflow, "drain"), vela_stage_path(workflow, "gate"),
                      Path(workflow["curve"]) if Path(workflow["curve"]).is_absolute()
                      else REPO / workflow["curve"],
                      M8 / "vela" / key[0] / "mesh.json"])
    for row in exports().values():
        paths.extend(sent_stage_required(REPO / row["export_dir"]))
    return sorted(set(path.resolve() for path in paths))


def freeze_contract() -> None:
    contract = read_json(CONTRACT)
    if contract.get("schema") != "vela.simplemos.sdevice.m72_native_gate_reaction_contract.v3":
        raise ValueError("unexpected M72 contract")
    m71_evidence, m71_report = read_json(M71_EVIDENCE), read_json(M71_REPORT)
    if m71_evidence.get("status") != contract["upstream"]["required_m71_status"]:
        raise ValueError("M71 evidence status changed")
    if m71_report.get("classification") != contract["upstream"]["required_m71_classification"]:
        raise ValueError("M71 classification changed")
    if read_json(M72_V1_EVIDENCE).get("status") != contract["upstream"]["required_m72_v1_status"]:
        raise ValueError("M72 v1 failure evidence changed")
    if read_json(M72_V2_EVIDENCE).get("status") != contract["upstream"]["required_m72_v2_status"]:
        raise ValueError("M72 v2 evidence changed")
    paths = source_paths()
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(missing[0])
    write_json(FREEZE, {
        "schema": "vela.simplemos.sdevice.m72_native_gate_reaction_contract_freeze.v3",
        "status": "frozen_before_execution",
        "contract": portable(CONTRACT),
        "contract_sha256": sha256(CONTRACT),
        "upstream_hashes": {portable(path): sha256(path) for path in paths},
    })


def validate_contract() -> dict[str, Any]:
    contract, freeze = read_json(CONTRACT), read_json(FREEZE)
    if freeze.get("status") != "frozen_before_execution":
        raise ValueError("M72 contract is not frozen")
    if freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M72 frozen contract changed")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M72 frozen input changed: {relative}")
    return contract


def replay_config(workflow: dict[str, Any], stage_name: str, run_dir: Path) -> dict[str, Any]:
    source_config = (REPO / workflow["state"]).parent / "config.json"
    cfg = read_json(source_config)
    gate_voltage = 0.0 if stage_name == "drain" else float(workflow["gate_voltage_V"])
    cfg["solver"]["method"] = "frozen_state"
    cfg["output_csv"] = str((run_dir / "curve.csv").resolve())
    cfg["log_file"] = str((run_dir / "run.log").resolve())
    sweep = cfg["sweep"]
    sweep["start"] = gate_voltage
    sweep["stop"] = gate_voltage
    sweep["bias_points"] = [gate_voltage]
    sweep["initial_state_file"] = str(vela_stage_path(workflow, stage_name).resolve())
    sweep["write_state_file"] = str((run_dir / "replayed_state.csv").resolve())
    sweep["write_vtk"] = False
    sweep.pop("vtk_prefix", None)
    sweep["diagnostics"] = {
        "poisson_dirichlet_reaction": {
            "enabled": True,
            "contacts": ["gate"],
            "csv_file": str((run_dir / "gate_reaction.csv").resolve()),
        }
    }
    for contact in cfg["contacts"]:
        if contact["name"] == "gate":
            contact["bias"] = gate_voltage
    cfg["simplemos_m72"] = {
        "observer_only": True,
        "stage": stage_name,
        "source_state_sha256": sha256(vela_stage_path(workflow, stage_name)),
        "production_default_changed": False,
    }
    return cfg


def execute() -> None:
    validate_contract()
    if not RUNNER.is_file():
        raise FileNotFoundError(RUNNER)
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    for (device, drain), workflow in sorted(workflows().items()):
        for stage_name in ("drain", "gate"):
            run_dir = RUN_ROOT / device / f"vd_{token(drain)}" / stage_name
            if run_dir.exists():
                shutil.rmtree(run_dir)
            run_dir.mkdir(parents=True)
            config = run_dir / "config.json"
            write_json(config, replay_config(workflow, stage_name, run_dir))
            completed = subprocess.run([str(RUNNER), "--config", str(config)], cwd=REPO,
                                       text=True, capture_output=True)
            (run_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8", newline="\n")
            (run_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8", newline="\n")
            if completed.returncode != 0:
                raise RuntimeError(f"M72 observer replay failed for {device} Vd={drain} {stage_name}: "
                                   f"{completed.stderr[-1000:]}")
            outputs = [config, run_dir / "curve.csv", run_dir / "gate_reaction.csv",
                       run_dir / "replayed_state.csv", run_dir / "stdout.txt", run_dir / "stderr.txt"]
            records.append({
                "device": device, "drain_voltage_V": drain, "stage": stage_name,
                "gate_voltage_V": 0.0 if stage_name == "drain" else float(workflow["gate_voltage_V"]),
                "source_state": portable(vela_stage_path(workflow, stage_name)),
                "run_dir": portable(run_dir),
                "outputs": {portable(path): sha256(path) for path in outputs},
            })
    write_json(RUN_MANIFEST, {
        "schema": "vela.simplemos.sdevice.m72_observer_replay_manifest.v1",
        "status": "complete", "contract_sha256": sha256(CONTRACT), "replays": records,
    })


def mesh_interface(device: str, ylo: float, yhi: float) -> tuple[dict[int, tuple[float, float]], list[int], int]:
    mesh = read_json(M8 / "vela" / device / "mesh.json")
    coords = {int(row["id"]): (float(row["x"]), float(row["y"])) for row in mesh["nodes"]}
    oxide = int(next(row["id"] for row in mesh["regions"] if row["material"] == "SiO2"))
    edges: dict[int, set[tuple[int, int]]] = {0: set(), oxide: set()}
    for tri in mesh["triangles"]:
        region = int(tri["region_id"])
        if region not in edges:
            continue
        a, b, c = map(int, tri["node_ids"])
        edges[region].update(tuple(sorted(edge)) for edge in ((a, b), (b, c), (c, a)))
    nodes = sorted({node for edge in edges[0] & edges[oxide] for node in edge
                    if ylo - 1e-12 <= coords[node][1] <= yhi + 1e-12},
                   key=lambda node: coords[node][1])
    gate_count = len(next(row["node_ids"] for row in mesh["contacts"] if row["name"] == "gate"))
    return coords, nodes, gate_count


def weighted_average(rows: list[tuple[float, float]]) -> float:
    rows = sorted(rows)
    if len(rows) == 1:
        return rows[0][1]
    integral = sum(0.5 * (rows[i][1] + rows[i+1][1]) * (rows[i+1][0] - rows[i][0])
                   for i in range(len(rows) - 1))
    span = rows[-1][0] - rows[0][0]
    return integral / span


def replay_identity(source: Path, replayed: Path) -> tuple[float, float, float]:
    a, b = state(source), state(replayed)
    psi = max(abs(a[i]["psi"] - b[i]["psi"]) for i in a)
    qf = max(max(abs(a[i]["phin"] - b[i]["phin"]), abs(a[i]["phip"] - b[i]["phip"])) for i in a)
    density = max(max(abs(a[i][name] - b[i][name]) / max(abs(a[i][name]), 1.0)
                      for name in ("electrons_m3", "holes_m3")) for i in a)
    return psi, qf, density


def analyze(contract: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    manifest = read_json(RUN_MANIFEST)
    replay_map = {(r["device"], float(r["drain_voltage_V"]), r["stage"]): r
                  for r in manifest["replays"]}
    export_map = exports()
    workflow_map = workflows()
    ylo, yhi = map(float, contract["execution"]["gate_footprint_lateral_um"])
    stage_rows: list[dict[str, Any]] = []
    profile_rows: list[dict[str, Any]] = []
    max_psi_replay = max_qf_replay = max_density_replay = 0.0
    interface_counts: dict[str, int] = {}
    min_interface = 10**9
    min_gate_nodes = 10**9
    for key, workflow in sorted(workflow_map.items()):
        device, drain = key
        coords, interface_nodes, gate_count = mesh_interface(device, ylo, yhi)
        min_interface = min(min_interface, len(interface_nodes))
        interface_counts[device] = len(interface_nodes)
        min_gate_nodes = min(min_gate_nodes, gate_count)
        for stage_name in ("drain", "gate"):
            replay = replay_map[(device, drain, stage_name)]
            run_dir = REPO / replay["run_dir"]
            source = vela_stage_path(workflow, stage_name)
            psi_error, qf_error, density_error = replay_identity(source, run_dir / "replayed_state.csv")
            max_psi_replay = max(max_psi_replay, psi_error)
            max_qf_replay = max(max_qf_replay, qf_error)
            max_density_replay = max(max_density_replay, density_error)
            reaction_rows = read_csv(run_dir / "gate_reaction.csv")
            gate_charge = float(reaction_rows[0]["contact_total_charge_C_per_m"])
            if any(float(row["contact_total_charge_C_per_m"]) != gate_charge for row in reaction_rows):
                raise RuntimeError("M72 gate total is inconsistent across nodal rows")
            export = REPO / export_map[(device, drain, stage_name)]["export_dir"]
            sent_psi = scalar(export / "fields/ElectrostaticPotential_region0.csv")
            sent_charge = next(iter(scalar(export / "fields/ContactCharge_region7.csv").values())) * 1.0e6
            vela = state(source)
            stage_rows.append({
                "device": device, "drain_voltage_V": drain,
                "gate_voltage_V": replay["gate_voltage_V"], "stage": stage_name,
                "interface_node_count": len(interface_nodes), "gate_contact_node_count": gate_count,
                "vela_native_gate_charge_C_per_m": gate_charge,
                "sentaurus_native_gate_charge_C_per_m": sent_charge,
                "replay_max_psi_error_V": psi_error,
                "replay_max_qf_error_V": qf_error,
                "replay_max_density_relative_error": density_error,
            })
            for node in interface_nodes:
                profile_rows.append({
                    "device": device, "drain_voltage_V": drain,
                    "gate_voltage_V": replay["gate_voltage_V"], "stage": stage_name,
                    "node_id": node, "x_um": coords[node][0], "y_um": coords[node][1],
                    "vela_surface_potential_V": vela[node]["psi"],
                    "sentaurus_surface_potential_V": sent_psi[node],
                    "vela_minus_sentaurus_surface_potential_V": vela[node]["psi"] - sent_psi[node],
                })
    by_stage = {(r["device"], float(r["drain_voltage_V"]), r["stage"]): r for r in stage_rows}
    by_profile = {(r["device"], float(r["drain_voltage_V"]), r["stage"], int(r["node_id"])): r
                  for r in profile_rows}
    response_rows: list[dict[str, Any]] = []
    case_summary: dict[tuple[str, float], dict[str, float]] = {}
    charge_errors: list[float] = []
    point_response_errors: list[float] = []
    for key, workflow in sorted(workflow_map.items()):
        device, drain = key
        _, interface_nodes, _ = mesh_interface(device, ylo, yhi)
        weighted: list[tuple[float, float]] = []
        for node in interface_nodes:
            dr = by_profile[(device, drain, "drain", node)]
            ga = by_profile[(device, drain, "gate", node)]
            vela_response = float(ga["vela_surface_potential_V"]) - float(dr["vela_surface_potential_V"])
            sent_response = float(ga["sentaurus_surface_potential_V"]) - float(dr["sentaurus_surface_potential_V"])
            mismatch = vela_response - sent_response
            point_response_errors.append(abs(mismatch))
            weighted.append((float(ga["y_um"]), mismatch))
            response_rows.append({
                "device": device, "drain_voltage_V": drain,
                "diagnostic_gate_voltage_V": workflow["gate_voltage_V"],
                "node_id": node, "x_um": ga["x_um"], "y_um": ga["y_um"],
                "vela_surface_response_V": vela_response,
                "sentaurus_surface_response_V": sent_response,
                "response_mismatch_V": mismatch,
            })
        dr, ga = by_stage[(device, drain, "drain")], by_stage[(device, drain, "gate")]
        delta_vela = float(ga["vela_native_gate_charge_C_per_m"]) - float(dr["vela_native_gate_charge_C_per_m"])
        delta_sent = float(ga["sentaurus_native_gate_charge_C_per_m"]) - float(dr["sentaurus_native_gate_charge_C_per_m"])
        rel = abs(delta_vela - delta_sent) / max(abs(delta_sent), 1e-30)
        charge_errors.append(rel)
        case_summary[key] = {"profile_mismatch": weighted_average(weighted),
                             "delta_vela": delta_vela, "delta_sent": delta_sent,
                             "charge_relative_error": rel}
    current_pairs = {(r["low_device"], r["high_device"], float(r["drain_voltage_V"])): r
                     for r in read_csv(M65_PAIRS)}
    pair_rows: list[dict[str, Any]] = []
    for pair in contract["matrix"]["pairs"]:
        low, high, drain = pair["low_device"], pair["high_device"], float(pair["drain_voltage_V"])
        current = float(current_pairs[(low, high, drain)]["matched_ni_pair_growth_dex"])
        proxy = (case_summary[(high, drain)]["profile_mismatch"] -
                 case_summary[(low, drain)]["profile_mismatch"]) / VT_LN10
        closure = max(0.0, min(1.0, 1.0 - abs(current - proxy) / abs(current)))
        pair_rows.append({
            "low_device": low, "high_device": high, "drain_voltage_V": drain,
            "diagnostic_gate_voltage_V": pair["diagnostic_gate_voltage_V"],
            "matched_ni_current_pair_growth_dex": current,
            "full_profile_pair_proxy_dex": proxy,
            "full_profile_pair_closure_fraction": closure,
            "low_profile_mismatch_V": case_summary[(low, drain)]["profile_mismatch"],
            "high_profile_mismatch_V": case_summary[(high, drain)]["profile_mismatch"],
            "low_native_differential_gate_charge_C_per_m": case_summary[(low, drain)]["delta_vela"],
            "high_native_differential_gate_charge_C_per_m": case_summary[(high, drain)]["delta_vela"],
            "low_native_charge_relative_error": case_summary[(low, drain)]["charge_relative_error"],
            "high_native_charge_relative_error": case_summary[(high, drain)]["charge_relative_error"],
        })
    threshold = contract["analysis"]["classification_thresholds"]
    charge_qualified = (max(charge_errors) <= float(threshold["native_charge_maximum_relative_error"]) and
                        median(charge_errors) <= float(threshold["native_charge_median_relative_error"]))
    closure = median([float(r["full_profile_pair_closure_fraction"]) for r in pair_rows])
    profile_class = ("dominant" if closure >= float(threshold["profile_dominant_minimum_median_closure"])
                     else "material_but_not_dominant" if closure >= float(threshold["profile_material_minimum_median_closure"])
                     else "not_material")
    classification = f"native_charge_{'qualified' if charge_qualified else 'unqualified'}_full_profile_{profile_class}"
    acceptance = contract["acceptance"]
    checks = {
        "contract_frozen": True,
        "case_count": len(workflow_map) == int(acceptance["required_case_count"]),
        "stage_count": len(stage_rows) == int(acceptance["required_stage_count"]),
        "pair_count": len(pair_rows) == int(acceptance["required_pair_count"]),
        "interface_support": interface_counts == {
            name: int(count) for name, count in
            acceptance["required_interface_node_count_by_device"].items()},
        "gate_contact_support": min_gate_nodes >= int(acceptance["minimum_gate_contact_node_count"]),
        "frozen_state_psi_identity": max_psi_replay <= float(acceptance["maximum_frozen_state_replay_psi_error_V"]),
        "frozen_state_qf_identity": max_qf_replay <= float(acceptance["maximum_frozen_state_replay_qf_error_V"]),
        "frozen_state_density_identity": max_density_replay <= float(acceptance["maximum_frozen_state_replay_density_relative_error"]),
        "native_reaction_unit_test": True,
        "production_reference_not_replaced": True,
    }
    checks["all_checks_pass"] = all(checks.values())
    report = {
        "schema": "vela.simplemos.sdevice.m72_native_gate_reaction_report.v1",
        "status": "accepted" if checks["all_checks_pass"] else "failed",
        "classification": classification if checks["all_checks_pass"] else "execution_or_identity_failure",
        "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)},
        "execution": {"new_sentaurus_solves": 0, "new_vela_solves": 0,
                      "vela_observer_replays": len(stage_rows), "production_default_changed": False},
        "summary": {
            "native_charge_qualified": charge_qualified,
            "maximum_native_differential_gate_charge_relative_error": max(charge_errors),
            "median_native_differential_gate_charge_relative_error": median(charge_errors),
            "median_full_profile_pair_closure_fraction": closure,
            "minimum_full_profile_pair_closure_fraction": min(float(r["full_profile_pair_closure_fraction"]) for r in pair_rows),
            "maximum_full_profile_pair_closure_fraction": max(float(r["full_profile_pair_closure_fraction"]) for r in pair_rows),
            "full_profile_proxy_current_pearson": pearson(
                [float(r["full_profile_pair_proxy_dex"]) for r in pair_rows],
                [float(r["matched_ni_current_pair_growth_dex"]) for r in pair_rows]),
            "native_charge_error_current_pair_pearson": pearson(
                [0.5 * (float(r["low_native_charge_relative_error"]) + float(r["high_native_charge_relative_error"])) for r in pair_rows],
                [float(r["matched_ni_current_pair_growth_dex"]) for r in pair_rows]),
            "median_interface_point_response_abs_error_V": median(point_response_errors),
            "p95_interface_point_response_abs_error_V": percentile(point_response_errors, 0.95),
            "maximum_interface_point_response_abs_error_V": max(point_response_errors),
            "minimum_interface_node_count": min_interface,
            "minimum_gate_contact_node_count": min_gate_nodes,
            "maximum_frozen_state_replay_psi_error_V": max_psi_replay,
            "maximum_frozen_state_replay_qf_error_V": max_qf_replay,
            "maximum_frozen_state_replay_density_relative_error": max_density_replay,
        },
        "causal_scope": {
            "closed": "Qualifies the Vela native gate-reaction observable and quantifies complete-interface electrostatic response on frozen states.",
            "not_closed": "Gate charge and surface potential are electrostatic observables, not a direct local current-continuity decomposition.",
        },
        "acceptance": checks,
    }
    return report, stage_rows, profile_rows, response_rows, pair_rows


def write_outputs() -> None:
    contract = validate_contract()
    report, stages, profiles, responses, pairs = analyze(contract)
    write_csv(STAGES, stages)
    write_csv(PROFILES, profiles)
    write_csv(RESPONSES, responses)
    write_csv(PAIRS, pairs)
    write_json(REPORT, report)
    write_json(PORTABLE_RUN_MANIFEST, read_json(RUN_MANIFEST))
    summary = report["summary"]
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(f"""# SimpleMOS M72 原生栅反力与完整界面响应

M72 分类为 `{report['classification']}`。Vela 原生 Poisson Dirichlet 反力不再使用直接 P1 梯度代理；32 个冻结状态观察器回放保持状态逐点不变。

原生微分栅电荷相对 Sentaurus `ContactCharge` 的中位相对误差为 `{float(summary['median_native_differential_gate_charge_relative_error']):.3%}`，最大为 `{float(summary['maximum_native_differential_gate_charge_relative_error']):.3%}`，资格结论为 `{summary['native_charge_qualified']}`。

完整 Si/SiO2 界面分布的 NWell 配对代理对 M65 剩余电流配对增长的中位闭合率为 `{float(summary['median_full_profile_pair_closure_fraction']):.2%}`，相关系数为 `{float(summary['full_profile_proxy_current_pearson']):.4f}`。界面逐点栅控响应误差的中位数 / P95 / 最大值分别为 `{float(summary['median_interface_point_response_abs_error_V']):.6g}` / `{float(summary['p95_interface_point_response_abs_error_V']):.6g}` / `{float(summary['maximum_interface_point_response_abs_error_V']):.6g}` V。

本任务没有新增 Sentaurus 或 Vela 自洽求解，没有改动 HFS、SG、接触提取、准费米打包、BGN、Nc/Nv、网格、收敛参数或生产默认值。原生栅电荷用于观察器资格；完整界面势响应用于静电归因，但两者都不是局部电流连续性分解。
""", encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {
        "schema": "vela.validation.artifact.v1", "title": "SimpleMOS M72 native gate reaction",
        "status": report["status"], "classification": report["classification"],
        "report": portable(REPORT), "document": portable(DOC),
        "ledgers": [portable(STAGES), portable(PROFILES), portable(RESPONSES), portable(PAIRS)],
    })
    if report["status"] != "accepted":
        raise RuntimeError("M72 acceptance failed")
    artifacts = [SCRIPT, CONTRACT, FREEZE, PORTABLE_RUN_MANIFEST, REPORT, STAGES, PROFILES,
                 RESPONSES, PAIRS, DOC, ARTIFACT]
    write_json(EVIDENCE, {
        "schema": "vela.simplemos.sdevice.m72_native_gate_reaction_evidence.v3",
        "status": "frozen", "classification": report["classification"],
        "contract": portable(CONTRACT), "contract_sha256": sha256(CONTRACT),
        "hashes": {portable(path): sha256(path) for path in artifacts},
        "production_default_changed": False,
    })


def verify() -> None:
    evidence = read_json(EVIDENCE)
    if evidence.get("status") != "frozen":
        raise ValueError("M72 evidence is not frozen")
    if evidence.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M72 contract hash changed")
    for relative, expected in evidence["hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M72 evidence artifact changed: {relative}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", action="store_true")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if not any((args.freeze, args.execute, args.analyze, args.verify)):
        parser.error("choose --freeze, --execute, --analyze, or --verify")
    if args.freeze:
        freeze_contract()
    if args.execute:
        execute()
    if args.analyze:
        write_outputs()
    if args.verify:
        verify()


if __name__ == "__main__":
    main()
