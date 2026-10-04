"""Freeze and analyze the SimpleMOS M71 gate electrostatic partition."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m71_gate_electrostatic_partition_contract_v2.json"
FREEZE = ROOT / "simplemos_m71_gate_electrostatic_partition_contract_freeze_v2.json"
M71_V1_EVIDENCE = ROOT / "simplemos_m71_gate_electrostatic_partition_evidence_v1_failed.json"
M69_EVIDENCE = ROOT / "simplemos_m69_stage_residual_localization_evidence.json"
M70_EVIDENCE = ROOT / "simplemos_m70_transport_support_partition_evidence.json"
M65_EVIDENCE = ROOT / "simplemos_m65_nobgn_intrinsic_density_attribution_evidence.json"
M69_EXPORTS = REPO / "build-release/m69_stage_v2/sentaurus_export_manifest.json"
M65_VELA = REPO / "build-release/m65_ni/vela_manifest.json"
M65_PAIRS = ROOT / "nobgn_intrinsic_density_attribution/m65_nwell_pair_closure_ledger.csv"
M67_NODES = ROOT / "residual_barrier_partition/m67_fixed_node_state_ledger.csv"
M8 = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics"
PORTABLE = ROOT / "gate_electrostatic_partition"
REPORT = PORTABLE / "m71_gate_electrostatic_partition_report.json"
STATES = PORTABLE / "m71_gate_electrostatic_state_ledger.csv"
PAIRS = PORTABLE / "m71_gate_electrostatic_pair_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m71_gate_electrostatic_partition_2026-09-03.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m71/artifact.json"
EVIDENCE = ROOT / "simplemos_m71_gate_electrostatic_partition_evidence.json"
SCRIPT = Path(__file__).resolve()
EPS0 = 8.8541878128e-12
Q = 1.602176634e-19
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
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def median(values: list[float]) -> float:
    values = sorted(values); n = len(values)
    return values[n // 2] if n % 2 else 0.5 * (values[n // 2 - 1] + values[n // 2])


def pearson(xs: list[float], ys: list[float]) -> float:
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    dx, dy = [x - mx for x in xs], [y - my for y in ys]
    den = math.sqrt(sum(x*x for x in dx) * sum(y*y for y in dy))
    return sum(x*y for x, y in zip(dx, dy, strict=True)) / den if den else 0.0


def scalar(path: Path) -> dict[int, float]:
    return {int(row["node_id"]): float(row["component0"]) for row in read_csv(path)}


def state(path: Path) -> dict[str, dict[int, float]]:
    rows = read_csv(path)
    return {
        "psi": {int(r["node_id"]): float(r["psi"]) for r in rows},
        "n": {int(r["node_id"]): float(r["electrons_m3"]) / 1e6 for r in rows},
        "p": {int(r["node_id"]): float(r["holes_m3"]) / 1e6 for r in rows},
    }


def vela_stage_path(workflow: dict[str, Any], stage_name: str) -> Path:
    final = REPO / workflow["state"]
    return final if stage_name == "gate" else final.parents[1] / "10_drain_ramp/state.csv"


def dielectric_regions(device: str) -> list[tuple[int, str]]:
    mesh = read_json(M8 / "vela" / device / "mesh.json")
    return [(int(r["id"]), r["material"]) for r in mesh["regions"]
            if r["material"] in ("SiO2", "Nitride")]


def required_state_files(export: Path, device: str) -> list[Path]:
    fields = export / "fields"
    paths = [export / "nodes.csv", export / "elements.csv", export / "contacts.csv",
            fields / "ElectrostaticPotential_region0.csv",
            fields / "eDensity_region0.csv", fields / "hDensity_region0.csv",
            fields / "ContactExternalVoltage_region7.csv",
            fields / "ContactCharge_region7.csv"]
    paths.extend(fields / f"ElectrostaticPotential_region{region}.csv"
                 for region, _ in dielectric_regions(device))
    return paths


def source_paths() -> list[Path]:
    paths = [M69_EVIDENCE, M70_EVIDENCE, M71_V1_EVIDENCE, M65_EVIDENCE, M69_EXPORTS, M65_VELA,
             M65_PAIRS, M67_NODES]
    workflows = {(row["device"], float(row["drain_voltage_V"])): row
                 for row in read_json(M65_VELA)["workflows"]}
    paths.append(REPO / read_json(M65_VELA)["materials"])
    for row in read_json(M69_EXPORTS)["states"]:
        if row["stage"] not in ("drain", "gate"):
            continue
        paths.extend(required_state_files(REPO / row["export_dir"], row["device"]))
    for key, workflow in workflows.items():
        paths.extend([vela_stage_path(workflow, "drain"), vela_stage_path(workflow, "gate"),
                      M8 / "vela" / key[0] / "mesh.json",
                      M8 / "neutral" / key[0] / "doping.csv"])
    return sorted(set(paths))


def freeze_contract() -> None:
    contract = read_json(CONTRACT)
    if contract.get("schema") != "vela.simplemos.sdevice.m71_gate_electrostatic_partition_contract.v2":
        raise ValueError("unexpected M71 contract")
    for evidence, name in ((M69_EVIDENCE, "m69"), (M70_EVIDENCE, "m70")):
        upstream = read_json(evidence)
        if upstream.get("status") != contract["upstream"][f"required_{name}_status"]:
            raise ValueError(f"{name.upper()} status changed")
        if upstream.get("classification") != contract["upstream"][f"required_{name}_classification"]:
            raise ValueError(f"{name.upper()} classification changed")
    if read_json(M71_V1_EVIDENCE).get("status") != contract["upstream"]["required_m71_v1_status"]:
        raise ValueError("M71 v1 failure evidence changed")
    paths = source_paths(); missing = [p for p in paths if not p.is_file()]
    if missing: raise FileNotFoundError(missing[0])
    write_json(FREEZE, {
        "schema": "vela.simplemos.sdevice.m71_gate_electrostatic_partition_contract_freeze.v2",
        "status": "frozen_before_execution", "contract": portable(CONTRACT),
        "contract_sha256": sha256(CONTRACT),
        "upstream_hashes": {portable(path): sha256(path) for path in paths},
    })


def validate_contract() -> dict[str, Any]:
    contract, freeze = read_json(CONTRACT), read_json(FREEZE)
    if freeze.get("status") != "frozen_before_execution" or freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M71 contract is not frozen")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M71 frozen input changed: {relative}")
    return contract


def mesh_data(device: str) -> tuple[dict[int, tuple[float, float]], list[dict[str, Any]], set[int], int, dict[int, float], dict[int, float]]:
    mesh = read_json(M8 / "vela" / device / "mesh.json")
    coords = {int(n["id"]): (float(n["x"]), float(n["y"])) for n in mesh["nodes"]}
    gate = set(next(c["node_ids"] for c in mesh["contacts"] if c["name"] == "gate"))
    oxide = int(next(r["id"] for r in mesh["regions"] if r["material"] == "SiO2"))
    doping = read_csv(M8 / "neutral" / device / "doping.csv")
    donors = {int(r["node_id"]): float(r["donors_cm3"]) for r in doping}
    acceptors = {int(r["node_id"]): float(r["acceptors_cm3"]) for r in doping}
    return coords, mesh["triangles"], gate, oxide, donors, acceptors


def region_edges(triangles: list[dict[str, Any]], region_id: int) -> tuple[dict[tuple[int, int], list[int]], set[tuple[int, int]]]:
    owners: dict[tuple[int, int], list[int]] = {}
    for tri in triangles:
        if int(tri["region_id"]) != region_id: continue
        a, b, c = map(int, tri["node_ids"])
        for edge in ((a, b), (b, c), (c, a)):
            owners.setdefault(tuple(sorted(edge)), []).append(int(tri["id"]))
    return owners, set(owners)


def triangle_gradient(points: list[tuple[float, float]], values: list[float]) -> tuple[float, float]:
    (x0, y0), (x1, y1), (x2, y2) = points
    det = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0)
    gx = ((values[1]-values[0])*(y2-y0) - (values[2]-values[0])*(y1-y0)) / det
    gy = ((x1-x0)*(values[2]-values[0]) - (x2-x0)*(values[1]-values[0])) / det
    return gx, gy


def reconstructed_gate_charge(coords: dict[int, tuple[float, float]], triangles: list[dict[str, Any]],
                              gate_nodes: set[int], oxide: int, psi: dict[int, float], eps_r: float) -> tuple[float, int]:
    owners, _ = region_edges(triangles, oxide)
    by_id = {int(t["id"]): t for t in triangles}
    boundary = [(edge, ids[0]) for edge, ids in owners.items()
                if len(ids) == 1 and edge[0] in gate_nodes and edge[1] in gate_nodes]
    charge = 0.0
    for (a, b), tri_id in boundary:
        nodes = list(map(int, by_id[tri_id]["node_ids"]))
        grad = triangle_gradient([coords[n] for n in nodes], [psi[n] for n in nodes])
        midpoint = ((coords[a][0]+coords[b][0])/2, (coords[a][1]+coords[b][1])/2)
        centroid = (sum(coords[n][0] for n in nodes)/3, sum(coords[n][1] for n in nodes)/3)
        tx, ty = coords[b][0]-coords[a][0], coords[b][1]-coords[a][1]
        edge_length_um = math.hypot(tx, ty)
        nx, ny = ty/edge_length_um, -tx/edge_length_um
        if nx*(midpoint[0]-centroid[0]) + ny*(midpoint[1]-centroid[1]) < 0:
            nx, ny = -nx, -ny
        length_m = edge_length_um * 1e-6
        grad_dot_n_V_per_m = (grad[0]*nx + grad[1]*ny) * 1e6
        charge += EPS0 * eps_r * grad_dot_n_V_per_m * length_m
    return charge, len(boundary)


def footprint_charge(coords: dict[int, tuple[float, float]], triangles: list[dict[str, Any]],
                     n: dict[int, float], p: dict[int, float], donors: dict[int, float],
                     acceptors: dict[int, float], ylo: float, yhi: float) -> dict[str, float]:
    result = {"electron": 0.0, "hole": 0.0, "dopant": 0.0}
    for tri in triangles:
        if int(tri["region_id"]) != 0: continue
        nodes = list(map(int, tri["node_ids"])); yc = sum(coords[i][1] for i in nodes) / 3
        if yc < ylo or yc > yhi: continue
        (x0,y0),(x1,y1),(x2,y2) = [coords[i] for i in nodes]
        area_m2 = abs((x1-x0)*(y2-y0)-(x2-x0)*(y1-y0)) * 0.5e-12
        avg = lambda values: sum(values[i] for i in nodes) / 3 * 1e6
        result["electron"] += -Q * avg(n) * area_m2
        result["hole"] += Q * avg(p) * area_m2
        result["dopant"] += Q * (avg(donors)-avg(acceptors)) * area_m2
    result["net"] = result["electron"] + result["hole"] + result["dopant"]
    return result


def analyze(contract: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    workflows = {(r["device"], float(r["drain_voltage_V"])): r for r in read_json(M65_VELA)["workflows"]}
    exports = {(r["device"], float(r["drain_voltage_V"]), r["stage"]): r
               for r in read_json(M69_EXPORTS)["states"] if r["stage"] in ("drain", "gate")}
    controls = {(r["device"], float(r["drain_voltage_V"])): int(r["sentaurus_control_node"])
                for r in read_csv(M67_NODES)}
    ylo, yhi = map(float, contract["execution"]["gate_footprint_lateral_um"])
    eps_by_material = {name: float(value) for name, value in
                       contract["execution"]["dielectric_relative_permittivity"].items()}
    state_rows: list[dict[str, Any]] = []
    coord_errors: list[float] = []
    for key, workflow in workflows.items():
        device, drain = key; coords, triangles, gate_nodes, oxide, donors, acceptors = mesh_data(device)
        dielectrics = dielectric_regions(device)
        _, si_edges = region_edges(triangles, 0); _, ox_edges = region_edges(triangles, oxide)
        interface_nodes = {node for edge in si_edges & ox_edges for node in edge if ylo-1e-12 <= coords[node][1] <= yhi+1e-12}
        control = controls[key]
        if control not in interface_nodes: raise RuntimeError(f"M71 control node is not on Si/oxide interface: {key}")
        partner = min(gate_nodes, key=lambda node: (abs(coords[node][1]-coords[control][1]), -coords[node][0]))
        for stage_name in ("drain", "gate"):
            export = REPO / exports[(device, drain, stage_name)]["export_dir"]
            exported_coords = {int(r["id"]): (float(r["x_um"]), float(r["y_um"])) for r in read_csv(export / "nodes.csv")}
            coord_error = max(max(abs(exported_coords[i][j]-coords[i][j]) for j in (0,1)) for i in coords)
            coord_errors.append(coord_error)
            sent_psi0 = scalar(export / "fields/ElectrostaticPotential_region0.csv")
            sent_dielectric_psi = {
                region: scalar(export / f"fields/ElectrostaticPotential_region{region}.csv")
                for region, _ in dielectrics}
            sent_psi1 = sent_dielectric_psi[oxide]
            sent_n = scalar(export / "fields/eDensity_region0.csv"); sent_p = scalar(export / "fields/hDensity_region0.csv")
            vela = state(vela_stage_path(workflow, stage_name))
            sent_oxide_q, oxide_edge_count = reconstructed_gate_charge(
                coords, triangles, gate_nodes, oxide, sent_psi1, eps_by_material["SiO2"])
            vela_oxide_q, vela_oxide_edge_count = reconstructed_gate_charge(
                coords, triangles, gate_nodes, oxide, vela["psi"], eps_by_material["SiO2"])
            sent_q = 0.0; vela_q = 0.0; edge_count = 0; vela_edge_count = 0
            for region, material in dielectrics:
                q_sent, n_sent = reconstructed_gate_charge(
                    coords, triangles, gate_nodes, region, sent_dielectric_psi[region], eps_by_material[material])
                q_vela, n_vela = reconstructed_gate_charge(
                    coords, triangles, gate_nodes, region, vela["psi"], eps_by_material[material])
                sent_q += q_sent; vela_q += q_vela; edge_count += n_sent; vela_edge_count += n_vela
            native_q = next(iter(scalar(export / "fields/ContactCharge_region7.csv").values())) * 1e6
            external_v = next(iter(scalar(export / "fields/ContactExternalVoltage_region7.csv").values()))
            sent_charge = footprint_charge(coords, triangles, sent_n, sent_p, donors, acceptors, ylo, yhi)
            vela_charge = footprint_charge(coords, triangles, vela["n"], vela["p"], donors, acceptors, ylo, yhi)
            state_rows.append({
                "device": device, "drain_voltage_V": drain, "gate_voltage_V": float(workflow["gate_voltage_V"]),
                "stage": stage_name, "control_node": control, "gate_partner_node": partner,
                "maximum_coordinate_error_um": coord_error,
                "oxide_gate_boundary_edge_count": oxide_edge_count,
                "full_dielectric_gate_boundary_edge_count": edge_count,
                "sentaurus_external_gate_voltage_V": external_v,
                "sentaurus_surface_potential_V": sent_psi0[control], "vela_surface_potential_V": vela["psi"][control],
                "sentaurus_gate_boundary_potential_V": sent_psi1[partner], "vela_gate_boundary_potential_V": vela["psi"][partner],
                "sentaurus_oxide_drop_V": sent_psi1[partner]-sent_psi0[control],
                "vela_oxide_drop_V": vela["psi"][partner]-vela["psi"][control],
                "sentaurus_oxide_only_gate_charge_C_per_m": sent_oxide_q,
                "vela_oxide_only_gate_charge_C_per_m": vela_oxide_q,
                "sentaurus_reconstructed_gate_charge_C_per_m": sent_q,
                "sentaurus_native_gate_charge_C_per_m": native_q,
                "vela_reconstructed_gate_charge_C_per_m": vela_q,
                "sentaurus_electron_charge_C_per_m": sent_charge["electron"],
                "sentaurus_hole_charge_C_per_m": sent_charge["hole"],
                "sentaurus_dopant_charge_C_per_m": sent_charge["dopant"],
                "sentaurus_net_semiconductor_charge_C_per_m": sent_charge["net"],
                "vela_electron_charge_C_per_m": vela_charge["electron"],
                "vela_hole_charge_C_per_m": vela_charge["hole"],
                "vela_dopant_charge_C_per_m": vela_charge["dopant"],
                "vela_net_semiconductor_charge_C_per_m": vela_charge["net"],
                "edge_count_identity": int(edge_count == vela_edge_count and
                                           oxide_edge_count == vela_oxide_edge_count),
            })
    by_state = {(r["device"], float(r["drain_voltage_V"]), r["stage"]): r for r in state_rows}
    current_pairs = {(r["low_device"], r["high_device"], float(r["drain_voltage_V"])): r for r in read_csv(M65_PAIRS)}
    pair_rows: list[dict[str, Any]] = []; partition_errors: list[float] = []
    for pair in contract["matrix"]["pairs"]:
        low, high, drain = pair["low_device"], pair["high_device"], float(pair["drain_voltage_V"])
        response: dict[str, dict[str, float]] = {}
        for device in (low, high):
            dr, ga = by_state[(device, drain, "drain")], by_state[(device, drain, "gate")]
            surface = {s: float(ga[f"{s}_surface_potential_V"])-float(dr[f"{s}_surface_potential_V"]) for s in ("sentaurus","vela")}
            boundary = {s: float(ga[f"{s}_gate_boundary_potential_V"])-float(dr[f"{s}_gate_boundary_potential_V"]) for s in ("sentaurus","vela")}
            oxide = {s: float(ga[f"{s}_oxide_drop_V"])-float(dr[f"{s}_oxide_drop_V"]) for s in ("sentaurus","vela")}
            for s in ("sentaurus","vela"): partition_errors.append(abs(boundary[s]-surface[s]-oxide[s]))
            response[device] = {
                "surface_mismatch": surface["vela"]-surface["sentaurus"],
                "oxide_mismatch": oxide["vela"]-oxide["sentaurus"],
                "boundary_mismatch": boundary["vela"]-boundary["sentaurus"],
                "gate_charge_mismatch": (float(ga["vela_reconstructed_gate_charge_C_per_m"])-float(dr["vela_reconstructed_gate_charge_C_per_m"])) - (float(ga["sentaurus_reconstructed_gate_charge_C_per_m"])-float(dr["sentaurus_reconstructed_gate_charge_C_per_m"])),
                "electron_charge_mismatch": (float(ga["vela_electron_charge_C_per_m"])-float(dr["vela_electron_charge_C_per_m"])) - (float(ga["sentaurus_electron_charge_C_per_m"])-float(dr["sentaurus_electron_charge_C_per_m"])),
                "hole_charge_mismatch": (float(ga["vela_hole_charge_C_per_m"])-float(dr["vela_hole_charge_C_per_m"])) - (float(ga["sentaurus_hole_charge_C_per_m"])-float(dr["sentaurus_hole_charge_C_per_m"])),
                "net_charge_mismatch": (float(ga["vela_net_semiconductor_charge_C_per_m"])-float(dr["vela_net_semiconductor_charge_C_per_m"])) - (float(ga["sentaurus_net_semiconductor_charge_C_per_m"])-float(dr["sentaurus_net_semiconductor_charge_C_per_m"])),
            }
        surface_proxy = (response[high]["surface_mismatch"]-response[low]["surface_mismatch"]) / VT_LN10
        current_growth = float(current_pairs[(low, high, drain)]["matched_ni_pair_growth_dex"])
        closure = max(0.0, min(1.0, 1.0-abs(current_growth-surface_proxy)/abs(current_growth)))
        pair_rows.append({
            "low_device": low, "high_device": high, "drain_voltage_V": drain,
            "diagnostic_gate_voltage_V": pair["diagnostic_gate_voltage_V"],
            "matched_ni_current_pair_growth_dex": current_growth,
            "surface_response_pair_proxy_dex": surface_proxy,
            "oxide_response_pair_proxy_dex": -(response[high]["oxide_mismatch"]-response[low]["oxide_mismatch"]) / VT_LN10,
            "gate_boundary_response_pair_proxy_dex": (response[high]["boundary_mismatch"]-response[low]["boundary_mismatch"]) / VT_LN10,
            "surface_pair_closure_fraction": closure,
            "gate_charge_response_mismatch_pair_C_per_m": response[high]["gate_charge_mismatch"]-response[low]["gate_charge_mismatch"],
            "electron_charge_response_mismatch_pair_C_per_m": response[high]["electron_charge_mismatch"]-response[low]["electron_charge_mismatch"],
            "hole_charge_response_mismatch_pair_C_per_m": response[high]["hole_charge_mismatch"]-response[low]["hole_charge_mismatch"],
            "net_semiconductor_charge_response_mismatch_pair_C_per_m": response[high]["net_charge_mismatch"]-response[low]["net_charge_mismatch"],
        })
    closures = [float(r["surface_pair_closure_fraction"]) for r in pair_rows]
    median_closure = median(closures); thresholds = contract["analysis"]["classification_thresholds"]
    charge_errors = [abs(float(r["sentaurus_reconstructed_gate_charge_C_per_m"])-float(r["sentaurus_native_gate_charge_C_per_m"])) / max(abs(float(r["sentaurus_native_gate_charge_C_per_m"])), 1e-30) for r in state_rows]
    charge_qualified = max(charge_errors) <= float(thresholds["gate_charge_proxy_maximum_native_relative_error"])
    charge_suffix = "charge_proxy_qualified" if charge_qualified else "charge_proxy_unqualified"
    if median_closure >= float(thresholds["surface_partition_dominant_minimum_median_closure"]):
        classification = f"gate_surface_partition_dominant_{charge_suffix}"
    elif median_closure >= float(thresholds["surface_partition_material_minimum_median_closure"]):
        classification = f"gate_surface_partition_material_but_not_dominant_{charge_suffix}"
    else:
        classification = f"gate_surface_partition_not_material_{charge_suffix}"
    gate_voltage_errors = [abs(float(r["sentaurus_external_gate_voltage_V"])-(float(r["gate_voltage_V"]) if r["stage"] == "gate" else 0.0)) for r in state_rows]
    checks = {
        "contract_frozen": read_json(FREEZE)["contract_sha256"] == sha256(CONTRACT),
        "case_count": len(workflows) == int(contract["acceptance"]["required_case_count"]),
        "stage_count": len(state_rows) == int(contract["acceptance"]["required_stage_count"]),
        "pair_count": len(pair_rows) == int(contract["acceptance"]["required_pair_count"]),
        "coordinate_identity": max(coord_errors) <= float(contract["acceptance"]["maximum_coordinate_error_um"]),
        "gate_partition_identity": max(partition_errors) <= float(contract["acceptance"]["maximum_gate_partition_identity_error_V"]),
        "sentaurus_gate_voltage_identity": max(gate_voltage_errors) <= float(contract["acceptance"]["maximum_sentaurus_gate_voltage_error_V"]),
        "full_dielectric_gate_boundary_support": min(int(r["full_dielectric_gate_boundary_edge_count"]) for r in state_rows) >= int(contract["acceptance"]["minimum_full_dielectric_gate_boundary_edge_count"]),
        "production_reference_not_replaced": True,
    }
    checks["all_checks_pass"] = all(checks.values())
    if not checks["all_checks_pass"]: classification = "execution_or_identity_failure"
    current = [float(r["matched_ni_current_pair_growth_dex"]) for r in pair_rows]
    surface = [float(r["surface_response_pair_proxy_dex"]) for r in pair_rows]
    charge = [float(r["gate_charge_response_mismatch_pair_C_per_m"]) for r in pair_rows]
    report = {
        "schema": "vela.simplemos.sdevice.m71_gate_electrostatic_partition_report.v2",
        "status": "accepted" if checks["all_checks_pass"] else "failed", "classification": classification,
        "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)},
        "execution": {"new_sentaurus_solves": 0, "new_vela_solves": 0, "production_default_changed": False},
        "summary": {
            "median_surface_pair_closure_fraction": median_closure,
            "minimum_surface_pair_closure_fraction": min(closures), "maximum_surface_pair_closure_fraction": max(closures),
            "median_current_pair_growth_dex": median(current), "median_surface_response_pair_proxy_dex": median(surface),
            "surface_proxy_current_pearson": pearson(surface, current),
            "gate_charge_pair_current_pearson": pearson(charge, current),
            "maximum_gate_partition_identity_error_V": max(partition_errors),
            "maximum_sentaurus_gate_voltage_error_V": max(gate_voltage_errors),
            "maximum_sentaurus_reconstructed_to_native_gate_charge_relative_error": max(charge_errors),
            "median_sentaurus_reconstructed_to_native_gate_charge_relative_error": median(charge_errors),
            "gate_charge_proxy_qualified": charge_qualified,
            "minimum_oxide_gate_boundary_edge_count": min(int(r["oxide_gate_boundary_edge_count"]) for r in state_rows),
            "minimum_full_dielectric_gate_boundary_edge_count": min(int(r["full_dielectric_gate_boundary_edge_count"]) for r in state_rows),
        },
        "causal_scope": {
            "closed": "Quantifies the gate-stage surface/oxide voltage partition and charge response on frozen matched-ni states.",
            "not_closed": "Voltage partition is an identity and integrated charge is not a local finite-volume current closure; no production parameter is tuned."
        },
        "acceptance": checks,
    }
    return report, state_rows, pair_rows


def freeze_results(report: dict[str, Any], state_rows: list[dict[str, Any]], pair_rows: list[dict[str, Any]]) -> None:
    write_csv(STATES, state_rows); write_csv(PAIRS, pair_rows); write_json(REPORT, report)
    s = report["summary"]
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(f"""# SimpleMOS M71 栅控静电分解

M71 分类为 `{report['classification']}`。在 M65 matched-ni、no-BGN 和 M69 冻结状态上，表面势配对代理对剩余 NWell 电流配对增长的中位闭合率为 `{float(s['median_surface_pair_closure_fraction']):.2%}`，代理与电流配对增长的 Pearson 相关系数为 `{float(s['surface_proxy_current_pearson']):.4f}`。

氧化层压降与硅表面势是同一栅压分配恒等式的互补项；最大恒等式误差为 `{float(s['maximum_gate_partition_identity_error_V']):.3e} V`。全介质 P1 栅边界位移通量相对 Sentaurus 原生 ContactCharge 的最大相对误差为 `{float(s['maximum_sentaurus_reconstructed_to_native_gate_charge_relative_error']):.2%}`，观察器资格为 `{s['gate_charge_proxy_qualified']}`；未通过时栅电荷及积分半导体电荷仅作描述，不作因果闭合。该任务不重新打开 HFS、SG、接触提取、BGN、Nc/Nv、网格或收敛调参。

机器报告：`{portable(REPORT)}`。
""", encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {"schema": "vela.validation.artifact.v1", "title": "SimpleMOS M71 gate electrostatic partition",
                          "status": report["status"], "classification": report["classification"],
                          "report": portable(REPORT), "ledgers": [portable(STATES), portable(PAIRS)]})
    artifacts = [REPORT, STATES, PAIRS, DOC, ARTIFACT]
    write_json(EVIDENCE, {"schema": "vela.simplemos.sdevice.m71_gate_electrostatic_partition_evidence.v2",
                          "status": "frozen" if report["acceptance"]["all_checks_pass"] else "failed",
                          "classification": report["classification"], "contract_sha256": sha256(CONTRACT),
                          "v1_failure_evidence": portable(M71_V1_EVIDENCE),
                          "implementation_hashes": {portable(SCRIPT): sha256(SCRIPT)},
                          "artifacts": {portable(p): sha256(p) for p in artifacts},
                          "new_sentaurus_execution": False, "new_vela_execution": False,
                          "production_reference_replaced": False})


def verify() -> None:
    validate_contract(); evidence = read_json(EVIDENCE)
    for relative, expected in evidence["implementation_hashes"].items():
        if sha256(REPO / relative) != expected: raise ValueError(f"implementation changed: {relative}")
    for relative, expected in evidence["artifacts"].items():
        if sha256(REPO / relative) != expected: raise ValueError(f"artifact changed: {relative}")
    if evidence["status"] != "frozen": raise ValueError("M71 evidence is not accepted")


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--freeze", action="store_true")
    parser.add_argument("--analyze", action="store_true"); parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.freeze: freeze_contract()
    if args.analyze:
        report, state_rows, pair_rows = analyze(validate_contract()); freeze_results(report, state_rows, pair_rows)
        if not report["acceptance"]["all_checks_pass"]: raise SystemExit(2)
    if args.verify: verify()
    if not (args.freeze or args.analyze or args.verify): parser.error("choose --freeze, --analyze, or --verify")


if __name__ == "__main__": main()
