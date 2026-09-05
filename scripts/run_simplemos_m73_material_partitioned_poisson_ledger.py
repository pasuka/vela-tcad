"""Freeze and execute the SimpleMOS M73 material-partitioned Poisson ledger."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import factorized


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m73_material_partitioned_poisson_ledger_contract_v2.json"
FREEZE = ROOT / "simplemos_m73_material_partitioned_poisson_ledger_contract_freeze_v2.json"
M73_V1_EVIDENCE = ROOT / "simplemos_m73_material_partitioned_poisson_ledger_evidence_v1_superseded.json"
M72_EVIDENCE = ROOT / "simplemos_m72_native_gate_reaction_evidence.json"
M69_EVIDENCE = ROOT / "simplemos_m69_stage_residual_localization_evidence.json"
M69_STATES = ROOT / "stage_residual_localization/m69_stage_state_ledger.csv"
M69_EXPORTS = REPO / "build-release/m69_stage_v2/sentaurus_export_manifest.json"
M65_MANIFEST = REPO / "build-release/m65_ni/vela_manifest.json"
M65_PAIRS = ROOT / "nobgn_intrinsic_density_attribution/m65_nwell_pair_closure_ledger.csv"
M65_MATERIALS = REPO / "build-release/m65_ni/vela/matched_materials.json"
M8 = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics"
PORTABLE_ROOT = ROOT / "material_partitioned_poisson_ledger"
REPORT = PORTABLE_ROOT / "m73_material_partitioned_poisson_ledger_report.json"
STAGES = PORTABLE_ROOT / "m73_poisson_stage_ledger.csv"
PAIRS = PORTABLE_ROOT / "m73_poisson_pair_ledger.csv"
SUMMARIES = PORTABLE_ROOT / "m73_variant_summary.csv"
DOC = REPO / "docs/validation/simplemos_m73_material_partitioned_poisson_ledger_2026-09-03.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m73/artifact.json"
EVIDENCE = ROOT / "simplemos_m73_material_partitioned_poisson_ledger_evidence.json"
SCRIPT = Path(__file__).resolve()

EPS0 = 8.8541878128e-12
Q = 1.602176634e-19
VT_LN10 = 8.617333262145e-5 * 300.0 * math.log(10.0)
VARIANTS = {
    "legacy_all_cell": ("legacy", "all_cell"),
    "region_local_all_cell": ("region_local", "all_cell"),
    "legacy_signed_si": ("legacy", "signed_si"),
    "region_local_barycentric_si": ("region_local", "barycentric_si"),
    "region_local_signed_si": ("region_local", "signed_si"),
}
COMPONENTS = ("dielectric", "electron", "hole", "dopant")


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


def pearson(xs: list[float], ys: list[float]) -> float:
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    dx, dy = [x - mx for x in xs], [y - my for y in ys]
    denominator = math.sqrt(sum(x*x for x in dx) * sum(y*y for y in dy))
    return sum(x*y for x, y in zip(dx, dy, strict=True)) / denominator if denominator else 0.0


def closure(target: float, predicted: float) -> float:
    return max(0.0, min(1.0, 1.0 - abs(target - predicted) / max(abs(target), 1e-300)))


def scalar(path: Path) -> dict[int, float]:
    return {int(row["node_id"]): float(row["component0"]) for row in read_csv(path)}


def state(path: Path) -> dict[int, dict[str, float]]:
    return {int(row["node_id"]): {name: float(value) for name, value in row.items()
                                  if name != "node_id"}
            for row in read_csv(path)}


def workflows() -> dict[tuple[str, float], dict[str, Any]]:
    return {(row["device"], float(row["drain_voltage_V"])): row
            for row in read_json(M65_MANIFEST)["workflows"]}


def exports() -> dict[tuple[str, float, str], Path]:
    return {(row["device"], float(row["drain_voltage_V"]), row["stage"]):
            REPO / row["export_dir"]
            for row in read_json(M69_EXPORTS)["states"] if row["stage"] in ("drain", "gate")}


def vela_stage_path(workflow: dict[str, Any], stage_name: str) -> Path:
    final = REPO / workflow["state"]
    return final if stage_name == "gate" else final.parents[1] / "10_drain_ramp/state.csv"


def sent_required(export: Path) -> list[Path]:
    fields = export / "fields"
    paths = [export / "nodes.csv", fields / "eDensity_region0.csv",
             fields / "hDensity_region0.csv"]
    paths.extend(sorted(fields.glob("ElectrostaticPotential_region*.csv")))
    return paths


def source_paths() -> list[Path]:
    paths = [M72_EVIDENCE, M69_EVIDENCE, M73_V1_EVIDENCE, M69_STATES, M69_EXPORTS,
             M65_MANIFEST, M65_PAIRS, M65_MATERIALS, SCRIPT]
    for key, workflow in workflows().items():
        paths.extend([vela_stage_path(workflow, "drain"),
                      vela_stage_path(workflow, "gate"),
                      M8 / "vela" / key[0] / "mesh.json",
                      M8 / "neutral" / key[0] / "doping.csv"])
    for export in exports().values():
        paths.extend(sent_required(export))
    return sorted(set(path.resolve() for path in paths))


def freeze_contract() -> None:
    contract = read_json(CONTRACT)
    if contract.get("schema") != "vela.simplemos.sdevice.m73_material_partitioned_poisson_ledger_contract.v2":
        raise ValueError("unexpected M73 contract")
    m72, m69 = read_json(M72_EVIDENCE), read_json(M69_EVIDENCE)
    if (m72.get("status") != contract["upstream"]["required_m72_status"] or
            m72.get("classification") != contract["upstream"]["required_m72_classification"]):
        raise ValueError("M72 qualification changed")
    if (m69.get("status") != contract["upstream"]["required_m69_status"] or
            m69.get("classification") != contract["upstream"]["required_m69_classification"]):
        raise ValueError("M69 qualification changed")
    if read_json(M73_V1_EVIDENCE).get("status") != contract["upstream"]["required_m73_v1_status"]:
        raise ValueError("M73 v1 supersession evidence changed")
    paths = source_paths()
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(missing[0])
    write_json(FREEZE, {
        "schema": "vela.simplemos.sdevice.m73_material_partitioned_poisson_ledger_contract_freeze.v2",
        "status": "frozen_before_execution",
        "contract": portable(CONTRACT),
        "contract_sha256": sha256(CONTRACT),
        "upstream_hashes": {portable(path): sha256(path) for path in paths},
    })


def validate_contract() -> dict[str, Any]:
    contract, freeze = read_json(CONTRACT), read_json(FREEZE)
    if freeze.get("status") != "frozen_before_execution":
        raise ValueError("M73 contract is not frozen")
    if freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M73 frozen contract changed")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M73 frozen input changed: {relative}")
    return contract


def triangle_area(points: list[tuple[float, float]]) -> float:
    (ax, ay), (bx, by), (cx, cy) = points
    return 0.5 * abs((bx-ax)*(cy-ay) - (cx-ax)*(by-ay))


def triangle_signed_measures(points: list[tuple[float, float]]) -> list[float]:
    def d2(a: int, b: int) -> float:
        return sum((points[b][k] - points[a][k])**2 for k in (0, 1))
    def cot(vertex: int, a: int, b: int) -> float:
        ux, uy = points[a][0]-points[vertex][0], points[a][1]-points[vertex][1]
        vx, vy = points[b][0]-points[vertex][0], points[b][1]-points[vertex][1]
        return (ux*vx + uy*vy) / abs(ux*vy - uy*vx)
    return [0.125 * (d2(i, (i+2)%3)*cot((i+1)%3, i, (i+2)%3)
                     + d2(i, (i+1)%3)*cot((i+2)%3, i, (i+1)%3))
            for i in range(3)]


def local_couple(points: list[tuple[float, float]], local_edge: tuple[int, int]) -> float:
    a, b = local_edge
    opposite = next(i for i in range(3) if i not in local_edge)
    ux, uy = points[a][0]-points[opposite][0], points[a][1]-points[opposite][1]
    vx, vy = points[b][0]-points[opposite][0], points[b][1]-points[opposite][1]
    cross = abs(ux*vy - uy*vx)
    length = math.hypot(points[b][0]-points[a][0], points[b][1]-points[a][1])
    if cross <= 1e-300 or length <= 1e-30:
        return 0.0
    cotangent = (ux*vx + uy*vy) / cross
    return triangle_area(points)/(3.0*length) if cotangent < 0.0 else 0.5*cotangent*length


class Geometry:
    def __init__(self, device: str):
        mesh = read_json(M8 / "vela" / device / "mesh.json")
        self.coords = {int(row["id"]): (float(row["x"]), float(row["y"]))
                       for row in mesh["nodes"]}
        self.count = len(self.coords)
        regions = {int(row["id"]): row["material"] for row in mesh["regions"]}
        eps = {row["name"]: float(row["eps_r"])*EPS0
               for row in read_json(M65_MATERIALS)["materials"]}
        all_volume = np.zeros(self.count)
        bary_si = np.zeros(self.count)
        signed_si = np.zeros(self.count)
        edge_parts: dict[tuple[int, int], list[tuple[str, float]]] = {}
        for tri in mesh["triangles"]:
            ids = list(map(int, tri["node_ids"]))
            points = [self.coords[node] for node in ids]
            material = regions[int(tri["region_id"])]
            area = triangle_area(points)
            signed = triangle_signed_measures(points)
            for local, node in enumerate(ids):
                all_volume[node] += area/3.0
                if material == "Si":
                    bary_si[node] += area/3.0
                    signed_si[node] += signed[local]
            for a, b in ((0, 1), (1, 2), (2, 0)):
                edge = tuple(sorted((ids[a], ids[b])))
                edge_parts.setdefault(edge, []).append((material, local_couple(points, (a, b))))
        self.volumes = {"all_cell": all_volume*1e-12,
                        "barycentric_si": bary_si*1e-12,
                        "signed_si": signed_si*1e-12}
        matrix_rows: dict[str, list[int]] = {"legacy": [], "region_local": []}
        matrix_cols: dict[str, list[int]] = {"legacy": [], "region_local": []}
        matrix_data: dict[str, list[float]] = {"legacy": [], "region_local": []}
        for (i, j), parts in edge_parts.items():
            length = math.hypot(self.coords[j][0]-self.coords[i][0],
                                self.coords[j][1]-self.coords[i][1])
            total_couple = sum(value for _, value in parts)
            legacy = (sum(eps[name] for name, _ in parts)/len(parts))*total_couple/length
            resolved = sum(eps[name]*value for name, value in parts)/length
            for kind, coefficient in (("legacy", legacy), ("region_local", resolved)):
                matrix_rows[kind].extend((i, i, j, j))
                matrix_cols[kind].extend((i, j, i, j))
                matrix_data[kind].extend((coefficient, -coefficient, -coefficient, coefficient))
        self.matrices = {kind: coo_matrix((matrix_data[kind],
                                           (matrix_rows[kind], matrix_cols[kind])),
                                          shape=(self.count, self.count)).tocsr()
                         for kind in matrix_rows}
        self.contact_nodes = sorted({int(node) for contact in mesh["contacts"]
                                     for node in contact["node_ids"]})
        contact = set(self.contact_nodes)
        self.free = np.array([node for node in range(self.count) if node not in contact], dtype=int)
        self.solvers = {kind: factorized(matrix[self.free][:, self.free].tocsc())
                        for kind, matrix in self.matrices.items()}

    def solve(self, kind: str, residual: np.ndarray) -> np.ndarray:
        answer = np.zeros(self.count)
        answer[self.free] = self.solvers[kind](-residual[self.free])
        return answer


def sentaurus_state(export: Path, geometry: Geometry) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    fields = export / "fields"
    potentials: dict[int, list[float]] = {}
    priority: dict[int, float] = {}
    for path in sorted(fields.glob("ElectrostaticPotential_region*.csv")):
        region = int(path.stem.rsplit("region", 1)[1])
        for node, value in scalar(path).items():
            potentials.setdefault(node, []).append(value)
            if node not in priority or region == 0:
                priority[node] = value
    if set(priority) != set(range(geometry.count)):
        raise ValueError(f"Sentaurus potential support is incomplete: {export}")
    spread = max((max(values)-min(values) for values in potentials.values()), default=0.0)
    psi = np.array([priority[node] for node in range(geometry.count)])
    n, p = np.zeros(geometry.count), np.zeros(geometry.count)
    for node, value in scalar(fields / "eDensity_region0.csv").items():
        n[node] = value*1e6
    for node, value in scalar(fields / "hDensity_region0.csv").items():
        p[node] = value*1e6
    return psi, n, p, spread


def coordinate_error(export: Path, geometry: Geometry) -> float:
    rows = read_csv(export / "nodes.csv")
    return max(max(abs(float(row[name])-geometry.coords[int(row["id"])][axis])
                   for axis, name in enumerate(("x_um", "y_um"))) for row in rows)


def doping(device: str, count: int) -> np.ndarray:
    result = np.zeros(count)
    for row in read_csv(M8 / "neutral" / device / "doping.csv"):
        node = int(row["node_id"])
        result[node] = (float(row["donors_cm3"])-float(row["acceptors_cm3"]))*1e6
    return result


def control_nodes() -> dict[tuple[str, float], int]:
    rows = read_csv(M69_STATES)
    return {(row["device"], float(row["drain_voltage_V"])): int(row["sentaurus_control_node"])
            for row in rows if row["stage"] == "gate"}


def analyze(contract: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    workflow_map, export_map, controls = workflows(), exports(), control_nodes()
    target_rows = read_csv(M65_PAIRS)
    target_map = {(row["low_device"], row["high_device"], float(row["drain_voltage_V"])): row
                  for row in target_rows}
    geometry_cache: dict[str, Geometry] = {}
    stage_rows: list[dict[str, Any]] = []
    max_coordinate = max_shared_spread = max_vela_replay = max_component_sum = 0.0
    for key, workflow in sorted(workflow_map.items()):
        device, drain = key
        geometry = geometry_cache.setdefault(device, Geometry(device))
        net_doping = doping(device, geometry.count)
        control = controls[key]
        for stage_name in ("drain", "gate"):
            export = export_map[(device, drain, stage_name)]
            max_coordinate = max(max_coordinate, coordinate_error(export, geometry))
            psi, n, p, spread = sentaurus_state(export, geometry)
            max_shared_spread = max(max_shared_spread, spread)
            vela = state(vela_stage_path(workflow, stage_name))
            vpsi = np.array([vela[node]["psi"] for node in range(geometry.count)])
            vn = np.array([vela[node]["electrons_m3"] for node in range(geometry.count)])
            vp = np.array([vela[node]["holes_m3"] for node in range(geometry.count)])
            legacy_charge = Q*(vn-vp-net_doping)*geometry.volumes["all_cell"]
            vela_residual = geometry.matrices["legacy"]@vpsi + legacy_charge
            vela_correction = geometry.solve("legacy", vela_residual)
            max_vela_replay = max(max_vela_replay, float(np.max(np.abs(vela_correction))))
            for variant, (matrix_kind, volume_kind) in VARIANTS.items():
                volume = geometry.volumes[volume_kind]
                components = {
                    "dielectric": geometry.matrices[matrix_kind]@psi,
                    "electron": Q*n*volume,
                    "hole": -Q*p*volume,
                    "dopant": -Q*net_doping*volume,
                }
                responses = {name: geometry.solve(matrix_kind, value)
                             for name, value in components.items()}
                total_residual = sum(components.values(), np.zeros(geometry.count))
                total_response = geometry.solve(matrix_kind, total_residual)
                summed_response = sum(responses.values(), np.zeros(geometry.count))
                identity = float(np.max(np.abs(total_response-summed_response)))
                max_component_sum = max(max_component_sum, identity)
                stage_rows.append({
                    "device": device, "drain_voltage_V": drain,
                    "gate_voltage_V": 0.0 if stage_name == "drain" else float(workflow["gate_voltage_V"]),
                    "stage": stage_name, "control_node": control, "variant": variant,
                    "dielectric_response_V": responses["dielectric"][control],
                    "electron_response_V": responses["electron"][control],
                    "hole_response_V": responses["hole"][control],
                    "dopant_response_V": responses["dopant"][control],
                    "total_response_V": total_response[control],
                    "component_sum_identity_error_V": identity,
                    "free_residual_l2_C_per_m": float(np.linalg.norm(total_residual[geometry.free])),
                    "vela_production_replay_max_correction_V": float(np.max(np.abs(vela_correction))),
                    "coordinate_error_um": coordinate_error(export, geometry),
                    "sentaurus_shared_node_potential_spread_V": spread,
                })
    stage_by = {(row["device"], float(row["drain_voltage_V"]), row["stage"], row["variant"]): row
                for row in stage_rows}
    case_response: dict[tuple[str, float, str], dict[str, float]] = {}
    for device, drain in sorted(workflow_map):
        for variant in VARIANTS:
            dr = stage_by[(device, drain, "drain", variant)]
            ga = stage_by[(device, drain, "gate", variant)]
            case_response[(device, drain, variant)] = {
                name: float(ga[f"{name}_response_V"])-float(dr[f"{name}_response_V"])
                for name in (*COMPONENTS, "total")}
    pair_rows: list[dict[str, Any]] = []
    for pair in contract["matrix"]["pairs"]:
        low, high, drain = pair["low_device"], pair["high_device"], float(pair["drain_voltage_V"])
        targets = target_map[(low, high, drain)]
        current_target = float(targets["matched_ni_pair_growth_dex"])
        electro_target = float(targets["matched_ni_electrostatic_barrier_pair_proxy_dex"])
        absolute: dict[str, dict[str, float]] = {}
        for variant in VARIANTS:
            lo, hi = case_response[(low, drain, variant)], case_response[(high, drain, variant)]
            absolute[variant] = {name: (hi[name]-lo[name])/VT_LN10
                                 for name in (*COMPONENTS, "total")}
        legacy = absolute["legacy_all_cell"]
        for variant in VARIANTS:
            proxies = absolute[variant]
            attributed = (legacy if variant == "legacy_all_cell" else
                          {name: legacy[name]-proxies[name] for name in (*COMPONENTS, "total")})
            pair_rows.append({
                "low_device": low, "high_device": high, "drain_voltage_V": drain,
                "diagnostic_gate_voltage_V": pair["diagnostic_gate_voltage_V"], "variant": variant,
                "current_pair_growth_target_dex": current_target,
                "electrostatic_pair_proxy_target_dex": electro_target,
                **{f"absolute_{name}_pair_proxy_dex": proxies[name] for name in (*COMPONENTS, "total")},
                **{f"attributed_{name}_pair_proxy_dex": attributed[name] for name in (*COMPONENTS, "total")},
                "same_sign_as_current": attributed["total"]*current_target > 0.0,
                "current_closure_fraction": closure(current_target, attributed["total"]),
                "same_sign_as_electrostatic": attributed["total"]*electro_target > 0.0,
                "electrostatic_closure_fraction": closure(electro_target, attributed["total"]),
            })
    summaries = []
    for variant in VARIANTS:
        rows = [row for row in pair_rows if row["variant"] == variant]
        predictions = [float(row["attributed_total_pair_proxy_dex"]) for row in rows]
        current = [float(row["current_pair_growth_target_dex"]) for row in rows]
        electro = [float(row["electrostatic_pair_proxy_target_dex"]) for row in rows]
        summaries.append({
            "variant": variant,
            "same_sign_current_pair_count": sum(bool(row["same_sign_as_current"]) for row in rows),
            "median_current_closure_fraction": median([float(row["current_closure_fraction"]) for row in rows]),
            "minimum_current_closure_fraction": min(float(row["current_closure_fraction"]) for row in rows),
            "prediction_current_pearson": pearson(predictions, current),
            "same_sign_electrostatic_pair_count": sum(bool(row["same_sign_as_electrostatic"]) for row in rows),
            "median_electrostatic_closure_fraction": median([float(row["electrostatic_closure_fraction"]) for row in rows]),
            "minimum_electrostatic_closure_fraction": min(float(row["electrostatic_closure_fraction"]) for row in rows),
            "prediction_electrostatic_pearson": pearson(predictions, electro),
            "median_attributed_total_pair_proxy_dex": median(predictions),
        })
    candidates = [row for row in summaries if row["variant"] != "legacy_all_cell"]
    best = max(candidates, key=lambda row: (float(row["median_current_closure_fraction"]),
                                            int(row["same_sign_current_pair_count"])))
    best_electro = max(candidates, key=lambda row: (
        float(row["median_electrostatic_closure_fraction"]),
        int(row["same_sign_electrostatic_pair_count"])))
    dominant = (int(best["same_sign_current_pair_count"]) >=
                int(contract["analysis"]["dominant_minimum_same_sign_pair_count"]) and
                float(best["median_current_closure_fraction"]) >=
                float(contract["analysis"]["dominant_minimum_median_current_closure_fraction"]))
    material = (float(best["median_current_closure_fraction"]) >=
                float(contract["analysis"]["material_minimum_median_current_closure_fraction"]))
    classification = ("material_partition_dominant" if dominant else
                      "material_partition_material_but_not_dominant" if material else
                      "material_partition_not_material")
    limits = contract["acceptance"]
    checks = {
        "contract_frozen": read_json(FREEZE)["contract_sha256"] == sha256(CONTRACT),
        "pair_count": len(pair_rows) == int(limits["required_pair_count"])*len(VARIANTS),
        "case_count": len(workflow_map) == int(limits["required_case_count"]),
        "stage_count": len(stage_rows) == int(limits["required_stage_count"])*len(VARIANTS),
        "coordinate_identity": max_coordinate <= float(limits["maximum_coordinate_error_um"]),
        "shared_node_potential_identity": max_shared_spread <= float(limits["maximum_sentaurus_shared_node_potential_spread_V"]),
        "vela_production_poisson_replay": max_vela_replay <= float(limits["maximum_vela_production_poisson_replay_correction_V"]),
        "component_sum_identity": max_component_sum <= float(limits["maximum_component_sum_identity_error_V"]),
        "production_reference_not_replaced": True,
    }
    checks["all_checks_pass"] = all(checks.values())
    if not checks["all_checks_pass"]:
        classification = "execution_or_identity_failure"
    legacy = next(row for row in summaries if row["variant"] == "legacy_all_cell")
    report = {
        "schema": "vela.simplemos.sdevice.m73_material_partitioned_poisson_ledger_report.v2",
        "status": "accepted" if checks["all_checks_pass"] else "failed",
        "classification": classification,
        "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)},
        "execution": {"new_sentaurus_solves": 0, "new_vela_self_consistent_solves": 0,
                      "offline_sparse_poisson_stage_responses": int(limits["required_stage_count"])*len(VARIANTS),
                      "production_default_changed": False},
        "summary": {
            "best_material_partition_variant": best["variant"],
            "best_same_sign_current_pair_count": best["same_sign_current_pair_count"],
            "best_median_current_closure_fraction": best["median_current_closure_fraction"],
            "best_current_variant_median_electrostatic_closure_fraction": best["median_electrostatic_closure_fraction"],
            "best_electrostatic_partition_variant": best_electro["variant"],
            "best_same_sign_electrostatic_pair_count": best_electro["same_sign_electrostatic_pair_count"],
            "best_median_electrostatic_closure_fraction": best_electro["median_electrostatic_closure_fraction"],
            "legacy_median_current_closure_fraction": legacy["median_current_closure_fraction"],
            "maximum_vela_production_poisson_replay_correction_V": max_vela_replay,
            "maximum_component_sum_identity_error_V": max_component_sum,
            "maximum_coordinate_error_um": max_coordinate,
            "maximum_sentaurus_shared_node_potential_spread_V": max_shared_spread,
        },
        "variant_summaries": summaries,
        "causal_scope": {
            "closed": "Tests material-local dielectric and semiconductor-only charge-volume hypotheses against all eight frozen high-NWell pairs.",
            "not_closed": "No diagnostic operator is promoted to a production default without a later self-consistent full-matrix qualification."
        },
        "acceptance": checks,
    }
    return report, stage_rows, pair_rows, summaries


def freeze_results(report: dict[str, Any], stage_rows: list[dict[str, Any]],
                   pair_rows: list[dict[str, Any]], summaries: list[dict[str, Any]]) -> None:
    write_csv(STAGES, stage_rows)
    write_csv(PAIRS, pair_rows)
    write_csv(SUMMARIES, summaries)
    write_json(REPORT, report)
    summary = report["summary"]
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(f"""# SimpleMOS M73 材料分区 Poisson 账本

M73 分类为 `{report['classification']}`。任务复用 M65/M69 的 32 个 drain/gate 冻结状态，没有新增 Sentaurus 或 Vela 自洽求解。

五种离线 Poisson 观察器中，对剩余电流配对增长解释最好的材料分区候选为 `{summary['best_material_partition_variant']}`：同号 `{summary['best_same_sign_current_pair_count']}/8`，中位闭合率 `{float(summary['best_median_current_closure_fraction']):.2%}`。对 M65 静电势垒配对代理解释最好的是 `{summary['best_electrostatic_partition_variant']}`：同号 `{summary['best_same_sign_electrostatic_pair_count']}/8`，中位闭合率 `{float(summary['best_median_electrostatic_closure_fraction']):.2%}`。生产 `legacy_all_cell` 冻结修正本身对电流目标的中位闭合率为 `{float(summary['legacy_median_current_closure_fraction']):.2%}`。

Vela 生产状态回放的最大 Poisson 线性修正为 `{float(summary['maximum_vela_production_poisson_replay_correction_V']):.3e}` V，分项线性叠加最大误差为 `{float(summary['maximum_component_sum_identity_error_V']):.3e}` V。上述结果只用于固定状态归因；未修改 HFS、SG、接触提取、准费米打包、BGN、Nc/Nv、网格、收敛参数或生产默认值。
""", encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {
        "schema": "vela.validation.artifact.v1", "title": "SimpleMOS M73 material-partitioned Poisson ledger",
        "status": report["status"], "classification": report["classification"],
        "report": portable(REPORT), "document": portable(DOC),
        "ledgers": [portable(STAGES), portable(PAIRS), portable(SUMMARIES)],
    })
    artifacts = [REPORT, STAGES, PAIRS, SUMMARIES, DOC, ARTIFACT]
    write_json(EVIDENCE, {
        "schema": "vela.simplemos.sdevice.m73_material_partitioned_poisson_ledger_evidence.v2",
        "status": "frozen" if report["acceptance"]["all_checks_pass"] else "failed",
        "classification": report["classification"], "contract_sha256": sha256(CONTRACT),
        "implementation_hashes": {portable(SCRIPT): sha256(SCRIPT)},
        "artifacts": {portable(path): sha256(path) for path in artifacts},
        "new_sentaurus_execution": False, "new_vela_self_consistent_execution": False,
        "production_reference_replaced": False, "acceptance": report["acceptance"],
    })


def verify() -> dict[str, Any]:
    validate_contract()
    evidence, report = read_json(EVIDENCE), read_json(REPORT)
    if evidence["status"] != "frozen" or report["status"] != "accepted":
        raise ValueError("M73 is not frozen")
    if evidence["implementation_hashes"][portable(SCRIPT)] != sha256(SCRIPT):
        raise ValueError("M73 implementation changed")
    for relative, expected in evidence["artifacts"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M73 artifact changed: {relative}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze-contract", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.freeze_contract:
        freeze_contract()
        print(json.dumps({"status": "frozen_before_execution", "contract_sha256": sha256(CONTRACT)}))
        return
    contract = validate_contract()
    if args.verify:
        print(json.dumps(verify()["summary"], indent=2))
        return
    if not args.analyze:
        parser.error("choose --freeze-contract, --analyze, or --verify")
    report, stage_rows, pair_rows, summaries = analyze(contract)
    freeze_results(report, stage_rows, pair_rows, summaries)
    print(json.dumps({"status": report["status"], "classification": report["classification"],
                      "summary": report["summary"], "acceptance": report["acceptance"]}, indent=2))


if __name__ == "__main__":
    main()
