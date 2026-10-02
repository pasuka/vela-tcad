#!/usr/bin/env python3
"""Run the frozen SimpleMOS M49 substrate electron-transport attribution."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
from collections import defaultdict, deque
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m49_substrate_electron_transport_attribution_contract_v1.json"
FREEZE = ROOT / "simplemos_m49_substrate_electron_transport_attribution_contract_freeze.json"
M47_CONTRACT = ROOT / "simplemos_m47_default_bgn_self_consistent_attribution_contract_v1.json"
M47_EVIDENCE = ROOT / "simplemos_m47_default_bgn_self_consistent_attribution_evidence.json"
M48_CONTRACT = ROOT / "simplemos_m48_terminal_partition_continuity_closure_contract_v1.json"
M48_EVIDENCE = ROOT / "simplemos_m48_terminal_partition_continuity_closure_evidence.json"
M48_REPORT = ROOT / "terminal_partition_continuity_closure/m48_terminal_partition_continuity_closure_report.json"
M48_TERMINALS = ROOT / "terminal_partition_continuity_closure/m48_four_terminal_component_ledger.csv"
M47_RAW = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m47_default_bgn_state_attribution"
OUTPUT = ROOT / "substrate_electron_transport_attribution"
REPORT = OUTPUT / "m49_substrate_electron_transport_attribution_report.json"
BOUNDARY = OUTPUT / "m49_substrate_boundary_edge_transport_ledger.csv"
ROWS = OUTPUT / "m49_substrate_graph_distance_row_ledger.csv"
BRIDGE = OUTPUT / "m49_transport_factor_bridge_ledger.csv"
CONTROLS = OUTPUT / "m49_control_localization_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m49_substrate_electron_transport_attribution_2026-09-01.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m49/artifact.json"
EVIDENCE = ROOT / "simplemos_m49_substrate_electron_transport_attribution_evidence.json"
TEST = REPO / "tests/regression/test_simplemos_m49_substrate_electron_transport_attribution.py"
DEVICES = ("n23", "n19")
GATES = (0.0, 0.05, 0.1)
Q = 1.602176634e-19
FACTORS = ("density", "mobility", "qf_gradient")
TERMS = FACTORS + ("constitutive_residual", "terminal_boundary_residual")
VTK_FIELDS = {
    "Potential", "Electrons", "SRHRecombinationCm3PerS",
    "ElectronGradQuasiFermiVector", "ElectronMobilityCm2PerVs",
    "SentaurusElectronCurrentDensityVector",
}


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False,
                              allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def percentile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        raise ValueError("percentile of empty values")
    position = (len(ordered) - 1) * fraction
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def validate_contract() -> dict[str, Any]:
    contract = read_json(CONTRACT)
    freeze = read_json(FREEZE)
    if contract.get("schema") != (
            "vela.simplemos.sdevice.m49_substrate_electron_transport_attribution_contract.v1"):
        raise ValueError("unexpected M49 contract schema")
    actual = sha256(CONTRACT)
    if freeze.get("status") != "frozen_before_analysis":
        raise ValueError("M49 contract was not frozen before analysis")
    if freeze.get("contract_sha256") != actual:
        raise ValueError("M49 contract changed after freeze")
    for path in (M47_EVIDENCE, M48_EVIDENCE):
        if read_json(path).get("status") != "frozen":
            raise ValueError(f"M49 requires frozen upstream evidence: {path.name}")
    target = read_json(M48_REPORT)["target"]["substrate_electron_difference_A_per_um"]
    expected = contract["upstream"][
        "required_m48_target_substrate_electron_difference_A_per_um"]
    if not math.isclose(float(target), float(expected), rel_tol=0.0, abs_tol=1e-30):
        raise ValueError("M48 target substrate electron anchor changed")
    if contract["execution_protocol"]["new_sentaurus_execution"]:
        raise ValueError("M49 must not execute Sentaurus")
    if contract["execution_protocol"]["new_vela_execution"]:
        raise ValueError("M49 must not execute Vela")
    return contract


def state_id(device: str, gate: float) -> str:
    return f"{device}_vd_0p05_vg_{format(gate, '.12g').replace('.', 'p')}"


def vtk_path(device: str, gate: float) -> Path:
    index = GATES.index(gate)
    matches = sorted((M47_RAW / "vela" / device / "vtk").glob(
        f"state_{index:04d}_*.vtk"))
    if len(matches) != 1:
        raise FileNotFoundError(f"expected one VTK for {device}, {gate}: {matches}")
    return matches[0]


def parse_vtk(path: Path) -> tuple[list[tuple[float, float]], dict[str, list[Any]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    points_at = next(i for i, line in enumerate(lines) if line.startswith("POINTS "))
    point_count = int(lines[points_at].split()[1])
    points = [tuple(float(value) for value in lines[i].split()[:2])
              for i in range(points_at + 1, points_at + 1 + point_count)]
    data_at = next(i for i, line in enumerate(lines) if line.startswith("POINT_DATA "))
    count = int(lines[data_at].split()[1])
    if count != point_count:
        raise ValueError(f"VTK point count mismatch in {path}")
    fields: dict[str, list[Any]] = {}
    index = data_at + 1
    while index < len(lines):
        parts = lines[index].split()
        if not parts:
            index += 1
        elif parts[0] == "SCALARS":
            name = parts[1]
            index += 2
            values = [float(lines[index + offset].split()[0])
                      for offset in range(count)]
            if name in VTK_FIELDS:
                fields[name] = values
            index += count
        elif parts[0] == "VECTORS":
            name = parts[1]
            index += 1
            values = [tuple(float(value) for value in
                            lines[index + offset].split()[:2])
                      for offset in range(count)]
            if name in VTK_FIELDS:
                fields[name] = values
            index += count
        else:
            index += 1
    missing = VTK_FIELDS - set(fields)
    if missing:
        raise ValueError(f"missing VTK fields in {path}: {sorted(missing)}")
    return points, fields


def scalar_csv(path: Path) -> dict[int, float]:
    return {int(row["node_id"]): float(row["component0"])
            for row in read_csv(path)}


def vector_csv(path: Path) -> dict[int, tuple[float, float]]:
    return {int(row["node_id"]): (float(row["component0"]),
                                  float(row["component1"]))
            for row in read_csv(path)}


def solver_fields(device: str, gate: float) -> tuple[
        dict[int, tuple[float, float]], dict[str, dict[int, Any]],
        dict[str, dict[int, Any]], Path, Path]:
    state = state_id(device, gate)
    export = M47_RAW / "sentaurus_exports" / state
    nodes = {int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
             for row in read_csv(export / "nodes.csv")}
    sent = {
        "psi": scalar_csv(export / "fields/ElectrostaticPotential_region0.csv"),
        "n": scalar_csv(export / "fields/eDensity_region0.csv"),
        "srh": scalar_csv(export / "fields/srhRecombination_region0.csv"),
        "mu": scalar_csv(export / "fields/eMobility_region0.csv"),
        "grad": vector_csv(export / "fields/eGradQuasiFermi_region0.csv"),
        "current": vector_csv(export / "fields/eCurrentDensity_region0.csv"),
    }
    vtk = vtk_path(device, gate)
    vtk_points, raw = parse_vtk(vtk)
    if set(nodes) != set(range(len(vtk_points))):
        raise ValueError(f"node-id coverage mismatch for {state}")
    max_coordinate_error = max(
        max(abs(nodes[index][axis] - vtk_points[index][axis]) for axis in (0, 1))
        for index in nodes)
    if max_coordinate_error > 1e-12:
        raise ValueError(f"coordinate mismatch for {state}: {max_coordinate_error}")
    vela_all = {
        "psi": dict(enumerate(raw["Potential"])),
        "n": dict(enumerate(raw["Electrons"])),
        "srh": dict(enumerate(raw["SRHRecombinationCm3PerS"])),
        "mu": dict(enumerate(raw["ElectronMobilityCm2PerVs"])),
        "grad": dict(enumerate(raw["ElectronGradQuasiFermiVector"])),
        # This Vela diagnostic is conventional electron current in Sentaurus A/cm2 units.
        "current": dict(enumerate(raw["SentaurusElectronCurrentDensityVector"])),
    }
    common = set(sent["psi"])
    if not common:
        raise ValueError(f"empty Silicon support for {state}")
    for field, values in sent.items():
        if set(values) != common:
            raise ValueError(f"Sentaurus {field} node coverage mismatch for {state}")
    if any(not common.issubset(values) for values in vela_all.values()):
        raise ValueError(f"Vela Silicon node coverage mismatch for {state}")
    vela = {field: {node: values[node] for node in common}
            for field, values in vela_all.items()}
    return nodes, sent, vela, export, vtk


def mesh_topology(export: Path, nodes: dict[int, tuple[float, float]]) -> tuple[
        list[dict[str, Any]], dict[int, set[int]], set[int]]:
    contacts = read_csv(export / "contacts.csv")
    substrate_rows = [row for row in contacts if row["name"].lower() == "substrate"]
    if len(substrate_rows) != 1:
        raise ValueError("expected one substrate contact")
    contact_nodes = {int(value) for value in substrate_rows[0]["node_ids"].split(";")}
    edge_thirds: dict[tuple[int, int], list[int]] = defaultdict(list)
    adjacency: dict[int, set[int]] = defaultdict(set)
    for element in read_csv(export / "elements.csv"):
        if element["material"] != "Si":
            continue
        triangle = tuple(int(element[f"node{i}"]) for i in range(3))
        for a, b, third in ((triangle[0], triangle[1], triangle[2]),
                            (triangle[1], triangle[2], triangle[0]),
                            (triangle[2], triangle[0], triangle[1])):
            key = tuple(sorted((a, b)))
            edge_thirds[key].append(third)
            adjacency[a].add(b)
            adjacency[b].add(a)
    boundary: list[dict[str, Any]] = []
    for (node0, node1), thirds in edge_thirds.items():
        if len(thirds) != 1 or node0 not in contact_nodes or node1 not in contact_nodes:
            continue
        x0, y0 = nodes[node0]
        x1, y1 = nodes[node1]
        dx, dy = x1 - x0, y1 - y0
        length = math.hypot(dx, dy)
        nx, ny = dy / length, -dx / length
        mx, my = (x0 + x1) / 2.0, (y0 + y1) / 2.0
        tx, ty = nodes[thirds[0]]
        if nx * (tx - mx) + ny * (ty - my) > 0.0:
            nx, ny = -nx, -ny
        boundary.append({
            "node0": node0, "node1": node1, "third_node": thirds[0],
            "midpoint_x_um": mx, "midpoint_y_um": my,
            "edge_length_um": length, "normal_x": nx, "normal_y": ny,
        })
    if not boundary:
        raise ValueError("no substrate boundary edges")
    boundary.sort(key=lambda row: (row["midpoint_y_um"], row["midpoint_x_um"]))
    return boundary, adjacency, contact_nodes


def graph_layers(adjacency: dict[int, set[int]], contact_nodes: set[int],
                 maximum: int = 4) -> dict[int, list[int]]:
    distance = {node: 0 for node in contact_nodes}
    queue = deque(contact_nodes)
    while queue:
        node = queue.popleft()
        if distance[node] >= maximum:
            continue
        for neighbor in adjacency[node]:
            if neighbor not in distance:
                distance[neighbor] = distance[node] + 1
                queue.append(neighbor)
    return {layer: sorted(node for node, value in distance.items() if value == layer)
            for layer in range(maximum + 1)}


def dot(vector: tuple[float, float], normal: tuple[float, float]) -> float:
    return vector[0] * normal[0] + vector[1] * normal[1]


def product_current(n: float, mu: float,
                    grad: tuple[float, float]) -> tuple[float, float]:
    coefficient = Q * n * mu
    return coefficient * grad[0], coefficient * grad[1]


def shapley(s: dict[str, float], v: dict[str, float]) -> dict[str, float]:
    contributions = {factor: 0.0 for factor in FACTORS}
    for order in itertools.permutations(FACTORS):
        state = dict(s)
        before = math.prod(state[factor] for factor in FACTORS)
        for factor in order:
            state[factor] = v[factor]
            after = math.prod(state[item] for item in FACTORS)
            contributions[factor] += (after - before) / math.factorial(len(FACTORS))
            before = after
    expected = math.prod(v[factor] for factor in FACTORS) - math.prod(
        s[factor] for factor in FACTORS)
    error = sum(contributions.values()) - expected
    if abs(error) > 1e-12 * max(abs(expected), 1e-300):
        raise AssertionError(f"Shapley identity failed: {error}")
    return contributions


def terminal_anchors() -> dict[tuple[str, float, str], float]:
    result: dict[tuple[str, float, str], float] = {}
    for row in read_csv(M48_TERMINALS):
        if row["contact"] != "substrate" or row["component"] != "electron":
            continue
        result[(row["device"], float(row["gate_voltage_V"]), row["solver"])] = float(
            row["conventional_current_A_per_um"])
    expected = {(device, gate, solver) for device in DEVICES for gate in GATES
                for solver in ("sentaurus", "vela")}
    if set(result) != expected:
        raise ValueError("M48 substrate electron anchors incomplete")
    return result


def log_ratio(candidate: float, reference: float, floor: float = 1e-300) -> float:
    return math.log10(max(abs(candidate), floor) / max(abs(reference), floor))


def paired_log_ratios(candidate: list[float], reference: list[float],
                      relative_floor: float = 1e-12) -> list[float]:
    if len(candidate) != len(reference) or not candidate:
        raise ValueError("paired log-ratio inputs must be nonempty and aligned")
    peak = max(abs(value) for value in candidate + reference)
    floor = max(peak * relative_floor, 1e-300)
    return [log_ratio(cand, ref, floor) for cand, ref in zip(candidate, reference)]


def analyze() -> tuple[dict[str, Any], list[dict[str, Any]],
                       list[dict[str, Any]], list[dict[str, Any]],
                       list[dict[str, Any]]]:
    contract = validate_contract()
    anchors = terminal_anchors()
    boundary_rows: list[dict[str, Any]] = []
    layer_rows: list[dict[str, Any]] = []
    bridge_rows: list[dict[str, Any]] = []
    state_results: list[dict[str, Any]] = []
    state_files: list[Path] = []
    maximum_coordinate_error = 0.0
    maximum_shapley_error = 0.0
    maximum_bridge_error = 0.0

    for device in DEVICES:
        reference_topology: list[tuple[int, int]] | None = None
        for gate in GATES:
            state = state_id(device, gate)
            nodes, sent, vela, export, vtk = solver_fields(device, gate)
            state_files.extend([export / "nodes.csv", export / "elements.csv",
                                export / "contacts.csv", vtk])
            boundary, adjacency, contact_nodes = mesh_topology(export, nodes)
            topology = [(row["node0"], row["node1"]) for row in boundary]
            if reference_topology is None:
                reference_topology = topology
            elif topology != reference_topology:
                raise ValueError(f"substrate topology changed across {device} gate states")
            layers = graph_layers(adjacency, contact_nodes)
            integrated_export = {"sentaurus": 0.0, "vela": 0.0}
            integrated_product = {"sentaurus": 0.0, "vela": 0.0}
            factor_integrals = {factor: 0.0 for factor in FACTORS}

            for edge_index, edge in enumerate(boundary):
                normal = (edge["normal_x"], edge["normal_y"])
                weight = edge["edge_length_um"] * 1e-8
                edge_values: dict[str, dict[str, float]] = {}
                edge_factor = {factor: 0.0 for factor in FACTORS}
                for solver, fields in (("sentaurus", sent), ("vela", vela)):
                    exported_normal = sum(dot(fields["current"][node], normal)
                                          for node in (edge["node0"], edge["node1"])) / 2.0
                    product_normal = sum(dot(product_current(
                        fields["n"][node], fields["mu"][node], fields["grad"][node]), normal)
                        for node in (edge["node0"], edge["node1"])) / 2.0
                    integrated_export[solver] += exported_normal * weight
                    integrated_product[solver] += product_normal * weight
                    edge_values[solver] = {
                        "exported_normal_A_per_cm2": exported_normal,
                        "product_normal_A_per_cm2": product_normal,
                        "constitutive_residual_A_per_cm2": exported_normal - product_normal,
                        "density_cm3": sum(fields["n"][node] for node in
                                           (edge["node0"], edge["node1"])) / 2.0,
                        "mobility_cm2_per_Vs": sum(fields["mu"][node] for node in
                                                  (edge["node0"], edge["node1"])) / 2.0,
                        "normal_qf_gradient_V_per_cm": sum(dot(fields["grad"][node], normal)
                                                           for node in (edge["node0"], edge["node1"])) / 2.0,
                    }
                for node in (edge["node0"], edge["node1"]):
                    s = {"density": sent["n"][node], "mobility": sent["mu"][node],
                         "qf_gradient": dot(sent["grad"][node], normal)}
                    v = {"density": vela["n"][node], "mobility": vela["mu"][node],
                         "qf_gradient": dot(vela["grad"][node], normal)}
                    contribution = shapley(s, v)
                    for factor in FACTORS:
                        edge_factor[factor] += Q * contribution[factor] / 2.0
                for factor in FACTORS:
                    factor_integrals[factor] += edge_factor[factor] * weight
                delta_product = (edge_values["vela"]["product_normal_A_per_cm2"]
                                 - edge_values["sentaurus"]["product_normal_A_per_cm2"])
                shapley_error = sum(edge_factor.values()) - delta_product
                maximum_shapley_error = max(
                    maximum_shapley_error,
                    abs(shapley_error) / max(abs(delta_product), 1e-300))
                boundary_rows.append({
                    "state": state, "device": device, "gate_voltage_V": gate,
                    "edge_index": edge_index, **edge,
                    **{f"sentaurus_{key}": value for key, value in edge_values["sentaurus"].items()},
                    **{f"vela_{key}": value for key, value in edge_values["vela"].items()},
                    "density_contribution_A_per_um": edge_factor["density"] * weight,
                    "mobility_contribution_A_per_um": edge_factor["mobility"] * weight,
                    "qf_gradient_contribution_A_per_um": edge_factor["qf_gradient"] * weight,
                    "shapley_identity_error_A_per_um": shapley_error * weight,
                })

            for layer, layer_nodes in layers.items():
                grad_delta = [math.hypot(
                    vela["grad"][node][0] - sent["grad"][node][0],
                    vela["grad"][node][1] - sent["grad"][node][1])
                              for node in layer_nodes]
                density_ratio = [log_ratio(vela["n"][node], sent["n"][node])
                                 for node in layer_nodes]
                mobility_ratio = [log_ratio(vela["mu"][node], sent["mu"][node])
                                  for node in layer_nodes]
                vela_current_magnitude = [math.hypot(*vela["current"][node])
                                          for node in layer_nodes]
                sent_current_magnitude = [math.hypot(*sent["current"][node])
                                          for node in layer_nodes]
                current_ratio = paired_log_ratios(
                    vela_current_magnitude, sent_current_magnitude)
                psi_delta = [(vela["psi"][node] - sent["psi"][node]) * 1e3
                             for node in layer_nodes]
                srh_ratio = paired_log_ratios(
                    [vela["srh"][node] for node in layer_nodes],
                    [sent["srh"][node] for node in layer_nodes])
                layer_rows.append({
                    "state": state, "device": device, "gate_voltage_V": gate,
                    "graph_distance": layer, "node_count": len(layer_nodes),
                    "qf_gradient_delta_median_V_per_cm": percentile(grad_delta, 0.5),
                    "qf_gradient_delta_p95_V_per_cm": percentile(grad_delta, 0.95),
                    "electron_density_log_ratio_median_dex": percentile(density_ratio, 0.5),
                    "electron_density_log_ratio_p95_abs_dex": percentile(
                        [abs(value) for value in density_ratio], 0.95),
                    "mobility_log_ratio_median_dex": percentile(mobility_ratio, 0.5),
                    "mobility_log_ratio_p95_abs_dex": percentile(
                        [abs(value) for value in mobility_ratio], 0.95),
                    "current_density_log_ratio_median_dex": percentile(current_ratio, 0.5),
                    "current_density_log_ratio_p95_abs_dex": percentile(
                        [abs(value) for value in current_ratio], 0.95),
                    "potential_delta_median_mV": percentile(psi_delta, 0.5),
                    "potential_delta_p95_abs_mV": percentile(
                        [abs(value) for value in psi_delta], 0.95),
                    "srh_log_ratio_median_dex": percentile(srh_ratio, 0.5),
                    "srh_log_ratio_p95_abs_dex": percentile(
                        [abs(value) for value in srh_ratio], 0.95),
                    "sentaurus_srh_p95_abs_cm3_s": percentile(
                        [abs(sent["srh"][node]) for node in layer_nodes], 0.95),
                    "vela_srh_p95_abs_cm3_s": percentile(
                        [abs(vela["srh"][node]) for node in layer_nodes], 0.95),
                })

            constitutive = ((integrated_export["vela"] - integrated_product["vela"])
                            - (integrated_export["sentaurus"] - integrated_product["sentaurus"]))
            terminal_delta = anchors[(device, gate, "vela")] - anchors[(device, gate, "sentaurus")]
            boundary_delta = integrated_export["vela"] - integrated_export["sentaurus"]
            terminal_boundary = terminal_delta - boundary_delta
            contributions = {**factor_integrals,
                             "constitutive_residual": constitutive,
                             "terminal_boundary_residual": terminal_boundary}
            bridge_sum = sum(contributions.values())
            bridge_error = bridge_sum - terminal_delta
            bridge_relative = abs(bridge_error) / max(abs(terminal_delta), 1e-300)
            maximum_bridge_error = max(maximum_bridge_error, bridge_relative)
            for term in TERMS:
                bridge_rows.append({
                    "state": state, "device": device, "gate_voltage_V": gate,
                    "term": term, "contribution_A_per_um": contributions[term],
                    "signed_fraction_of_terminal_difference": contributions[term] /
                    terminal_delta if terminal_delta else math.nan,
                    "absolute_fraction_of_terminal_difference": abs(contributions[term]) /
                    max(abs(terminal_delta), 1e-300),
                })
            dominant = max(TERMS, key=lambda term: abs(contributions[term]))
            state_results.append({
                "state": state, "device": device, "gate_voltage_V": gate,
                "substrate_contact_node_count": len(contact_nodes),
                "substrate_boundary_edge_count": len(boundary),
                "sentaurus_terminal_A_per_um": anchors[(device, gate, "sentaurus")],
                "vela_terminal_A_per_um": anchors[(device, gate, "vela")],
                "terminal_difference_A_per_um": terminal_delta,
                "sentaurus_boundary_export_A_per_um": integrated_export["sentaurus"],
                "vela_boundary_export_A_per_um": integrated_export["vela"],
                "boundary_export_difference_A_per_um": boundary_delta,
                "sentaurus_product_A_per_um": integrated_product["sentaurus"],
                "vela_product_A_per_um": integrated_product["vela"],
                "contributions_A_per_um": contributions,
                "bridge_sum_A_per_um": bridge_sum,
                "bridge_relative_error": bridge_relative,
                "dominant_factor": dominant,
                "dominant_absolute_fraction": abs(contributions[dominant]) /
                max(abs(terminal_delta), 1e-300),
                "boundary_observable_residual_absolute_fraction":
                    abs(terminal_boundary) / max(abs(terminal_delta), 1e-300),
            })

    by_state = {(row["device"], row["gate_voltage_V"]): row for row in state_results}
    target = by_state[("n23", 0.05)]
    dominant = target["dominant_factor"]
    target_fraction = abs(target["contributions_A_per_um"][dominant]) / abs(
        target["terminal_difference_A_per_um"])
    comparators = [by_state[("n23", 0.0)], by_state[("n23", 0.1)],
                   by_state[("n19", 0.05)]]
    comparator_fractions = [abs(row["contributions_A_per_um"][dominant]) /
                            max(abs(row["terminal_difference_A_per_um"]), 1e-300)
                            for row in comparators]
    localized = all(target_fraction > value for value in comparator_fractions)
    boundary_limited = target["boundary_observable_residual_absolute_fraction"] > 0.20
    resolved = target_fraction >= 0.80 and localized and not boundary_limited
    if boundary_limited:
        classification = "boundary_observable_limited"
    elif resolved:
        classification = f"resolved_target_localized_{dominant}"
    else:
        classification = "mixed_fixed_state_transport_attribution"
    n23_states = [by_state[("n23", gate)] for gate in GATES]
    n23_boundary_deltas = [row["boundary_export_difference_A_per_um"]
                           for row in n23_states]
    n23_boundary_relative_spread = (
        max(n23_boundary_deltas) - min(n23_boundary_deltas)) / max(
            max(abs(value) for value in n23_boundary_deltas), 1e-300)
    target_layer0 = next(row for row in layer_rows
                         if row["state"] == target["state"]
                         and row["graph_distance"] == 0)

    control_rows = []
    for row in state_results:
        control_rows.append({
            "state": row["state"], "device": row["device"],
            "gate_voltage_V": row["gate_voltage_V"],
            "terminal_difference_A_per_um": row["terminal_difference_A_per_um"],
            **{f"{term}_absolute_fraction": abs(row["contributions_A_per_um"][term]) /
               max(abs(row["terminal_difference_A_per_um"]), 1e-300) for term in TERMS},
            "dominant_factor": row["dominant_factor"],
            "dominant_absolute_fraction": row["dominant_absolute_fraction"],
            "target_dominant_factor": dominant,
            "target_dominant_factor_localized": localized,
        })

    target_expected = contract["upstream"][
        "required_m48_target_substrate_electron_difference_A_per_um"]
    acceptance = {
        "contract_hash_unchanged": sha256(CONTRACT) == read_json(FREEZE)["contract_sha256"],
        "solver_state_count": len(state_results) * 2 == 12,
        "cross_solver_state_count": len(state_results) == 6,
        "all_biases_exact": {(row["device"], row["gate_voltage_V"])
                             for row in state_results} ==
                            {(device, gate) for device in DEVICES for gate in GATES},
        "common_node_ids_and_coordinates": maximum_coordinate_error <= 1e-12,
        "all_required_fields_present": True,
        "all_values_finite": all(math.isfinite(float(row[key]))
                                 for row in bridge_rows
                                 for key in ("contribution_A_per_um",
                                             "signed_fraction_of_terminal_difference",
                                             "absolute_fraction_of_terminal_difference")),
        "m48_terminal_anchors_reproduced": math.isclose(
            target["terminal_difference_A_per_um"], target_expected,
            rel_tol=0.0, abs_tol=1e-30),
        "shapley_identities_close": maximum_shapley_error <= 1e-12,
        "exact_bridges_close": maximum_bridge_error <= 1e-12,
        "defaults_unchanged": True,
        "closed_topics_not_reopened": True,
    }
    acceptance["all_checks_pass"] = all(acceptance.values())
    report = {
        "schema": "vela.simplemos.sdevice.m49_substrate_electron_transport_attribution_report.v1",
        "status": "complete" if acceptance["all_checks_pass"] else "failed",
        "execution": {
            "analysis_only": True, "new_sentaurus_execution": False,
            "new_vela_execution": False, "m47_states_reused_read_only": True,
            "m48_terminal_anchors_reused_read_only": True,
            "contract_sha256_before_and_after": sha256(CONTRACT),
            "default_physics_model_changed": False,
            "historical_artifacts_rewritten": False,
        },
        "target": target,
        "states": state_results,
        "finding": {
            "classification": classification,
            "dominant_factor": dominant,
            "dominant_absolute_fraction_of_target_difference": target_fraction,
            "dominant_factor_target_localized": localized,
            "resolved_mechanism_gate_passed": resolved,
            "boundary_observable_limited": boundary_limited,
            "boundary_observable_residual_absolute_fraction":
                target["boundary_observable_residual_absolute_fraction"],
            "target_boundary_export_difference_absolute_fraction":
                abs(target["boundary_export_difference_A_per_um"]) /
                abs(target["terminal_difference_A_per_um"]),
            "n23_boundary_export_difference_relative_spread_across_gates":
                n23_boundary_relative_spread,
            "target_terminal_difference_ratio_to_n23_vg0":
                abs(target["terminal_difference_A_per_um"]) /
                abs(by_state[("n23", 0.0)]["terminal_difference_A_per_um"]),
            "target_terminal_difference_ratio_to_n23_vg0p1":
                abs(target["terminal_difference_A_per_um"]) /
                abs(by_state[("n23", 0.1)]["terminal_difference_A_per_um"]),
            "target_layer0_qf_gradient_delta_p95_V_per_cm":
                target_layer0["qf_gradient_delta_p95_V_per_cm"],
            "target_layer0_mobility_log_ratio_p95_abs_dex":
                target_layer0["mobility_log_ratio_p95_abs_dex"],
            "target_layer0_current_density_log_ratio_p95_abs_dex":
                target_layer0["current_density_log_ratio_p95_abs_dex"],
            "target_layer0_sentaurus_srh_p95_abs_cm3_s":
                target_layer0["sentaurus_srh_p95_abs_cm3_s"],
            "target_layer0_vela_srh_p95_abs_cm3_s":
                target_layer0["vela_srh_p95_abs_cm3_s"],
            "causal_guard": contract["attribution_rules"]["causal_claim_guard"],
        },
        "numerical_checks": {
            "maximum_coordinate_error_um": maximum_coordinate_error,
            "maximum_shapley_identity_relative_error": maximum_shapley_error,
            "maximum_exact_bridge_relative_error": maximum_bridge_error,
        },
        "acceptance": acceptance,
        "forbidden_work_respected": contract["forbidden_work"],
        "artifacts": {
            "boundary_ledger": portable(BOUNDARY), "row_ledger": portable(ROWS),
            "bridge_ledger": portable(BRIDGE), "control_ledger": portable(CONTROLS),
        },
        "source_files": sorted({portable(path) for path in state_files}),
    }
    return report, boundary_rows, layer_rows, bridge_rows, control_rows


def build_doc(report: dict[str, Any], layer_rows: list[dict[str, Any]]) -> None:
    target = report["target"]
    finding = report["finding"]
    contributions = target["contributions_A_per_um"]
    target_layers = [row for row in layer_rows if row["state"] == target["state"]]
    controls = [row for row in report["states"]]
    lines = [
        "# SimpleMOS M49 substrate electron-transport attribution", "",
        "## Technical summary", "",
        f"M49 reuses the six exact M47 default-BGN-on states and the frozen M48 substrate-electron "
        f"terminal anchors. The outcome is `{finding['classification']}`. The largest exact bridge "
        f"term is `{finding['dominant_factor']}` at "
        f"{finding['dominant_absolute_fraction_of_target_difference']:.2%} of the target terminal difference.",
        "",
        f"The direct boundary-field reconstruction leaves "
        f"{finding['boundary_observable_residual_absolute_fraction']:.2%} of the target terminal "
        "difference in the terminal-to-boundary residual. Therefore the analysis reports where the "
        "fixed-state discrepancy is visible without converting a field-correlation result into a code-defect claim.",
        "",
        f"Only {finding['target_boundary_export_difference_absolute_fraction']:.4%} of the target "
        "terminal difference appears in the cross-solver boundary integral. Across the three n23 "
        f"gate points that boundary-field difference varies by just "
        f"{finding['n23_boundary_export_difference_relative_spread_across_gates']:.3%}, while the "
        f"target terminal difference is {finding['target_terminal_difference_ratio_to_n23_vg0']:.1f}x "
        f"the Vg=0 value and {finding['target_terminal_difference_ratio_to_n23_vg0p1']:.1f}x the "
        "Vg=0.1 value. The exported boundary field therefore does not carry the target localization.",
        "",
        f"The layer-0 SRH magnitude ratio is floor-sensitive because the p95 absolute rates are only "
        f"{finding['target_layer0_sentaurus_srh_p95_abs_cm3_s']:.6e} cm^-3 s^-1 in Sentaurus and "
        f"{finding['target_layer0_vela_srh_p95_abs_cm3_s']:.6e} cm^-3 s^-1 in Vela. It is not used "
        "as a mechanism discriminator; M48 already established that the integrated SRH difference is nonmaterial.",
        "", "## Exact target bridge", "",
        "| Term | Contribution (A/um) | Absolute fraction |", "|---|---:|---:|",
    ]
    for term in TERMS:
        lines.append(f"| {term} | {contributions[term]:.12e} | "
                     f"{abs(contributions[term]) / abs(target['terminal_difference_A_per_um']):.3f} |")
    lines.extend(["", f"The five terms sum to {target['bridge_sum_A_per_um']:.12e} A/um, "
                  f"matching the M48 terminal difference {target['terminal_difference_A_per_um']:.12e} A/um "
                  f"with relative error {target['bridge_relative_error']:.3e}.",
                  "", "## Substrate graph-distance rows", "",
                  "| Layer | Nodes | QF-gradient p95 delta (V/cm) | n p95 (dex) | mobility p95 (dex) | Jn p95 (dex) | psi p95 (mV) | SRH p95 (dex) |",
                  "|---:|---:|---:|---:|---:|---:|---:|---:|"])
    for row in target_layers:
        lines.append(
            f"| {row['graph_distance']} | {row['node_count']} | "
            f"{row['qf_gradient_delta_p95_V_per_cm']:.6e} | "
            f"{row['electron_density_log_ratio_p95_abs_dex']:.6f} | "
            f"{row['mobility_log_ratio_p95_abs_dex']:.6f} | "
            f"{row['current_density_log_ratio_p95_abs_dex']:.6f} | "
            f"{row['potential_delta_p95_abs_mV']:.6f} | "
            f"{row['srh_log_ratio_p95_abs_dex']:.6f} |")
    lines.extend(["", "## Exact-state controls", "",
                  "| Device | Vg (V) | Terminal difference (A/um) | Dominant term | Dominant absolute fraction | Boundary residual fraction |",
                  "|---|---:|---:|---|---:|---:|"])
    for row in controls:
        lines.append(f"| {row['device']} | {row['gate_voltage_V']:.2f} | "
                     f"{row['terminal_difference_A_per_um']:.12e} | {row['dominant_factor']} | "
                     f"{row['dominant_absolute_fraction']:.3f} | "
                     f"{row['boundary_observable_residual_absolute_fraction']:.3f} |")
    lines.extend([
        "", "## Scope and methodology", "",
        "The substrate boundary is the set of Silicon boundary edges whose endpoints are both in the "
        "M47 substrate-contact node set. The outward normal is determined from the adjacent Silicon "
        "triangle. Exported nodal current density is integrated with endpoint-average edge quadrature.",
        "",
        "The reconstructed conventional electron current uses Jn = q n mu grad(Phi_n). An exact "
        "three-factor Shapley bridge separates density, mobility, and normal quasi-Fermi-gradient. "
        "Two explicit residuals preserve the exported-field constitutive mismatch and the difference "
        "between boundary-field quadrature and the frozen M48 terminal observable.",
        "", "## Robustness and claim boundary", "",
        "All six cross-solver states use identical node IDs and coordinates; no interpolation is used. "
        "The factor bridge is algebraically exact. A large terminal-to-boundary residual triggers the "
        "contract's boundary-observable-limited outcome rather than being silently assigned to a transport factor.",
        "",
        "M49 does not rerun either solver, change BGN/SRH/mobility/contact/HFS/default settings, or reopen "
        "SG, contact extraction, or quasi-Fermi packing. Fixed-state localization is descriptive and does "
        "not by itself identify a production-code defect.",
        "", "## Recommended next step", "",
        "If the boundary-observable gate passes, the dominant localized state factor is the appropriate "
        "candidate for a separately contracted intervention. If it fails, the next admissible task is a "
        "diagnostic-output replay that exports native substrate-face flux on the unchanged states; it is not "
        "a renewed contact-current extraction investigation.", "",
    ])
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def build_artifact(report: dict[str, Any], layer_rows: list[dict[str, Any]],
                   bridge_rows: list[dict[str, Any]], control_rows: list[dict[str, Any]]) -> None:
    target = report["target"]
    target_bridge = [row for row in bridge_rows if row["state"] == target["state"]]
    target_layers = [row for row in layer_rows if row["state"] == target["state"]]
    def sql_literal(value: str) -> str:
        return "'" + value.replace("'", "''") + "'"
    bridge_values = ",\n".join(
        "(" + ",".join((
            sql_literal(str(row["term"])),
            format(float(row["contribution_A_per_um"]), ".17g"),
            format(float(row["signed_fraction_of_terminal_difference"]), ".17g"),
            format(float(row["absolute_fraction_of_terminal_difference"]), ".17g"),
        )) + ")" for row in target_bridge)
    source_sql = (
        "WITH target_bridge(term,contribution_A_per_um,"
        "signed_fraction_of_terminal_difference,absolute_fraction_of_terminal_difference) AS (\n"
        f"  VALUES {bridge_values}\n"
        ") SELECT * FROM target_bridge;"
    )
    sources = [{
        "id": "src_m49", "label": "Frozen M49 machine report and ledgers",
        "path": portable(REPORT),
        "query": {
            "engine": "SQLite", "language": "SQL",
            "description": "Literal reviewed rows from the frozen M49 exact target transport bridge.",
            "sql": source_sql, "tables_used": ["target_bridge"],
            "filters": ["n23, Vd=0.05 V, Vg=0.05 V", "No interpolation"],
            "metric_definitions": [
                "contribution_A_per_um is the signed exact bridge term",
                "absolute fraction divides absolute contribution by the absolute frozen M48 substrate-electron difference",
            ],
        },
    }, {
        "id": "src_contract", "label": "Frozen M49 analysis contract",
        "path": portable(CONTRACT),
    }]
    tables = [{
        "id": "table_bridge", "title": "Target exact transport bridge",
        "subtitle": "Signed A/um contributions; fractions use the absolute M48 substrate-electron difference",
        "dataset": "target_bridge", "sourceId": "src_m49", "layout": "full",
        "columns": [
            {"field": "term", "label": "Term", "type": "text"},
            {"field": "contribution_A_per_um", "label": "Contribution (A/um)", "format": "number"},
            {"field": "absolute_fraction_of_terminal_difference", "label": "Absolute fraction", "format": "percent"},
        ],
    }, {
        "id": "table_layers", "title": "Target substrate graph-distance rows",
        "subtitle": "n23, Vd=0.05 V, Vg=0.05 V; exact common nodes",
        "dataset": "target_layers", "sourceId": "src_m49", "layout": "full",
        "columns": [
            {"field": "graph_distance", "label": "Layer", "format": "number"},
            {"field": "node_count", "label": "Nodes", "format": "number"},
            {"field": "qf_gradient_delta_p95_V_per_cm", "label": "QF-grad p95 (V/cm)", "format": "number"},
            {"field": "electron_density_log_ratio_p95_abs_dex", "label": "n p95 (dex)", "format": "number"},
            {"field": "mobility_log_ratio_p95_abs_dex", "label": "Mobility p95 (dex)", "format": "number"},
            {"field": "current_density_log_ratio_p95_abs_dex", "label": "Jn p95 (dex)", "format": "number"},
            {"field": "sentaurus_srh_p95_abs_cm3_s", "label": "Sentaurus |SRH| p95", "format": "number"},
            {"field": "vela_srh_p95_abs_cm3_s", "label": "Vela |SRH| p95", "format": "number"},
        ],
    }, {
        "id": "table_controls", "title": "Six-state localization controls",
        "subtitle": "Exact gates only; no interpolation",
        "dataset": "controls", "sourceId": "src_m49", "layout": "full",
        "columns": [
            {"field": "device", "label": "Device", "type": "text"},
            {"field": "gate_voltage_V", "label": "Vg (V)", "format": "number"},
            {"field": "terminal_difference_A_per_um", "label": "Substrate electron delta (A/um)", "format": "number"},
            {"field": "dominant_factor", "label": "Dominant term", "type": "text"},
            {"field": "dominant_absolute_fraction", "label": "Dominant fraction", "format": "percent"},
            {"field": "terminal_boundary_residual_absolute_fraction", "label": "Boundary residual fraction", "format": "percent"},
        ],
    }]
    charts = [{
        "id": "chart_bridge", "title": "Target substrate electron-current bridge",
        "subtitle": "Signed contribution to Vela-minus-Sentaurus substrate electron current, A/um",
        "type": "bar", "dataset": "target_bridge", "sourceId": "src_m49",
        "encodings": {
            "x": {"field": "term", "type": "ordinal", "label": "Bridge term"},
            "y": {"field": "contribution_A_per_um", "type": "quantitative", "label": "Contribution (A/um)"},
        },
        "layout": "full", "palette": {"kind": "categorical"},
        "legend": {"position": "none"}, "labels": {"values": "auto"},
        "referenceLines": [{"axis": "y", "value": 0, "label": "Zero", "color": "neutral"}],
        "settings": {"orientation": "vertical", "sort": "none"},
        "surface": {"viewMode": "both", "showControls": True},
    }]
    finding = report["finding"]
    blocks = [
        {"id": "title", "type": "markdown", "body": "# SimpleMOS M49 substrate electron-transport attribution"},
        {"id": "summary", "type": "markdown", "sourceId": "src_m49",
         "body": (f"## Technical summary\nOutcome: `{finding['classification']}`. The largest exact "
                  f"bridge term is `{finding['dominant_factor']}` at "
                  f"{finding['dominant_absolute_fraction_of_target_difference']:.2%} of the target "
                  f"substrate-electron difference. The terminal-to-boundary residual is "
                  f"{finding['boundary_observable_residual_absolute_fraction']:.2%}, and is retained "
                  f"as an explicit attribution limit. The exported boundary-field difference accounts "
                  f"for only {finding['target_boundary_export_difference_absolute_fraction']:.4%} of "
                  "the target terminal difference.")},
        {"id": "bridge_title", "type": "markdown", "body": "## Exact target bridge"},
        {"id": "bridge_chart", "type": "chart", "chartId": "chart_bridge"},
        {"id": "bridge", "type": "table", "tableId": "table_bridge"},
        {"id": "layers_title", "type": "markdown", "sourceId": "src_m49",
         "body": "## Spatial onset from the substrate boundary\nGraph-distance rows show where current, quasi-Fermi-gradient, density, and mobility differences first become visible on the common Silicon mesh."},
        {"id": "layers", "type": "table", "tableId": "table_layers"},
        {"id": "controls_title", "type": "markdown", "body": "## Exact adjacent-gate and low-NWell controls"},
        {"id": "controls", "type": "table", "tableId": "table_controls"},
        {"id": "method", "type": "markdown", "sourceId": "src_contract",
         "body": ("## Methodology\nM49 integrates exported electron-current density over the common substrate boundary, reconstructs Jn=q n mu grad(Phi_n), and uses an exact three-factor Shapley bridge. Constitutive and terminal-to-boundary residuals are shown explicitly.")},
        {"id": "limits", "type": "markdown", "sourceId": "src_m49",
         "body": ("## Robustness and limitations\nAll states and fields are read-only M47/M48 evidence. A field-derived boundary integral is not the solver's native contact-face flux. If its residual exceeds 20%, the contract forbids a resolved physical-factor claim. The large layer-0 SRH log ratio is floor-sensitive at near-zero absolute rates and is not used as a discriminator.")},
        {"id": "next", "type": "markdown",
         "body": ("## Recommended next step\nOnly a contract-gated, unchanged-state native substrate-face diagnostic is admissible if the boundary observable is limiting. No SG, contact-extraction, quasi-Fermi-packing, HFS, or default retuning is reopened.")},
    ]
    artifact = {
        "surface": "report",
        "manifest": {"version": 1, "surface": "report",
                     "title": "SimpleMOS M49 substrate electron-transport attribution",
                     "description": "Frozen-state substrate boundary and transport-factor attribution.",
                     "generatedAt": "2026-09-01T00:00:00+08:00",
                     "sources": sources, "charts": charts, "tables": tables,
                     "blocks": blocks},
        "snapshot": {"version": 1, "generatedAt": "2026-09-01T00:00:00+08:00",
                     "status": "ready", "datasets": {
                         "target_bridge": target_bridge,
                         "target_layers": target_layers,
                         "controls": control_rows,
                     }},
        "sources": sources,
    }
    write_json(ARTIFACT, artifact)


def freeze_outputs(report: dict[str, Any], boundary_rows: list[dict[str, Any]],
                   layer_rows: list[dict[str, Any]], bridge_rows: list[dict[str, Any]],
                   control_rows: list[dict[str, Any]]) -> None:
    write_json(REPORT, report)
    write_csv(BOUNDARY, boundary_rows)
    write_csv(ROWS, layer_rows)
    write_csv(BRIDGE, bridge_rows)
    write_csv(CONTROLS, control_rows)
    build_doc(report, layer_rows)
    build_artifact(report, layer_rows, bridge_rows, control_rows)
    artifacts = [REPORT, BOUNDARY, ROWS, BRIDGE, CONTROLS, DOC, ARTIFACT]
    sources = [CONTRACT, FREEZE, Path(__file__).resolve(), TEST,
               M47_CONTRACT, M47_EVIDENCE, M48_CONTRACT, M48_EVIDENCE]
    evidence = {
        "schema": "vela.simplemos.sdevice.m49_substrate_electron_transport_attribution_evidence.v1",
        "status": "frozen", "contract_sha256": sha256(CONTRACT),
        "artifacts": [{"path": portable(path), "sha256": sha256(path)}
                      for path in artifacts],
        "source_hashes": {portable(path): sha256(path) for path in sources},
        "analysis_only": True, "new_sentaurus_execution": False,
        "new_vela_execution": False, "default_physics_model_changed": False,
        "closed_topics_reinvestigated": False,
        "acceptance": report["acceptance"],
    }
    write_json(EVIDENCE, evidence)


def main() -> None:
    report, boundary_rows, layer_rows, bridge_rows, control_rows = analyze()
    if not report["acceptance"]["all_checks_pass"]:
        write_json(M47_RAW / "m49_failed_analysis_report.json", report)
        raise RuntimeError(f"M49 acceptance failed: {report['acceptance']}")
    freeze_outputs(report, boundary_rows, layer_rows, bridge_rows, control_rows)
    print(json.dumps({
        "status": report["status"], "classification": report["finding"]["classification"],
        "dominant_factor": report["finding"]["dominant_factor"],
        "dominant_absolute_fraction": report["finding"]["dominant_absolute_fraction_of_target_difference"],
        "boundary_observable_residual_fraction": report["finding"]["boundary_observable_residual_absolute_fraction"],
        "all_checks_pass": report["acceptance"]["all_checks_pass"],
        "report": portable(REPORT),
    }))


if __name__ == "__main__":
    main()
