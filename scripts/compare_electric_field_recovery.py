#!/usr/bin/env python3
"""Compare P1 Tri3 electric-field recovery methods against Sentaurus node fields."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any


METHODS = ["area", "ls1d", "ls1d2", "spr"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--bias", type=float, default=0.0)
    parser.add_argument("--out-dir", type=Path, default=Path("build/diagnostics"))
    parser.add_argument("--potential-field", default="ElectrostaticPotential,Potential")
    parser.add_argument("--electric-field", default="ElectricField")
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def parse_node_list(raw: str) -> list[int]:
    return [int(item) for item in raw.replace(";", ",").replace(" ", ",").split(",") if item]


def read_nodes(root: Path) -> dict[int, tuple[float, float]]:
    nodes: dict[int, tuple[float, float]] = {}
    for row in read_csv(root / "nodes.csv"):
        nodes[int(row["id"])] = (float(row["x_um"]), float(row["y_um"]))
    return nodes


def read_elements(root: Path) -> list[dict[str, Any]]:
    elements: list[dict[str, Any]] = []
    region_ids: dict[str, int] = {}
    for row in read_csv(root / "elements.csv"):
        region_name = row.get("region", "region0") or "region0"
        if region_name not in region_ids:
            region_ids[region_name] = len(region_ids)
        elements.append({
            "id": int(row.get("id", len(elements))),
            "nodes": [int(row["node0"]), int(row["node1"]), int(row["node2"])],
            "region": region_ids[region_name],
            "region_name": region_name,
            "material": row.get("material", ""),
        })
    return elements


def read_contacts(root: Path) -> set[int]:
    path = root / "contacts.csv"
    if not path.exists():
        return set()
    nodes: set[int] = set()
    for row in read_csv(path):
        nodes.update(parse_node_list(row.get("node_ids", "")))
    return nodes


def field_paths(root: Path, field_names: list[str]) -> list[Path]:
    paths: list[Path] = []
    for name in field_names:
        paths.extend(sorted((root / "fields").glob(f"{name}_region*.csv")))
        paths.extend(sorted(root.glob(f"{name}_region*.csv")))
    return paths


def region_from_path(path: Path) -> int:
    stem = path.stem
    marker = "_region"
    if marker not in stem:
        return 0
    return int(stem.rsplit(marker, 1)[1])


def read_scalar_field(root: Path, names: list[str]) -> dict[int, float]:
    paths = field_paths(root, names)
    if not paths:
        raise FileNotFoundError(f"missing scalar field {names} under {root}")
    values: dict[int, float] = {}
    for path in paths:
        for row in read_csv(path):
            cols = [name for name in row if name != "node_id"]
            if not cols:
                continue
            values[int(row["node_id"])] = float(row[cols[0]])
    if not values:
        raise FileNotFoundError(f"empty scalar field {names} under {root}")
    return values


def read_sentaurus_electric(root: Path, names: list[str]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in field_paths(root, names):
        region = region_from_path(path)
        for row in read_csv(path):
            cols = [name for name in row if name != "node_id"]
            comps = [float(row[name]) for name in cols if row.get(name, "") != ""]
            ex = comps[0] if len(comps) >= 2 else math.nan
            ey = comps[1] if len(comps) >= 2 else math.nan
            emag = math.sqrt(ex * ex + ey * ey) if len(comps) >= 2 else (comps[0] if comps else math.nan)
            rows.append({"node_id": int(row["node_id"]), "region": region, "ex": ex, "ey": ey, "emag": emag})
    if rows:
        return rows
    raise FileNotFoundError(f"missing electric field {names} under {root}")


def cell_gradient(nodes: dict[int, tuple[float, float]], tri: list[int], values: dict[int, float]) -> tuple[float, float, float] | None:
    n0, n1, n2 = tri
    x0, y0 = nodes[n0]
    x1, y1 = nodes[n1]
    x2, y2 = nodes[n2]
    dx1, dy1 = x1 - x0, y1 - y0
    dx2, dy2 = x2 - x0, y2 - y0
    det = dx1 * dy2 - dy1 * dx2
    if abs(det) <= 1.0e-30:
        return None
    dv1 = values[n1] - values[n0]
    dv2 = values[n2] - values[n0]
    gx = (dv1 * dy2 - dy1 * dv2) / det
    gy = (dx1 * dv2 - dv1 * dx2) / det
    area = 0.5 * abs(det)
    return gx, gy, area


def build_patch(elements: list[dict[str, Any]]) -> tuple[dict[int, list[int]], dict[tuple[int, int], set[int]], dict[int, set[int]]]:
    node_cells: dict[int, list[int]] = defaultdict(list)
    region_neighbors: dict[tuple[int, int], set[int]] = defaultdict(set)
    node_regions: dict[int, set[int]] = defaultdict(set)
    for index, cell in enumerate(elements):
        region = int(cell["region"])
        tri = list(cell["nodes"])
        for node in tri:
            node_cells[node].append(index)
            node_regions[node].add(region)
            for other in tri:
                if other != node:
                    region_neighbors[(node, region)].add(other)
    return node_cells, region_neighbors, node_regions



def solve3(matrix: list[list[float]], rhs: list[float]) -> list[float] | None:
    a = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(3):
        pivot = max(range(col, 3), key=lambda row: abs(a[row][col]))
        if abs(a[pivot][col]) <= 1.0e-30:
            return None
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
        scale = a[col][col]
        for j in range(col, 4):
            a[col][j] /= scale
        for row in range(3):
            if row == col:
                continue
            factor = a[row][col]
            for j in range(col, 4):
                a[row][j] -= factor * a[col][j]
    return [a[i][3] for i in range(3)]

def recover_fields(nodes: dict[int, tuple[float, float]], elements: list[dict[str, Any]], potential: dict[int, float]) -> dict[str, dict[tuple[int, int], tuple[float, float]]]:
    node_cells, region_neighbors, node_regions = build_patch(elements)
    cell_fields: list[dict[str, Any]] = []
    for cell in elements:
        grad = cell_gradient(nodes, list(cell["nodes"]), potential)
        if grad is None:
            cell_fields.append({"valid": False})
            continue
        gx, gy, area = grad
        cell_fields.append({"valid": True, "region": int(cell["region"]), "nodes": list(cell["nodes"]), "ex": -gx, "ey": -gy, "area": area})

    area: dict[tuple[int, int], tuple[float, float]] = {}
    accum: dict[tuple[int, int], list[float]] = defaultdict(lambda: [0.0, 0.0, 0.0])
    for field in cell_fields:
        if not field.get("valid"):
            continue
        for node in field["nodes"]:
            key = (node, int(field["region"]))
            accum[key][0] += float(field["area"]) * float(field["ex"])
            accum[key][1] += float(field["area"]) * float(field["ey"])
            accum[key][2] += float(field["area"])
    for key, (sx, sy, sa) in accum.items():
        if sa > 0.0:
            area[key] = (sx / sa, sy / sa)

    def ls(node: int, region: int, power: int) -> tuple[float, float]:
        cx, cy = nodes[node]
        cv = potential[node]
        sxx = sxy = syy = sxv = syv = 0.0
        for other in region_neighbors.get((node, region), set()):
            ox, oy = nodes[other]
            dx, dy = ox - cx, oy - cy
            distance = math.hypot(dx, dy)
            if distance <= 1.0e-30:
                continue
            weight = 1.0 / (distance ** power)
            dv = potential[other] - cv
            sxx += weight * dx * dx
            sxy += weight * dx * dy
            syy += weight * dy * dy
            sxv += weight * dx * dv
            syv += weight * dy * dv
        det = sxx * syy - sxy * sxy
        if abs(det) <= 1.0e-30:
            return area.get((node, region), (0.0, 0.0))
        gx = (sxv * syy - syv * sxy) / det
        gy = (sxx * syv - sxy * sxv) / det
        return -gx, -gy

    ls1d = {(node, region): ls(node, region, 1) for node, regions in node_regions.items() for region in regions}
    ls1d2 = {(node, region): ls(node, region, 2) for node, regions in node_regions.items() for region in regions}

    def spr(node: int, region: int) -> tuple[float, float]:
        samples: list[tuple[float, float, float, float]] = []
        for cell_index in node_cells.get(node, []):
            cell = elements[cell_index]
            field = cell_fields[cell_index]
            if int(cell["region"]) != region or not field.get("valid"):
                continue
            xs = [nodes[n][0] for n in cell["nodes"]]
            ys = [nodes[n][1] for n in cell["nodes"]]
            samples.append((sum(xs) / 3.0, sum(ys) / 3.0, float(field["ex"]), float(field["ey"])))
        if len(samples) < 3:
            return ls1d.get((node, region), area.get((node, region), (0.0, 0.0)))
        normal = [[0.0, 0.0, 0.0] for _ in range(3)]
        rhs_ex = [0.0, 0.0, 0.0]
        rhs_ey = [0.0, 0.0, 0.0]
        for x, y, ex, ey in samples:
            row = [1.0, x, y]
            for i in range(3):
                rhs_ex[i] += row[i] * ex
                rhs_ey[i] += row[i] * ey
                for j in range(3):
                    normal[i][j] += row[i] * row[j]
        coeff_ex = solve3(normal, rhs_ex)
        coeff_ey = solve3(normal, rhs_ey)
        if coeff_ex is None or coeff_ey is None:
            return ls1d.get((node, region), area.get((node, region), (0.0, 0.0)))
        x0, y0 = nodes[node]
        return (
            coeff_ex[0] + coeff_ex[1] * x0 + coeff_ex[2] * y0,
            coeff_ey[0] + coeff_ey[1] * x0 + coeff_ey[2] * y0,
        )

    spr_fields = {(node, region): spr(node, region) for node, regions in node_regions.items() for region in regions}
    return {"area": area, "ls1d": ls1d, "ls1d2": ls1d2, "spr": spr_fields}


def classify_nodes(elements: list[dict[str, Any]], contacts: set[int]) -> dict[int, str]:
    edge_counts: dict[tuple[int, int], int] = defaultdict(int)
    for cell in elements:
        tri = list(cell["nodes"])
        for a, b in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])]:
            edge_counts[tuple(sorted((a, b)))] += 1
    boundary_nodes: set[int] = set()
    for (a, b), count in edge_counts.items():
        if count == 1:
            boundary_nodes.update([a, b])
    classes: dict[int, str] = {}
    for cell in elements:
        for node in cell["nodes"]:
            if node in contacts:
                classes[node] = "contact"
            elif node in boundary_nodes:
                classes[node] = "boundary"
            else:
                classes[node] = "interior"
    return classes


def error(sent: dict[str, Any], cand: tuple[float, float]) -> float:
    ex, ey = cand
    if math.isfinite(sent["ex"]) and math.isfinite(sent["ey"]):
        return math.hypot(ex - sent["ex"], ey - sent["ey"])
    return abs(math.hypot(ex, ey) - sent["emag"])


def median(values: list[float]) -> float:
    return statistics.median(values) if values else math.nan


def finite(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    return value



def write_cell_and_edge_outputs(out_dir: Path,
                                bias: float,
                                nodes: dict[int, tuple[float, float]],
                                elements: list[dict[str, Any]],
                                potential: dict[int, float],
                                electron_qf: dict[int, float] | None,
                                hole_qf: dict[int, float] | None) -> None:
    cell_path = out_dir / "electric_field_recovery_cells.csv"
    with cell_path.open("w", newline="") as handle:
        fields = [
            "bias_V", "cell_id", "node0", "node1", "node2", "region",
            "cell_Ex", "cell_Ey", "cell_Emag",
            "cell_grad_eQF_x", "cell_grad_eQF_y", "cell_grad_eQF_mag",
            "cell_grad_hQF_x", "cell_grad_hQF_y", "cell_grad_hQF_mag",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for cell in elements:
            tri = list(cell["nodes"])
            grad = cell_gradient(nodes, tri, potential)
            if grad is None:
                ex = ey = emag = math.nan
            else:
                gx, gy, _ = grad
                ex = -gx * 1.0e4
                ey = -gy * 1.0e4
                emag = math.hypot(ex, ey)
            egx = egy = egm = math.nan
            if electron_qf is not None:
                egrad = cell_gradient(nodes, tri, electron_qf)
                if egrad is not None:
                    egx, egy, _ = egrad
                    egx *= 1.0e4
                    egy *= 1.0e4
                    egm = math.hypot(egx, egy)
            hgx = hgy = hgm = math.nan
            if hole_qf is not None:
                hgrad = cell_gradient(nodes, tri, hole_qf)
                if hgrad is not None:
                    hgx, hgy, _ = hgrad
                    hgx *= 1.0e4
                    hgy *= 1.0e4
                    hgm = math.hypot(hgx, hgy)
            writer.writerow({
                "bias_V": bias,
                "cell_id": cell["id"],
                "node0": tri[0],
                "node1": tri[1],
                "node2": tri[2],
                "region": cell["region_name"],
                "cell_Ex": finite(ex),
                "cell_Ey": finite(ey),
                "cell_Emag": finite(emag),
                "cell_grad_eQF_x": finite(egx),
                "cell_grad_eQF_y": finite(egy),
                "cell_grad_eQF_mag": finite(egm),
                "cell_grad_hQF_x": finite(hgx),
                "cell_grad_hQF_y": finite(hgy),
                "cell_grad_hQF_mag": finite(hgm),
            })

    edge_path = out_dir / "electric_field_recovery_edges.csv"
    edges: set[tuple[int, int]] = set()
    for cell in elements:
        tri = list(cell["nodes"])
        edges.update(tuple(sorted(edge)) for edge in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])])
    with edge_path.open("w", newline="") as handle:
        fields = ["bias_V", "edge_id", "node0", "node1", "edge_Ex", "edge_Ey", "edge_projected", "edge_Emag"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for edge_id, (n0, n1) in enumerate(sorted(edges)):
            x0, y0 = nodes[n0]
            x1, y1 = nodes[n1]
            length = math.hypot(x1 - x0, y1 - y0)
            if length <= 1.0e-30:
                ex = ey = projected = emag = math.nan
            else:
                projected = -(potential[n1] - potential[n0]) / length * 1.0e4
                ex = projected * (x1 - x0) / length
                ey = projected * (y1 - y0) / length
                emag = abs(projected)
            writer.writerow({
                "bias_V": bias,
                "edge_id": edge_id,
                "node0": n0,
                "node1": n1,
                "edge_Ex": finite(ex),
                "edge_Ey": finite(ey),
                "edge_projected": finite(projected),
                "edge_Emag": finite(emag),
            })

def main() -> None:
    args = parse_args()
    nodes = read_nodes(args.sentaurus_root)
    elements = read_elements(args.sentaurus_root)
    contacts = read_contacts(args.sentaurus_root)
    potential = read_scalar_field(args.sentaurus_root, [item.strip() for item in args.potential_field.split(",") if item.strip()])
    try:
        electron_qf = read_scalar_field(args.sentaurus_root, ["eQuasiFermiPotential", "ElectronQuasiFermi"])
    except FileNotFoundError:
        electron_qf = None
    try:
        hole_qf = read_scalar_field(args.sentaurus_root, ["hQuasiFermiPotential", "HoleQuasiFermi"])
    except FileNotFoundError:
        hole_qf = None
    sentaurus_rows = read_sentaurus_electric(args.sentaurus_root, [item.strip() for item in args.electric_field.split(",") if item.strip()])
    recovered = recover_fields(nodes, elements, potential)
    classes = classify_nodes(elements, contacts)
    region_names = {int(cell["region"]): str(cell["region_name"]) for cell in elements}

    args.out_dir.mkdir(parents=True, exist_ok=True)
    compare_path = args.out_dir / "electric_field_recovery_compare.csv"
    summary_path = args.out_dir / "electric_field_recovery_compare_summary.md"
    write_cell_and_edge_outputs(args.out_dir, args.bias, nodes, elements, potential, electron_qf, hole_qf)

    fieldnames = [
        "bias_V", "node_id", "x_um", "y_um", "region", "node_class",
        "sentaurus_Ex", "sentaurus_Ey", "sentaurus_Emag",
        "area_Ex", "area_Ey", "area_Emag",
        "ls1d_Ex", "ls1d_Ey", "ls1d_Emag",
        "ls1d2_Ex", "ls1d2_Ey", "ls1d2_Emag",
        "spr_Ex", "spr_Ey", "spr_Emag",
        "error_area", "error_ls1d", "error_ls1d2", "error_spr",
    ]
    stats: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    with compare_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for sent in sentaurus_rows:
            node = int(sent["node_id"])
            region = int(sent["region"])
            x, y = nodes[node]
            row: dict[str, Any] = {
                "bias_V": args.bias,
                "node_id": node,
                "x_um": x,
                "y_um": y,
                "region": region_names.get(region, f"region{region}"),
                "node_class": classes.get(node, "boundary"),
                "sentaurus_Ex": sent["ex"],
                "sentaurus_Ey": sent["ey"],
                "sentaurus_Emag": sent["emag"],
            }
            for method in METHODS:
                ex, ey = recovered[method].get((node, region), recovered[method].get((node, 0), (math.nan, math.nan)))
                # Coordinates are in um, so recovered gradients are V/um. Convert to V/cm.
                ex *= 1.0e4
                ey *= 1.0e4
                emag = math.hypot(ex, ey) if math.isfinite(ex) and math.isfinite(ey) else math.nan
                row[f"{method}_Ex"] = ex
                row[f"{method}_Ey"] = ey
                row[f"{method}_Emag"] = emag
                err = error(sent, (ex, ey))
                row[f"error_{method}"] = err
                if math.isfinite(err):
                    stats[row["node_class"]][method].append(err)
            writer.writerow({name: finite(row.get(name, "")) for name in fieldnames})

    lines = [
        "# Electric Field Recovery Compare Summary",
        "",
        "RecoveredNodeElectricField is a post-processing field, not the conservative FVM edge flux.",
        "",
        f"Bias: {args.bias:g} V",
        "",
    ]
    for node_class in ["interior", "boundary", "contact"]:
        lines.append(f"## {node_class}")
        class_stats = stats.get(node_class, {})
        best_method = "none"
        best_mean = math.inf
        lines.append("| method | mean_error | median_error | max_error | samples |")
        lines.append("| --- | ---: | ---: | ---: | ---: |")
        for method in METHODS:
            values = class_stats.get(method, [])
            mean = sum(values) / len(values) if values else math.nan
            med = median(values)
            max_value = max(values) if values else math.nan
            if math.isfinite(mean) and mean < best_mean:
                best_mean = mean
                best_method = method
            lines.append(f"| {method} | {mean:.17g} | {med:.17g} | {max_value:.17g} | {len(values)} |")
        lines.append("")
        lines.append(f"Closest to Sentaurus by mean error: {best_method}")
        lines.append("")
    summary_path.write_text("\n".join(lines), encoding="ascii")


if __name__ == "__main__":
    main()
