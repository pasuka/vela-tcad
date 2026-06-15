#!/usr/bin/env python3
"""Report left-right symmetry of PN2D 0V electric-field magnitudes."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any

from PIL import Image, ImageDraw, ImageFont


REPORT_NAME = "pn2d_0v_electric_field_symmetry.json"
MARKDOWN_NAME = "pn2d_0v_electric_field_symmetry.md"
PAIR_CSV_NAME = "pn2d_0v_electric_field_symmetry_pairs.csv"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_nodes(path: Path) -> dict[int, tuple[float, float]]:
    nodes: dict[int, tuple[float, float]] = {}
    for row in read_csv_rows(path):
        nodes[int(float(row["id"]))] = (
            float(row.get("x_um", row.get("x", 0.0))),
            float(row.get("y_um", row.get("y", 0.0))),
        )
    return nodes


def read_field_nodes(path: Path) -> tuple[dict[int, float], dict[int, float]]:
    sentaurus: dict[int, float] = {}
    vela: dict[int, float] = {}
    for row in read_csv_rows(path):
        node_id = int(float(row["node_id"]))
        sentaurus[node_id] = float(row["sentaurus_E"])
        vela[node_id] = float(row["vela_E"])
    return sentaurus, vela


def mesh_triangles(mesh: dict[str, Any]) -> list[tuple[int, list[int]]]:
    triangles: list[tuple[int, list[int]]] = []
    for index, entry in enumerate(mesh.get("triangles", [])):
        ids = [int(node) for node in entry.get("node_ids", [])]
        if len(ids) == 3:
            triangles.append((int(entry.get("id", index)), ids))
    return triangles


def find_mirror_pairs(nodes: dict[int, tuple[float, float]],
                      axis_x_um: float,
                      tolerance_um: float) -> tuple[list[tuple[int, int, float]], list[int]]:
    ids = sorted(nodes)
    used: set[int] = set()
    pairs: list[tuple[int, int, float]] = []
    unmatched: list[int] = []
    for node_id in ids:
        if node_id in used:
            continue
        x, y = nodes[node_id]
        target = (2.0 * axis_x_um - x, y)
        mirror_id = min(
            ids,
            key=lambda other: (nodes[other][0] - target[0]) ** 2 + (nodes[other][1] - target[1]) ** 2,
        )
        distance = math.hypot(nodes[mirror_id][0] - target[0], nodes[mirror_id][1] - target[1])
        if distance > tolerance_um:
            unmatched.append(node_id)
            used.add(node_id)
            continue
        a, b = sorted((node_id, mirror_id))
        if a not in used or b not in used:
            pairs.append((a, b, distance))
            used.add(a)
            used.add(b)
    return pairs, unmatched


def quantile(values: list[float], p: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    return ordered[int(round((len(ordered) - 1) * p))]


def solver_stats(values: dict[int, float],
                 pairs: list[tuple[int, int, float]],
                 gate: float) -> tuple[dict[str, Any], dict[int, float], list[dict[str, Any]]]:
    diffs: list[float] = []
    node_error: dict[int, float] = {}
    rows: list[dict[str, Any]] = []
    for a, b, distance in pairs:
        if a not in values or b not in values:
            continue
        diff = abs(values[a] - values[b])
        rel = diff / max(abs(values[a]), abs(values[b]), 1.0e-300)
        diffs.append(diff)
        node_error[a] = diff
        node_error[b] = diff
        rows.append({
            "left_node": a,
            "right_node": b,
            "mirror_distance_um": distance,
            "left_value_V_per_cm": values[a],
            "right_value_V_per_cm": values[b],
            "abs_diff_V_per_cm": diff,
            "rel_diff": rel,
        })
    rows.sort(key=lambda row: row["abs_diff_V_per_cm"], reverse=True)
    stats = {
        "pairs_compared": len(diffs),
        "mean_abs_diff_V_per_cm": sum(diffs) / len(diffs) if diffs else None,
        "median_abs_diff_V_per_cm": statistics.median(diffs) if diffs else None,
        "p95_abs_diff_V_per_cm": quantile(diffs, 0.95),
        "max_abs_diff_V_per_cm": max(diffs) if diffs else None,
        "classification": "symmetric" if diffs and max(diffs) <= gate else "asymmetric",
        "top_pairs": rows[:12],
    }
    return stats, node_error, rows


def load_font(size: int, bold: bool = False) -> ImageFont.ImageFont:
    candidates = [
        "C:/Windows/Fonts/segoeuib.ttf" if bold else "C:/Windows/Fonts/segoeui.ttf",
        "C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            pass
    return ImageFont.load_default()


def color_for(value: float, lo: float, hi: float) -> tuple[int, int, int]:
    palette = [
        (36, 60, 135),
        (35, 137, 142),
        (113, 190, 109),
        (253, 231, 37),
        (210, 45, 45),
    ]
    z = (math.log10(max(value, 1.0e-12)) - lo) / max(hi - lo, 1.0e-12)
    z = max(0.0, min(1.0, z))
    pos = z * (len(palette) - 1)
    idx = min(int(pos), len(palette) - 2)
    frac = pos - idx
    return tuple(round(palette[idx][i] * (1.0 - frac) + palette[idx + 1][i] * frac) for i in range(3))


def draw_heatmap(path: Path,
                 title: str,
                 nodes: dict[int, tuple[float, float]],
                 triangles: list[tuple[int, list[int]]],
                 node_error: dict[int, float],
                 stats: dict[str, Any],
                 axis_x_um: float) -> None:
    width, height = 1280, 560
    pad_l, pad_r, pad_t, pad_b = 80, 210, 74, 88
    xs = [point[0] for point in nodes.values()]
    ys = [point[1] for point in nodes.values()]
    minx, maxx = min(xs), max(xs)
    miny, maxy = min(ys), max(ys)
    plot_w = width - pad_l - pad_r
    plot_h = height - pad_t - pad_b
    scale = min(plot_w / (maxx - minx), plot_h / (maxy - miny))
    actual_w = scale * (maxx - minx)
    actual_h = scale * (maxy - miny)
    off_x = pad_l + (plot_w - actual_w) / 2.0
    off_y = pad_t + (plot_h - actual_h) / 2.0

    def sx(x: float) -> float:
        return off_x + (x - minx) * scale

    def sy(y: float) -> float:
        return off_y + actual_h - (y - miny) * scale

    values = [
        sum(node_error.get(node, 0.0) for node in ids) / 3.0
        for _, ids in triangles
        if all(node in nodes for node in ids)
    ]
    nonzero = [value for value in values if value > 0.0]
    lo = math.log10(max(min(nonzero), 1.0e-12)) if nonzero else -12.0
    hi = math.log10(max(nonzero)) if nonzero else 0.0
    if hi - lo < 1.0:
        hi = lo + 1.0

    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image, "RGBA")
    title_font = load_font(28, True)
    sub_font = load_font(15)
    axis_font = load_font(13)
    small_font = load_font(12)
    label_font = load_font(13, True)
    draw.text((pad_l, 18), title, font=title_font, fill=(0, 0, 0))
    draw.text(
        (pad_l, 51),
        f"Mirror axis x={axis_x_um:g} um; cell color = average mirrored-node |E(x,y)-E(2*xj-x,y)|, V/cm",
        font=sub_font,
        fill=(45, 45, 45),
    )
    for _, ids in triangles:
        if not all(node in nodes for node in ids):
            continue
        value = sum(node_error.get(node, 0.0) for node in ids) / 3.0
        points = [(sx(nodes[node][0]), sy(nodes[node][1])) for node in ids]
        draw.polygon(points, fill=color_for(value, lo, hi))
    for _, ids in triangles:
        if not all(node in nodes for node in ids):
            continue
        points = [(sx(nodes[node][0]), sy(nodes[node][1])) for node in ids]
        draw.line([points[0], points[1], points[2], points[0]], fill=(255, 255, 255, 60), width=1)
    draw.rectangle((off_x, off_y, off_x + actual_w, off_y + actual_h), outline=(35, 35, 35), width=1)
    axis_px = sx(axis_x_um)
    draw.line((axis_px, off_y, axis_px, off_y + actual_h), fill=(0, 0, 0, 210), width=2)
    draw.text((axis_px + 5, off_y + 8), "mirror axis", font=small_font, fill=(0, 0, 0))

    for i in range(5):
        x = minx + (maxx - minx) * i / 4.0
        px = sx(x)
        draw.line((px, off_y + actual_h, px, off_y + actual_h + 6), fill=(50, 50, 50), width=1)
        text = f"{x:.1f}"
        draw.text((px - draw.textlength(text, font=axis_font) / 2.0, off_y + actual_h + 10), text, font=axis_font, fill=(50, 50, 50))
    for i in range(3):
        y = miny + (maxy - miny) * i / 2.0
        py = sy(y)
        draw.line((off_x - 6, py, off_x, py), fill=(50, 50, 50), width=1)
        text = f"{y:.2f}"
        draw.text((off_x - 10 - draw.textlength(text, font=axis_font), py - 7), text, font=axis_font, fill=(50, 50, 50))
    draw.text((off_x + actual_w / 2.0 - 25, height - 31), "x (um)", font=axis_font, fill=(50, 50, 50))
    draw.text((18, off_y + actual_h / 2.0 - 7), "y (um)", font=axis_font, fill=(50, 50, 50))

    lx, ly = width - pad_r + 52, pad_t + 4
    bar_w, bar_h = 26, 292
    for j in range(bar_h):
        z = 1.0 - j / max(bar_h - 1, 1)
        logv = lo + (hi - lo) * z
        draw.rectangle((lx, ly + j, lx + bar_w, ly + j + 1), fill=color_for(10.0 ** logv, lo, hi))
    draw.rectangle((lx, ly, lx + bar_w, ly + bar_h), outline=(50, 50, 50), width=1)
    draw.text((lx - 10, ly - 24), "log10(err)", font=label_font, fill=(30, 30, 30))
    for frac in [0.0, 0.25, 0.5, 0.75, 1.0]:
        logv = hi - (hi - lo) * frac
        y = ly + bar_h * frac
        draw.line((lx + bar_w, y, lx + bar_w + 6, y), fill=(50, 50, 50), width=1)
        draw.text((lx + bar_w + 11, y - 7), f"1e{logv:.1f}", font=axis_font, fill=(50, 50, 50))
    y0 = ly + bar_h + 42
    lines = [
        f"class {stats['classification']}",
        f"median {stats['median_abs_diff_V_per_cm']:.3g}",
        f"p95    {stats['p95_abs_diff_V_per_cm']:.3g}",
        f"max    {stats['max_abs_diff_V_per_cm']:.3g}",
    ]
    for line in lines:
        draw.text((width - pad_r + 22, y0), line, font=small_font, fill=(40, 40, 40))
        y0 += 16
    draw.text((pad_l, height - 58), "Zero error is drawn with the lowest color; use the JSON/CSV for exact zero values.", font=small_font, fill=(45, 45, 45))
    draw.text((pad_l, height - 40), f"Pairs compared: {stats['pairs_compared']}; values are ElectricField magnitudes in V/cm.", font=small_font, fill=(45, 45, 45))
    image.save(path)


def write_pair_csv(path: Path, sentaurus_rows: list[dict[str, Any]], vela_rows: list[dict[str, Any]]) -> None:
    by_pair: dict[tuple[int, int], dict[str, Any]] = {}
    for row in sentaurus_rows:
        by_pair[(row["left_node"], row["right_node"])] = {
            "left_node": row["left_node"],
            "right_node": row["right_node"],
            "mirror_distance_um": row["mirror_distance_um"],
            "sentaurus_left": row["left_value_V_per_cm"],
            "sentaurus_right": row["right_value_V_per_cm"],
            "sentaurus_abs_diff_V_per_cm": row["abs_diff_V_per_cm"],
        }
    for row in vela_rows:
        entry = by_pair.setdefault((row["left_node"], row["right_node"]), {
            "left_node": row["left_node"],
            "right_node": row["right_node"],
            "mirror_distance_um": row["mirror_distance_um"],
        })
        entry.update({
            "vela_left": row["left_value_V_per_cm"],
            "vela_right": row["right_value_V_per_cm"],
            "vela_abs_diff_V_per_cm": row["abs_diff_V_per_cm"],
        })
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        fieldnames = [
            "left_node", "right_node", "mirror_distance_um",
            "sentaurus_left", "sentaurus_right", "sentaurus_abs_diff_V_per_cm",
            "vela_left", "vela_right", "vela_abs_diff_V_per_cm",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in sorted(by_pair.values(), key=lambda item: max(item.get("sentaurus_abs_diff_V_per_cm", 0.0), item.get("vela_abs_diff_V_per_cm", 0.0)), reverse=True):
            writer.writerow(row)


def markdown_report(report: dict[str, Any]) -> str:
    lines = [
        "# PN2D 0V Electric Field Symmetry",
        "",
        f"Status: {report.get('status')}",
        f"Mirror axis: x={report.get('axis_x_um')} um",
        "",
        "| Solver | Class | Pairs | Median abs diff | P95 abs diff | Max abs diff |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for name, stats in report.get("solvers", {}).items():
        lines.append(
            f"| {name} | {stats.get('classification')} | {stats.get('pairs_compared')} | "
            f"{stats.get('median_abs_diff_V_per_cm')} | {stats.get('p95_abs_diff_V_per_cm')} | "
            f"{stats.get('max_abs_diff_V_per_cm')} |"
        )
    return "\n".join(lines) + "\n"


def default_field_node_csv(reference_root: Path) -> Path:
    return reference_root / "reports" / "0v_electric_field" / "pn2d_0v_electric_field_nodes.csv"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--field-node-csv", type=Path)
    parser.add_argument("--axis-x-um", type=float, default=1.0)
    parser.add_argument("--mirror-tolerance-um", type=float, default=1.0e-6)
    parser.add_argument("--symmetry-abs-gate-v-per-cm", type=float, default=1.0e-6)
    args = parser.parse_args()

    report_path = args.output_dir / REPORT_NAME
    md_path = args.output_dir / MARKDOWN_NAME
    try:
        field_csv = args.field_node_csv or default_field_node_csv(args.reference_root)
        nodes = read_nodes(args.reference_root / "nodes.csv")
        mesh = read_json(args.reference_root / "vela" / "mesh.json")
        sentaurus, vela = read_field_nodes(field_csv)
        pairs, unmatched = find_mirror_pairs(nodes, args.axis_x_um, args.mirror_tolerance_um)
        sentaurus_stats, sentaurus_error, sentaurus_rows = solver_stats(sentaurus, pairs, args.symmetry_abs_gate_v_per_cm)
        vela_stats, vela_error, vela_rows = solver_stats(vela, pairs, args.symmetry_abs_gate_v_per_cm)
        triangles = mesh_triangles(mesh)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        sentaurus_png = args.output_dir / "pn2d_0v_sentaurus_electric_field_symmetry.png"
        vela_png = args.output_dir / "pn2d_0v_vela_electric_field_symmetry.png"
        draw_heatmap(sentaurus_png, "Sentaurus PN2D 0V Electric-Field Mirror Error", nodes, triangles, sentaurus_error, sentaurus_stats, args.axis_x_um)
        draw_heatmap(vela_png, "Vela PN2D 0V Electric-Field Mirror Error", nodes, triangles, vela_error, vela_stats, args.axis_x_um)
        pair_csv = args.output_dir / PAIR_CSV_NAME
        write_pair_csv(pair_csv, sentaurus_rows, vela_rows)
        report = {
            "status": "pass" if pairs else "fail",
            "axis_x_um": args.axis_x_um,
            "mirror_tolerance_um": args.mirror_tolerance_um,
            "symmetry_abs_gate_V_per_cm": args.symmetry_abs_gate_v_per_cm,
            "field_node_csv": str(field_csv),
            "pairing": {
                "nodes": len(nodes),
                "pairs": len(pairs),
                "unmatched_nodes": len(unmatched),
                "unmatched_node_ids": unmatched[:20],
                "max_mirror_distance_um": max((distance for _, _, distance in pairs), default=None),
            },
            "solvers": {
                "sentaurus": sentaurus_stats,
                "vela": vela_stats,
            },
            "outputs": {
                "pair_csv": str(pair_csv),
                "sentaurus_png": str(sentaurus_png),
                "vela_png": str(vela_png),
            },
        }
        write_json(report_path, report)
        md_path.write_text(markdown_report(report))
        return 0 if report["status"] == "pass" else 1
    except Exception as exc:
        report = {"status": "error", "error": str(exc)}
        write_json(report_path, report)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        md_path.write_text(markdown_report(report))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
