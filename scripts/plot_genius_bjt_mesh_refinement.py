#!/usr/bin/env python3
"""Plot the accepted and locally refined Genius BJT meshes side by side."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
REFINED_ROOT = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "genius_bjt_sentaurus2022"
    / "mesh_sensitivity"
    / "local_refined"
)
DEFAULT_OUTPUT = FIXTURE / "figures" / "mesh_sensitivity"

WINDOWS = {
    "base_collector": {
        "label": "基区—基集结加密区",
        "bounds": (0.75, 0.42, 5.25, 0.88),
        "color": "#E68613",
    },
    "emitter_base_left": {
        "label": "左侧发射结加密区",
        "bounds": (2.52, 0.00, 2.98, 0.55),
        "color": "#8C5AAE",
    },
    "emitter_base_right": {
        "label": "右侧发射结加密区",
        "bounds": (4.02, 0.00, 4.48, 0.55),
        "color": "#8C5AAE",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_mesh(path: Path) -> tuple[list[tuple[float, float]], list[tuple[int, int, int]]]:
    mesh = json.loads(path.read_text(encoding="utf-8"))
    nodes = [(float(node["x"]), float(node["y"])) for node in mesh["nodes"]]
    triangles = [tuple(int(value) for value in cell["node_ids"]) for cell in mesh["triangles"]]
    return nodes, triangles


def unique_segments(
    nodes: list[tuple[float, float]], triangles: list[tuple[int, int, int]]
) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    edges: set[tuple[int, int]] = set()
    for a, b, c in triangles:
        for first, second in ((a, b), (b, c), (c, a)):
            edges.add((first, second) if first < second else (second, first))
    return [(nodes[first], nodes[second]) for first, second in sorted(edges)]


def within(point: tuple[float, float], bounds: tuple[float, float, float, float]) -> bool:
    x, y = point
    xmin, ymin, xmax, ymax = bounds
    return xmin <= x <= xmax and ymin <= y <= ymax


def window_counts(
    nodes: list[tuple[float, float]], triangles: list[tuple[int, int, int]]
) -> dict[str, dict[str, int]]:
    result = {}
    for name, spec in WINDOWS.items():
        bounds = spec["bounds"]
        node_count = sum(within(point, bounds) for point in nodes)
        triangle_count = 0
        for triangle in triangles:
            centroid = (
                sum(nodes[node][0] for node in triangle) / 3.0,
                sum(nodes[node][1] for node in triangle) / 3.0,
            )
            triangle_count += within(centroid, bounds)
        result[name] = {"nodes": node_count, "triangles_by_centroid": triangle_count}
    return result


def add_mesh(
    axis: plt.Axes,
    segments: list[tuple[tuple[float, float], tuple[float, float]]],
    *,
    color: str,
    bounds: tuple[float, float, float, float],
    linewidth: float,
    alpha: float,
) -> None:
    collection = LineCollection(segments, colors=color, linewidths=linewidth, alpha=alpha)
    axis.add_collection(collection)
    xmin, ymin, xmax, ymax = bounds
    axis.set_xlim(xmin, xmax)
    axis.set_ylim(ymax, ymin)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("横向位置 x（μm）")
    axis.set_ylabel("深度 y（μm）")
    axis.grid(False)
    for spine in axis.spines.values():
        spine.set_color("#606B78")
        spine.set_linewidth(0.8)


def add_windows(axis: plt.Axes, *, labels: bool) -> None:
    for name, spec in WINDOWS.items():
        xmin, ymin, xmax, ymax = spec["bounds"]
        axis.add_patch(
            Rectangle(
                (xmin, ymin),
                xmax - xmin,
                ymax - ymin,
                fill=False,
                edgecolor=spec["color"],
                linewidth=1.8,
                linestyle="--",
                zorder=5,
            )
        )
        if labels:
            if name == "base_collector":
                axis.text(
                    0.88,
                    0.98,
                    spec["label"],
                    color=spec["color"],
                    fontsize=9,
                    weight="bold",
                    ha="left",
                    va="bottom",
                )
            elif name == "emitter_base_left":
                axis.text(
                    2.50,
                    0.05,
                    "发射结左侧",
                    color=spec["color"],
                    fontsize=8,
                    weight="bold",
                    ha="right",
                    va="top",
                )
            else:
                axis.text(
                    4.50,
                    0.05,
                    "发射结右侧",
                    color=spec["color"],
                    fontsize=8,
                    weight="bold",
                    ha="left",
                    va="top",
                )


def configure_fonts() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Microsoft YaHei", "SimHei", "Arial Unicode MS", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "text.color": "#28313B",
            "axes.labelcolor": "#384452",
            "xtick.color": "#566473",
            "ytick.color": "#566473",
        }
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-mesh", type=Path, default=FIXTURE / "vela" / "input" / "mesh.json")
    parser.add_argument("--refined-mesh", type=Path, default=REFINED_ROOT / "vela" / "input" / "mesh.json")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    configure_fonts()
    baseline_nodes, baseline_triangles = load_mesh(args.baseline_mesh)
    refined_nodes, refined_triangles = load_mesh(args.refined_mesh)
    baseline_segments = unique_segments(baseline_nodes, baseline_triangles)
    refined_segments = unique_segments(refined_nodes, refined_triangles)
    baseline_counts = window_counts(baseline_nodes, baseline_triangles)
    refined_counts = window_counts(refined_nodes, refined_triangles)

    fig, axes = plt.subplots(
        3,
        2,
        figsize=(16, 9.2),
        constrained_layout=False,
        gridspec_kw={"height_ratios": [2.2, 0.9, 1.5]},
    )
    fig.patch.set_facecolor("white")
    fig.suptitle("Genius NPN BJT器件网格局部加密前后对比", fontsize=19, weight="bold", y=0.995)
    fig.text(
        0.5,
        0.952,
        "几何、接触和掺杂保持不变；彩色虚线框为SDE新增局部加密窗口",
        ha="center",
        fontsize=10.5,
        color="#566473",
    )

    panels = (
        (axes[0, 0], baseline_segments, "#677483", (0.0, 0.0, 6.0, 2.0), 0.19, 0.66),
        (axes[0, 1], refined_segments, "#1F5B8F", (0.0, 0.0, 6.0, 2.0), 0.15, 0.62),
        (axes[1, 0], baseline_segments, "#677483", (0.65, 0.35, 5.35, 0.95), 0.32, 0.72),
        (axes[1, 1], refined_segments, "#1F5B8F", (0.65, 0.35, 5.35, 0.95), 0.24, 0.67),
        (axes[2, 0], baseline_segments, "#677483", (2.35, 0.0, 4.65, 0.62), 0.34, 0.74),
        (axes[2, 1], refined_segments, "#1F5B8F", (2.35, 0.0, 4.65, 0.62), 0.23, 0.68),
    )
    for axis, segments, color, bounds, linewidth, alpha in panels:
        add_mesh(axis, segments, color=color, bounds=bounds, linewidth=linewidth, alpha=alpha)

    axes[0, 0].set_title(f"(a) 加密前：{len(baseline_nodes)}节点 / {len(baseline_triangles)}三角形")
    axes[0, 1].set_title(f"(b) 加密后：{len(refined_nodes)}节点 / {len(refined_triangles)}三角形")
    axes[1, 0].set_title(
        f"(c) 加密前：基区—基集结局部（{baseline_counts['base_collector']['nodes']}节点）"
    )
    axes[1, 1].set_title(
        f"(d) 加密后：基区—基集结局部（{refined_counts['base_collector']['nodes']}节点）"
    )
    baseline_emitter = sum(
        baseline_counts[name]["nodes"] for name in ("emitter_base_left", "emitter_base_right")
    )
    refined_emitter = sum(
        refined_counts[name]["nodes"] for name in ("emitter_base_left", "emitter_base_right")
    )
    axes[2, 0].set_title(f"(e) 加密前：发射结两侧（窗口节点合计{baseline_emitter}）")
    axes[2, 1].set_title(f"(f) 加密后：发射结两侧（窗口节点合计{refined_emitter}）")

    add_windows(axes[0, 0], labels=False)
    add_windows(axes[0, 1], labels=True)
    for axis in axes[1, :]:
        xmin, ymin, xmax, ymax = WINDOWS["base_collector"]["bounds"]
        axis.add_patch(
            Rectangle(
                (xmin, ymin), xmax - xmin, ymax - ymin,
                fill=False, edgecolor=WINDOWS["base_collector"]["color"], linewidth=2.0, linestyle="--"
            )
        )
    for axis in axes[2, :]:
        for name in ("emitter_base_left", "emitter_base_right"):
            xmin, ymin, xmax, ymax = WINDOWS[name]["bounds"]
            axis.add_patch(
                Rectangle(
                    (xmin, ymin), xmax - xmin, ymax - ymin,
                    fill=False, edgecolor=WINDOWS[name]["color"], linewidth=2.0, linestyle="--"
                )
            )

    legend = [
        Line2D([0], [0], color="#677483", lw=1.4, label="加密前网格"),
        Line2D([0], [0], color="#1F5B8F", lw=1.4, label="加密后网格"),
        Line2D([0], [0], color="#E68613", lw=2.0, ls="--", label="基区—基集结加密窗口"),
        Line2D([0], [0], color="#8C5AAE", lw=2.0, ls="--", label="发射结左右加密窗口"),
    ]
    fig.legend(handles=legend, loc="lower center", ncol=4, frameon=False, fontsize=10, bbox_to_anchor=(0.5, 0.012))
    fig.subplots_adjust(left=0.065, right=0.982, top=0.885, bottom=0.082, hspace=0.34, wspace=0.15)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "mesh_refinement_comparison.png"
    fig.savefig(output, dpi=240, facecolor="white", bbox_inches="tight")
    plt.close(fig)

    manifest = {
        "schema": "vela.genius_bjt_mesh_refinement_figure.v1",
        "baseline": {
            "nodes": len(baseline_nodes),
            "triangles": len(baseline_triangles),
            "window_counts": baseline_counts,
            "mesh_sha256": sha256(args.baseline_mesh),
        },
        "refined": {
            "nodes": len(refined_nodes),
            "triangles": len(refined_triangles),
            "window_counts": refined_counts,
            "mesh_sha256": sha256(args.refined_mesh),
        },
        "windows": WINDOWS,
        "output": str(output.resolve()),
        "output_sha256": sha256(output),
    }
    manifest_path = args.output_dir / "mesh_refinement_figure_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
