"""Render the frozen SLOT daily-report snapshot; --prepare reads local evidence.

This is visualization only: no solver or remote job is started. Coordinates are
swapped for display (horizontal=original y, depth=original x), without rescaling.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib import font_manager
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.patches import Patch
import numpy as np


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def prepare(repo, out):
    w = repo / "build/slot_energy_reference_20260921"
    inputs = []

    def record(path):
        path = path.resolve()
        inputs.append({"path": path.relative_to(repo).as_posix(),
                       "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        return path

    def js(path):
        return json.loads(record(path).read_text(encoding="utf-8"))

    cfg = js(w / "av_boundary_projection/config.json")
    mesh = js(Path(cfg["mesh_file"]))
    nodes = sorted(mesh["nodes"], key=lambda r: r["id"])
    assert [r["id"] for r in nodes] == list(range(len(nodes)))
    xy = np.array([[r["y"], r["x"]] for r in nodes])
    tri = np.array([r["node_ids"] for r in mesh["triangles"]])
    regions = np.array([r["region_id"] for r in mesh["triangles"]])
    arrays = {"xy": xy, "tri": tri, "regions": regions}
    doping = read_csv(record(Path(cfg["node_doping_file"])))
    arrays["net_doping_cm3"] = np.zeros(len(nodes))
    for r in doping:
        arrays["net_doping_cm3"][int(r["node_id"])] = float(r["donors_cm3"]) - float(r["acceptors_cm3"])
    controls = js(w / "av_boundary_controls.json")
    projected = next(r for r in controls if r["case"] == "projection")
    assert projected["qualified"]
    state = read_csv(record(Path(projected["state_file"])))
    arrays["psi_10V"] = np.zeros(len(nodes))
    for r in state:
        arrays["psi_10V"][int(r["node_id"])] = float(r["psi"])
    budget = js(w / "native_av10_local_source/integrated_comparison.json")
    factor = budget["local_term_to_A_per_um"]
    for label, relative in [("hole_before", "native_av10_local_source/local_audit/terms.csv"),
                            ("hole_projected", "av_boundary_projection/local_audit/terms.csv")]:
        arrays[label] = np.zeros(len(nodes))
        for r in read_csv(record(w / relative)):
            arrays[label][int(r["node_id"])] = float(r["impact_hole_source"]) * factor
    geometry = js(repo / "build/slot_element_ialmob_20260921/local_control_r2/native_geometry.json")
    volumes = np.zeros(len(nodes))
    for c in geometry["cells"]:
        for node, volume in zip(c["node_ids"], c["vertex_measure_m2"]):
            volumes[node] += volume
    arrays["hole_native"] = np.zeros(len(nodes))
    for r in read_csv(record(w / "native_av10/export/fields/hImpactIonization_region0.csv")):
        node = int(r["node_id"])
        arrays["hole_native"][node] = float(r["component0"]) * volumes[node] * 1.602176634e-19
    for key in ("hole_before", "hole_projected", "hole_native"):
        assert np.all(np.isfinite(arrays[key])) and np.min(arrays[key]) >= 0
    assert np.isclose(arrays["hole_native"].sum(), 1.6409740694375286e-17, rtol=1e-10, atol=0)
    high = js(w / "qualified_noav_high.json")
    low = js(w / "qualified_projection_av_low.json")
    assert all(r["qualified"] for r in high + low)
    native = js(repo / "build/slot_element_ialmob_20260921/oracle/avalanche_current_control.json")
    metrics = {
        "date": "2026-09-21", "implementation_commit": "076ab59c",
        "nodes": len(nodes), "triangles": len(tri),
        "contacts": [{"name": c["name"], "node_ids": c["node_ids"]} for c in mesh["contacts"]],
        "materials": {str(r["id"]): r["material"] for r in mesh["regions"]},
        "noav": [{k: r[k] for k in ("bias_V", "current_A_per_um", "native_current_A_per_um", "relative_current_error_percent")} for r in high],
        "av_projection": [{k: r[k] for k in ("bias_V", "current_A_per_um", "native_current_A_per_um", "relative_current_error_percent")} for r in low] + [{
            "bias_V": 10, "current_A_per_um": projected["current_A_per_um"],
            "native_current_A_per_um": native["grad"]["10"],
            "relative_current_error_percent": projected["relative_native_error_percent"]}],
        "av10_before_percent": 4.357689270181, "av10_native_sg_A_per_um": native["sg"]["10"],
        "source_sums_A_per_um": {key: float(arrays[key].sum()) for key in ("hole_native", "hole_before", "hole_projected")},
        "test_count": 888, "av_jacobian_relative_difference": 7.466074069896746e-9,
        "scope": "Qualified SG-source projection controls; native AvalDensGradQF parity and BVDS remain incomplete. Higher-voltage work is paused.",
        "coordinate_convention": "horizontal=original y, depth=original x; micrometres; equal geometric scale",
        "source_map_convention": "Nonnegative nodal hole-avalanche source integrals in A/um at the SAME native 10 V state; display interpolates nodal values, not volumetric generation density.",
        "inputs": inputs,
    }
    (out / "data").mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "data/figure_data.npz", **arrays)
    (out / "data/summary.json").write_text(json.dumps(metrics, ensure_ascii=False, indent=2)+"\n", encoding="utf-8", newline="\n")


def render(out):
    data = np.load(out / "data/figure_data.npz", allow_pickle=False)
    meta = json.loads((out / "data/summary.json").read_text(encoding="utf-8"))
    for font in (Path("C:/Windows/Fonts/msyh.ttc"), Path("C:/Windows/Fonts/msyhbd.ttc")):
        if font.exists(): font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({"font.family": ["Microsoft YaHei", "DejaVu Sans"], "font.size": 12,
                         "axes.unicode_minus": False, "axes.spines.top": False,
                         "axes.spines.right": False, "svg.fonttype": "path"})
    xy, tri, regions = (data[k] for k in ("xy", "tri", "regions"))
    triang = mtri.Triangulation(xy[:, 0], xy[:, 1], tri)
    silicon = mtri.Triangulation(xy[:, 0], xy[:, 1], tri, mask=regions != 0)
    colors = {"Si": "#dce7f1", "SiO2": "#f6d68b", "PolySilicon": "#528ec5", "Nitride": "#b1a1cc"}
    pairs = {}
    for t, region in zip(tri, regions):
        for i, j in zip(t, np.roll(t, -1)):
            pairs.setdefault(tuple(sorted((int(i), int(j)))), []).append(int(region))
    boundary = [xy[list(edge)] for edge, rs in pairs.items() if len(rs) == 1 or len(set(rs)) > 1]

    def materials(ax):
        ax.add_collection(PolyCollection(xy[tri], facecolors=[colors[meta["materials"][str(r)]] for r in regions], edgecolors="none"))
        ax.add_collection(LineCollection(boundary, colors="#576775", linewidths=.6))

    def axes(ax, ylim, title):
        ax.set(xlim=(0, 4.69), ylim=ylim, xlabel="横向位置 y（µm）", ylabel="深度 x（µm）")
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(title, loc="left", fontsize=14, fontweight="bold", pad=12)

    def contacts(ax, labels=False):
        for c in meta["contacts"]:
            pts = xy[c["node_ids"]]
            ax.plot(pts[:, 0], pts[:, 1], ".", color="#132d46", markersize=2.5)
            if labels and c["name"] != "substrate":
                pos = pts[np.argmin(np.linalg.norm(pts - pts.mean(axis=0), axis=1))]
                target = {"drain": (.45, -.37), "source": (4.37, -.36), "gate": (4.06, -.44), "SLOT": (2.53, -.46)}[c["name"]]
                ax.annotate(c["name"], pos, xytext=target, ha="center", fontsize=11,
                            arrowprops={"arrowstyle": "-", "color": "#132d46", "lw": .9})

    def canvas(title, subtitle):
        fig = plt.figure(figsize=(16, 9), facecolor="white")
        fig.text(.05, .94, title, fontsize=24, fontweight="bold", color="#18364f")
        fig.text(.05, .895, subtitle, fontsize=12, color="#516778")
        fig.text(.05, .026, "Vela TCAD · SLOT LDMOS · 2026-09-21  |  Release / UMFPACK / 300 K  |  验收门槛未放宽", color="#516778", fontsize=11)
        return fig

    def save(fig, stem):
        (out / "figures").mkdir(exist_ok=True)
        # Dense mesh fields remain high-resolution raster layers in SVG;
        # axes and annotation outlines stay vector without 50 MB Gouraud files.
        for ax in fig.axes:
            for collection in ax.collections:
                collection.set_rasterized(True)
        for ext in ("png", "svg"):
            path = out / "figures" / (stem + "." + ext)
            fig.savefig(path, dpi=220, facecolor="white")
            if ext == "svg":
                path.write_bytes(b"\n".join(line.rstrip() for line in path.read_bytes().splitlines()) + b"\n")
        plt.close(fig)

    fig = canvas("真实器件结构与验证范围", "13,881 个节点 · 27,275 个三角形  |  原始网格坐标换轴显示，几何比例保持一致")
    ax = fig.add_axes((.055, .18, .29, .64)); materials(ax)
    im = ax.tripcolor(silicon, data["net_doping_cm3"], shading="gouraud", cmap="RdBu_r", norm=SymLogNorm(linthresh=1e14, vmin=-1e20, vmax=1e20))
    axes(ax, (8.1, -.52), "A  器件截面 / 净掺杂"); contacts(ax)
    ax.text(2.35, 7.65, "substrate = 0 V", ha="center", fontsize=10)
    cb = fig.colorbar(im, cax=fig.add_axes((.072, .10, .25, .018)), orientation="horizontal", ticks=[-1e20, -1e16, 0, 1e16, 1e20])
    cb.ax.tick_params(labelsize=9); cb.set_label("ND − NA（cm⁻³）；红：n 型，蓝：p 型", fontsize=10)
    ax = fig.add_axes((.405, .56, .55, .26)); materials(ax); axes(ax, (.65, -.55), "B  表面材料与电极"); contacts(ax, True)
    ax.legend(handles=[Patch(facecolor=c, label=n) for n, c in colors.items()], loc="lower center", bbox_to_anchor=(.5, -.49), ncol=4, frameon=False, fontsize=10)
    fig.text(.415, .41, "已验证", color="#168071", fontsize=18, fontweight="bold")
    for y, text in [( .354, "无雪崩：从本地平衡态分段推进到 10 V"),
                    (.306, "自洽雪崩：投影配置的 2 / 5 / 10 V 检查点通过"),
                    (.258, "完整偏导：迁移率—电流—雪崩源保持一致"),
                    (.210, "Release 全量测试：888 / 888 通过")]:
        fig.text(.415, y, text, fontsize=14, color="#18364f")
    fig.text(.415, .126, "未完成：高压连续扫描、外电阻击穿与击穿后支路", fontsize=14, color="#aa572a")
    fig.text(.415, .079, "Vg / Vs / Vsubstrate / VSLOT = 0 V；当前为关断漏极偏压扫描。", fontsize=11, color="#516778")
    save(fig, "01_device_structure")

    fig = canvas("10 V 已验收状态与电流对照", "电势图取自本地自洽雪崩投影检查点；下方仅比较已合格检查点，不代表完整击穿曲线")
    ax = fig.add_axes((.07, .58, .60, .24)); materials(ax)
    im = ax.tripcolor(silicon, data["psi_10V"], shading="gouraud", cmap="viridis", vmin=0, vmax=11)
    axes(ax, (1.05, -.50), "A  表面区域静电势"); contacts(ax)
    fig.colorbar(im, cax=fig.add_axes((.69, .60, .014, .20))).set_label("ψ（V）")
    fig.text(.755, .76, "10 V 自洽雪崩", fontsize=17, color="#18364f", fontweight="bold")
    fig.text(.755, .693, "+4.36% → +0.783%", fontsize=19, color="#168071", fontweight="bold")
    fig.text(.755, .626, "仅增加绝缘界面驱动力投影\n与原生 SG 控制相差 +0.104%", fontsize=12, color="#516778", linespacing=1.8)
    for left, records, title, color in [(.075, meta["noav"], "B  无雪崩：本地续算合格点", "#227f9d"),
                                        (.565, meta["av_projection"], "C  自洽雪崩：投影独立检查点", "#168071")]:
        ax = fig.add_axes((left, .145, .365, .27))
        x = np.arange(len(records)); vals = [r["relative_current_error_percent"] for r in records]
        ax.axhline(0, color="#91a2b0", lw=1)
        ax.bar(x, vals, color=color, width=.48)
        ax.set_xticks(x, [f'{r["bias_V"]} V' for r in records]); ax.set_ylabel("相对原始 Sentaurus 的电流差（%）")
        ax.set_title(title, loc="left", fontweight="bold", fontsize=14, pad=15)
        ax.set_ylim((-.07, .07) if left < .5 else (-.18, .94))
        ax.grid(axis="y", alpha=.15); ax.set_axisbelow(True)
        for xx, val in zip(x, vals): ax.annotate(f"{val:+.4f}%", (xx, val), xytext=(0, 6 if val >= 0 else -15), textcoords="offset points", ha="center", fontsize=12)
    fig.text(.08, .077, "误差定义：(IVela / ISentaurus − 1) × 100%；各图纵轴范围不同。SG 源项控制尚未严格还原 AvalDensGradQF。", fontsize=11, color="#516778")
    save(fig, "02_potential_and_progress")

    fig = canvas("雪崩偏差定位：绝缘界面驱动力", "固定同一个原生 10 V 状态，仅比较源项算子；三图使用相同色标，排除自洽场变化影响")
    for left, key, title in [(.055, "hole_native", "原生 Sentaurus"), (.377, "hole_before", "Vela · 全梯度"), (.699, "hole_projected", "Vela · 界面切向投影")]:
        ax = fig.add_axes((left, .48, .263, .32))
        materials(ax)
        im = ax.tripcolor(silicon, np.maximum(data[key], 1e-25), shading="gouraud", cmap="magma", norm=LogNorm(1e-23, 3e-18))
        ax.add_collection(LineCollection(boundary, colors="#8392a1", linewidths=.55))
        ax.set(xlim=(3.45, 4.04), ylim=(.30, -.035), xlabel="横向 y（µm）", ylabel="深度 x（µm）")
        ax.set_aspect("equal"); ax.set_title(title, fontsize=15, fontweight="bold", pad=14)
        p = xy[10493]; ax.plot(*p, marker="o", ms=7, mfc="none", mec="#38d9cb", mew=1.4)
    cb = fig.colorbar(im, cax=fig.add_axes((.19, .405, .62, .018)), orientation="horizontal")
    cb.set_label("空穴雪崩节点源积分（A/µm，对数色标）；低于 10⁻²³ 按色标下限显示", fontsize=11)
    ax = fig.add_axes((.08, .175, .39, .155))
    sums = [meta["source_sums_A_per_um"][key] / 1e-17 for key in ("hole_native", "hole_before", "hole_projected")]
    ax.barh([2, 1, 0], sums, color=["#536f88", "#c77950", "#168071"], height=.56)
    ax.set_yticks([2, 1, 0], ["原生", "全梯度", "投影"]); ax.set_xlim(0, 11.5)
    ax.set_xlabel("全 Si 空穴雪崩源积分（×10⁻¹⁷ A/µm）", fontsize=11)
    for y, val in zip([2, 1, 0], sums): ax.text(val+.16, y, f"{val:.3f}", va="center", fontsize=12)
    ax.grid(axis="x", alpha=.15); ax.set_axisbelow(True)
    fig.text(.57, .30, "圆圈：界面节点 10493", fontsize=15, fontweight="bold", color="#18364f")
    fig.text(.57, .247, "电离系数：全梯度约为原生的 33.35 倍", fontsize=14, color="#18364f")
    fig.text(.57, .196, "切向投影后，该节点相对差约 5.65×10⁻⁶", fontsize=14, color="#168071")
    fig.text(.57, .138, "仍有剩余源项差异；不能据此宣称全域或 BVDS 对齐。", fontsize=11, color="#aa572a")
    fig.text(.08, .077, "图中插值用于显示节点积分量，不是体产生率；柱图按全 Si 节点求和，三幅局部图不覆盖全部积分区域。", fontsize=10, color="#516778")
    save(fig, "03_avalanche_boundary")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, default=Path("docs/reports/slot_ldmos_daily_2026-09-21"))
    parser.add_argument("--prepare", action="store_true", help="extract a frozen plotting snapshot from local evidence")
    args = parser.parse_args()
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=True)
    if args.prepare: prepare(args.repo_root.resolve(), out)
    render(out)
    print(out)


if __name__ == "__main__":
    main()
