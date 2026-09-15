"""Export the 2026-09-13 LDMOS daily-report figures from frozen evidence.

Run from the worktree root with UCRT64 Python. No solver execution is required.
Generated PNG/PDF files and chart-ready CSVs are explicitly requested artifacts.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.collections import LineCollection
from matplotlib.colors import TwoSlopeNorm
from matplotlib.tri import Triangulation
import numpy as np

from analyze_templates_ldmos_stage4_d5 import read_curve


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, default=Path('reference_staging/templates_ldmos_native_poisson_20260913'))
    parser.add_argument('--baseline', type=Path, default=Path('reference_staging/templates_ldmos_generation_alignment_20260913'))
    parser.add_argument('--output', type=Path, default=Path('docs/validation/figures/ldmos_2026-09-13'))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    used = {}

    def source(path, expected=None):
        path = Path(path)
        actual = sha(path)
        if expected is not None:
            assert actual == expected, f'Changed evidence: {path}'
        used[str(path.resolve())] = actual
        return read(path)

    final_path = args.evidence / 'p3_full_final.json'
    final = source(final_path)
    old = source(args.baseline / 'b4_r3_full_final.json')
    assert final['status'] == old['status'] == 'pass'
    assert len(final['thermal']) == len(old['thermal']) == 62
    assert all(t['qualification_pass'] for t in final['thermal'] + old['thermal'])
    hashes = {str(Path(p).resolve()): h for p, h in final['sources_sha256'].items()}

    def checked(path):
        path = Path(path)
        return source(path, hashes[str(path.resolve())])

    work = source(args.evidence / 'full_work_comparison.json')
    benchmark = source(args.evidence / 'benchmark_serial_summary.json')
    local = source(args.evidence / 'field_decomposition_full_summary.json')
    old_local = source(args.baseline / 'field_decomposition_r3_summary.json')
    assert len(benchmark['runs']) == 8
    assert all(r['exact_equal'] and r['gate']['pass_gate'] for r in benchmark['runs'])
    native_root = Path('reference_staging/templates_ldmos_d0_electrothermal_20260912/native_fields_r1')
    curves, chart_rows = {}, []
    for g in (4, 8):
        ledger = checked(args.evidence / f'p3_full_vg{g}/ledger.json')
        assert ledger['status'] == 'complete' and len(ledger['exact_points']) == 31
        assert sha(args.evidence / f'p3_full_vg{g}/ledger.json') == work['versions']['p3'][str(g)]['ledger_sha256']
        path = native_root / 'normalized' / f'IdVd_Vg{1 if g == 4 else 2}_n4_des_drain_curve.csv'
        assert sha(path) == hashes[str(path.resolve())]
        used[str(path.resolve())] = sha(path)
        native = read_curve(path)
        thermal = sorted((t for t in final['thermal'] if t['gate'] == g), key=lambda t: t['bias_V'])
        points = ledger['exact_points']
        bias = np.array([p['bias_V'] for p in points])
        assert len(native) == 31
        assert np.allclose(bias, [v for v, _ in native], atol=1e-9, rtol=0)
        assert np.allclose(bias, [t['bias_V'] for t in thermal], atol=1e-9, rtol=0)
        current = []
        for p in points:
            state = checked(p['result'])
            current.append(next(c['total_outflow_A_per_m'] for c in state['contacts'] if c['contact'] == 'drain'))
        # Only the last state is retained; the other 60 large result files are released.
        reference = np.array([i for _, i in native]) * 1e6
        current = np.array(current)
        error = np.full(31, np.nan)
        error[1:] = 100 * np.abs(current[1:] - reference[1:]) / np.abs(reference[1:])
        assert np.isclose(error[1:].max(), final['metrics'][f'Vg{g}']['relative_error_percent']['max'], atol=1e-9)
        peak = np.array([t['candidate_peak_K'] for t in thermal])
        native_peak = np.array([t['reference_peak_K'] for t in thermal])
        curves[g] = dict(bias=bias, current=current, reference=reference, error=error,
                         peak=peak, native_peak=native_peak, temperature=state['temperature_K'])
        for k, vd in enumerate(bias):
            chart_rows.append(dict(gate_V=g, drain_V=vd, vela_current_A_per_m=current[k],
                native_current_A_per_m=reference[k], current_error_percent=None if k == 0 else error[k],
                vela_peak_K=peak[k], native_peak_K=native_peak[k]))
        print(f'Checked gate {g}: all 31 exact states and reference curves', flush=True)

    font = Path('C:/Windows/Fonts/msyh.ttc')
    assert font.exists(), 'Chinese font required'
    font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({'font.family': font_manager.FontProperties(fname=str(font)).get_name(),
        'font.size': 11, 'axes.titlesize': 13, 'axes.labelsize': 11, 'axes.unicode_minus': False,
        'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42,
        'savefig.facecolor': 'white', 'axes.edgecolor': '#a5adb5', 'text.color': '#243342'})
    blue, orange, grey = '#2466a4', '#c47722', '#8996a3'
    colors = {4: blue, 8: orange}
    exports = []

    def finish(fig, stem, note, left=.075, bottom=.19, top=.82, wspace=.32):
        fig.text(left, .035, note, fontsize=9.2, color='#526170', va='bottom', linespacing=1.65)
        fig.subplots_adjust(left=left, right=.97, top=top, bottom=bottom, wspace=wspace)
        for ax in fig.axes:
            ax.set_axisbelow(True)
        for ext in ('png', 'pdf'):
            p = args.output / f'{stem}.{ext}'
            fig.savefig(p, dpi=240)
            exports.append(p)
        plt.close(fig)
        print(f'Exported {stem}', flush=True)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.8))
    fig.suptitle('D0 双栅压完整曲线：电流与自热温升', x=.075, y=.97, ha='left', fontsize=19)
    fig.text(.075, .89, 'Vd = 0–40 V · 每栅压 31 个精确点 · 原电学与批准热学门限全部通过', fontsize=11)
    for g, c in curves.items():
        for ax, candidate, reference in ((axes[0], c['current'], c['reference']),
                                          (axes[1], c['peak']-300, c['native_peak']-300)):
            ax.plot(c['bias'], reference, color=colors[g], lw=1.9, label=f'Sentaurus · Vg={g} V')
            ax.plot(c['bias'], candidate, 'o', ms=4, mfc='white', mew=1.15, color=colors[g], label=f'Vela P3 · Vg={g} V')
            ax.set(xlim=(0, 40), xlabel='漏极电压 Vd (V)')
            ax.grid(alpha=.17)
    axes[0].set(title='单位宽度漏极电流', ylabel='Id (A/m)')
    axes[1].set(title='全域峰值温升', ylabel='Tmax − 300 K (K)')
    current_limit = max(max(c['current'].max(), c['reference'].max()) for c in curves.values())*1.10
    rise_limit = max(max(c['peak'].max(), c['native_peak'].max())-300 for c in curves.values())*1.10
    axes[0].set_ylim(0, current_limit)
    axes[1].set_ylim(0, rise_limit)
    assert all(c['current'].max() < current_limit and c['reference'].max() < current_limit for c in curves.values())
    assert all(c['peak'].max()-300 < rise_limit and c['native_peak'].max()-300 < rise_limit for c in curves.values())
    axes[0].legend(frameon=False, fontsize=9, loc='lower right')
    finish(fig, '01_d0_full_curves', '实线：Sentaurus；空心点：Vela P3。曲线高度重合，误差量级另见图 02。\n同一 10241 节点网格；有限空穴接触，关闭 Auger 生成；独立 D0 电热入口。')

    error_rows = []
    for label, data, fields in (('B4 R3', old, old_local), ('N1 / P3', final, local)):
        for g in (4, 8):
            t = [x for x in data['thermal'] if x['gate'] == g]
            error_rows.append(dict(version=label, gate_V=g,
                current_max_error_percent=data['metrics'][f'Vg{g}']['relative_error_percent']['max'],
                peak_max_error_K=max(x['gates']['peak_temperature_rise']['error_K'] for x in t),
                conduction_band_max_error_eV=abs(fields['gates'][str(g)]['Ec_max_node']['Ec_difference_eV'])))
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.8))
    fig.suptitle('原生 Poisson 几何对齐前后的误差', x=.07, y=.97, ha='left', fontsize=19)
    fig.text(.07, .89, 'B4 R3：原几何；N1 / P3：原生边系数 + 硅电荷体积，性能优化保持接受状态一致', fontsize=10.5)
    for ax, key, title, unit in zip(axes,
            ('current_max_error_percent', 'peak_max_error_K', 'conduction_band_max_error_eV'),
            ('完整曲线最大电流误差', '完整曲线最大峰温误差', '40 V 最大导带差'), ('误差 (%)', '误差 (K)', '绝对差 (eV)')):
        for j, label in enumerate(('B4 R3', 'N1 / P3')):
            values = [r[key] for r in error_rows if r['version'] == label]
            bars = ax.bar(np.arange(2)+(j-.5)*.32, values, width=.30, color=grey if j == 0 else blue,
                          hatch='//' if j == 0 else None, label=label)
            ax.bar_label(bars, labels=[f'{v:.5f}' if v < .1 else f'{v:.3f}' for v in values], padding=5, fontsize=10)
        ax.set(xticks=[0, 1], xticklabels=['Vg=4 V', 'Vg=8 V'], title=title, ylabel=unit)
        ax.set_ylim(0, max(r[key] for r in error_rows)*1.30)
        ax.grid(axis='y', alpha=.17)
    axes[0].legend(frameon=False, fontsize=9)
    finish(fig, '02_geometry_alignment_errors', '电流相对误差统计非零漏压 30 点；零压另验严格零电流。峰温误差覆盖每栅压 31 点。\n导带差统计 40 V 自洽状态的硅节点，最大偏差节点可随版本改变；局部余差仍存在。', left=.07, wspace=.35)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.8))
    fig.suptitle('同输入高压点：串行性能对照', x=.075, y=.97, ha='left', fontsize=19)
    fig.text(.075, .89, 'Vg=8 V / Vd=40 V · 每次均 7 次 Newton 更新 · 全部状态、残差与轨迹逐项一致', fontsize=11)
    names = ['N1\n新几何基线', 'P1\n+物性准备', 'P2\n+结构复用', 'P3\n+屏蔽根复用']
    for ax, key, title in zip(axes, ('wall_seconds', 'cpu_seconds'), ('子进程墙钟', '子进程 CPU 时间')):
        values = [benchmark['medians'][v][key] for v in ('n1', 'p1', 'p2', 'p3')]
        bars = ax.bar(range(4), values, width=.56, color=[grey, grey, grey, blue])
        for j, v in enumerate(('n1', 'p1', 'p2', 'p3')):
            runs = [r[key] for r in benchmark['runs'] if r['label'] == v]
            ax.scatter([j-.10, j+.10], runs, marker='D', s=23, facecolors='white', edgecolors='#273746', zorder=4)
            ax.annotate(f'{values[j]:.2f}', (j, max(runs)), xytext=(0, 12), textcoords='offset points', ha='center', fontsize=11)
        ax.set(xticks=range(4), xticklabels=names, ylabel='时间 (s)', title=title)
        ax.set_ylim(0, max(r[key] for r in benchmark['runs'])*1.23)
        ax.grid(axis='y', alpha=.17)
    finish(fig, '03_serial_performance', '柱与数值：每版两次测量的中位数；空心菱形：原始测量（并非置信区间）。顺序 N1→P1→P2→P3→P3→P2→P1→N1。\nP3 相对 N1：墙钟减少 62.4%，CPU 减少 66.0%。P1/P2 单独未显示明确短点加速。', bottom=.23)

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.8))
    fig.suptitle('完整曲线的求解工作量变化', x=.07, y=.97, ha='left', fontsize=19)
    fig.text(.07, .89, '同一新几何 · 原门限不变 · N1 与 P3 的全部接受状态一致', fontsize=11)
    metrics = [('newton_updates', 'Newton 更新次数'), ('rejected_updates', '拒绝尝试内的更新次数'), ('symbolic_analyses', 'SparseLU 符号分析次数')]
    count_rows = []
    for ax, (key, title) in zip(axes, metrics):
        for j, version in enumerate(('n1', 'p3')):
            values = [work['versions'][version][str(g)][key] for g in (4, 8)]
            bars = ax.bar(np.arange(2)+(j-.5)*.32, values, width=.30, color=grey if j == 0 else blue,
                          hatch='//' if j == 0 else None, label=version.upper())
            ax.bar_label(bars, padding=4, fontsize=11)
            count_rows.extend(dict(version=version, gate_V=g, metric=key, count=v) for g, v in zip((4, 8), values))
        ax.set(xticks=[0, 1], xticklabels=['Vg=4 V', 'Vg=8 V'], title=title, ylabel='次数')
        ax.set_ylim(0, max(work['versions'][v][str(g)][key] for v in ('n1', 'p3') for g in (4, 8))*1.25)
        ax.grid(axis='y', alpha=.17)
    axes[0].legend(frameon=False, fontsize=9)
    finish(fig, '04_solver_work', '停滞尝试在第 12 次更新后拒绝并减步，替代耗尽 60 次预算；没有把停滞状态判为收敛。\n只复用完全相同稀疏结构的符号分析；数值分解仍逐次执行，次数等于 Newton 更新。', left=.07)

    cfg = source(args.evidence / 'p3_vg8_zero.json')
    mesh_path = Path(cfg['mesh_file'])
    manifest = checked(native_root / 'manifest.json')
    mesh_hashes = {str(Path(p).resolve()): h for p, h in manifest['sources_sha256'].items()}
    mesh = source(mesh_path, mesh_hashes[str(mesh_path.resolve())])
    nodes = sorted(mesh['nodes'], key=lambda n: n['id'])
    assert [n['id'] for n in nodes] == list(range(10241))
    xy = np.array([[n['x'], n['y']] for n in nodes])
    triangles = np.array([t['node_ids'] for t in mesh['triangles']])
    triangulation = Triangulation(xy[:, 0], xy[:, 1], triangles)
    field = next(f for f in manifest['fields'] if f['gate_V'] == 8 and f['point_index'] == 30)
    native = checked(field['temperature_file'])
    assert native['node_id'] == list(range(10241))
    reference = np.array(native['temperature_K'])
    candidate = np.array(curves[8]['temperature'])
    difference = candidate-reference
    limit = max(.01, float(np.ceil(np.max(np.abs(difference))*100)/100))
    edges = {}
    for t in mesh['triangles']:
        ids = t['node_ids']
        for a, b in zip(ids, ids[1:]+ids[:1]):
            edges.setdefault(tuple(sorted((a, b))), []).append(t['region_id'])
    outlines = [xy[list(e)] for e, regions in edges.items() if len(regions) == 1 or len(set(regions)) > 1]
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.8))
    fig.suptitle('高压自热温度场：Vg=8 V，Vd=40 V', x=.055, y=.97, ha='left', fontsize=19)
    fig.text(.055, .88, '同一 10241 节点网格及原始坐标；前两图共用温标，差值图使用独立温差标尺', fontsize=10.5)
    for ax, values, title in zip(axes, (reference, candidate, difference), ('Sentaurus', 'Vela P3', '差值：Vela − Sentaurus')):
        norm = TwoSlopeNorm(vcenter=0, vmin=-limit, vmax=limit) if ax is axes[2] else plt.Normalize(300, 520)
        pc = ax.tripcolor(triangulation, values, shading='gouraud', cmap='RdBu_r' if ax is axes[2] else 'inferno', norm=norm, rasterized=True)
        ax.add_collection(LineCollection(outlines, colors='#a6adb4', linewidths=.45))
        ax.set(title=title, xlabel='网格 x (μm)', xlim=(-10.32, 10), ylim=(11, 0), aspect='equal')
        ax.set_xticks([-10, 0, 10])
        fig.colorbar(pc, ax=ax, orientation='horizontal', pad=.23, fraction=.08, label='温差 (K)' if ax is axes[2] else '晶格温度 (K)')
    axes[0].set_ylabel('网格 y (μm)')
    assert min(reference.min(), candidate.min()) >= 300-1e-9
    assert max(reference.max(), candidate.max()) <= 520
    assert np.abs(difference).max() <= limit
    t = next(t for t in final['thermal'] if t['gate'] == 8 and abs(t['bias_V']-40) < 1e-9)
    rms = t['gates']['temperature_rise_field']['rms_error_K']
    peak_error = abs(candidate.max()-reference.max())
    finish(fig, '05_temperature_field', f'峰温：Sentaurus {reference.max():.3f} K；Vela {candidate.max():.3f} K；峰温误差 {peak_error:.5f} K。\n体积加权温度场 RMS 误差 {rms:.5f} K；灰线为材料界面/外边界，色标覆盖全部节点极值。', left=.055, bottom=.30, top=.80, wspace=.22)

    def csv_export(stem, rows):
        path = args.output / f'{stem}.csv'
        with path.open('w', encoding='utf-8-sig', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        exports.append(path)

    csv_export('curves', chart_rows)
    csv_export('alignment_errors', error_rows)
    csv_export('serial_runs', [dict(index=r['index'], version=r['label'], wall_seconds=r['wall_seconds'], cpu_seconds=r['cpu_seconds']) for r in benchmark['runs']])
    csv_export('solver_work', count_rows)
    used[str(Path(__file__).resolve())] = sha(__file__)
    provenance = dict(scope='2026-09-13 daily report, frozen D0 P3 and B4 R3; existing simulation data only',
        input_sha256=used, input_commit='28350b6',
        field_max_absolute_difference_K=float(np.abs(difference).max()),
        field_peak_error_K=float(peak_error), field_volume_weighted_rms_error_K=rms,
        validation=dict(exact_points=62, all_original_gates_pass=True, serial_runs=8,
            serial_results_exact=True, current_errors_recomputed=True,
            source_hashes_verified=True, node_ids_matched=True, color_scales_unclipped=True),
        outputs_sha256={p.name: sha(p) for p in exports})
    (args.output / 'figure_provenance.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in provenance.items() if k not in ('input_sha256', 'outputs_sha256')}, indent=2))


if __name__ == '__main__':
    main()
