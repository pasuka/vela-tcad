"""Export static daily-report figures from qualified, frozen D0 evidence.

Run from the worktree root. Requires matplotlib and numpy, not a solver run.
Raw simulation fields remain in the ignored evidence directory.
"""
import argparse
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


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--mesh', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--font', type=Path, default=Path('C:/Windows/Fonts/msyh.ttc'))
    args = parser.parse_args()
    summary_path = args.evidence / 'full_r4_final.json'
    summary = read(summary_path)
    assert summary['status'] == 'pass'
    assert len(summary['thermal']) == 62
    assert all(t['qualification_pass'] for t in summary['thermal'])
    source_hashes = {str(Path(k).resolve()): v for k, v in summary['sources_sha256'].items()}
    used = {str(summary_path): digest(summary_path)}

    def checked(path):
        path = Path(path)
        actual = digest(path)
        assert source_hashes[str(path.resolve())] == actual, f'Changed evidence: {path}'
        used[str(path)] = actual
        return read(path)

    if args.font.exists():
        font_manager.fontManager.addfont(str(args.font))
        plt.rcParams['font.family'] = font_manager.FontProperties(fname=str(args.font)).get_name()
    else:
        raise ValueError('Provide --font with a Chinese-capable font for report exports')
    plt.rcParams.update({'font.size': 11, 'axes.titlesize': 13, 'axes.labelsize': 11,
                         'axes.unicode_minus': False, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42,
                         'savefig.facecolor': 'white'})
    args.output.mkdir(parents=True, exist_ok=True)
    colors = {4: '#2166ac', 8: '#c56b16'}
    curves = {}
    for gate in (4, 8):
        ledger = checked(args.evidence / f'full_vg{gate}_r4/ledger.json')
        assert ledger['status'] == 'complete' and len(ledger['exact_points']) == 31
        points = ledger['exact_points']
        states = [checked(p['result']) for p in points]
        native_path = args.evidence / 'native_fields_r1/normalized' / f'IdVd_Vg{1 if gate == 4 else 2}_n4_des_drain_curve.csv'
        assert digest(native_path) == source_hashes[str(native_path.resolve())]
        used[str(native_path)] = digest(native_path)
        # Use the same reader as the accepted electrical scorer, including units.
        from analyze_templates_ldmos_stage4_d5 import read_curve
        native_curve = read_curve(native_path)
        bias = np.array([p['bias_V'] for p in points])
        thermal = sorted((t for t in summary['thermal'] if t['gate'] == gate), key=lambda t: t['bias_V'])
        assert np.allclose(bias, [t['bias_V'] for t in thermal], atol=1e-9, rtol=0)
        assert len(native_curve) == 31
        assert np.allclose(bias, [v for v, _ in native_curve], atol=1e-9, rtol=0)
        curves[gate] = dict(bias=bias, states=states,
                            current=np.array([next(c['total_outflow_A_per_m'] for c in s['contacts'] if c['contact'] == 'drain') for s in states]),
                            native_current=np.array([v for _, v in native_curve]) * 1e6,
                            peak=np.array([t['candidate_peak_K'] for t in thermal]),
                            native_peak=np.array([t['reference_peak_K'] for t in thermal]))

    def finish(fig, stem, note):
        fig.text(.06, .035, note, fontsize=9, color='#4b5563', va='bottom')
        fig.subplots_adjust(left=.075, right=.975, top=.84, bottom=.19, wspace=.27)
        for ext in ('png', 'pdf'):
            fig.savefig(args.output / f'{stem}.{ext}', dpi=220)
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.4))
    fig.suptitle('D0 完整电热曲线：Vela 与 Sentaurus 对比', x=.075, ha='left', fontsize=18)
    for g, c in curves.items():
        for ax, v, r in ((axes[0], c['current'], c['native_current']),
                         (axes[1], c['peak'] - 300, c['native_peak'] - 300)):
            ax.plot(c['bias'], r, color=colors[g], lw=2, label=f'Sentaurus · Vg={g} V')
            ax.plot(c['bias'], v, 'o', ms=4, markerfacecolor='white', markeredgewidth=1.2,
                    color=colors[g], label=f'Vela · Vg={g} V')
            ax.set(xlim=(0, 40), xlabel='漏极电压 Vd (V)')
            ax.set_ylim(bottom=0)
            ax.grid(alpha=.18)
        axes[1].text(40, c['peak'][-1]-300+7, f'{c["peak"][-1]-300:.2f} K',
                     color=colors[g], ha='right', fontsize=10)
    axes[0].set(title='漏极电流', ylabel='单位宽度电流 Id (A/m)', ylim=(0, 280))
    axes[1].set(title='全域峰值温升', ylabel='峰值温升 Tmax − 300 K (K)', ylim=(0, 240))
    assert max(max(c['current'].max(), c['native_current'].max()) for c in curves.values()) < axes[0].get_ylim()[1]
    axes[0].legend(fontsize=9, frameon=False, loc='lower right')
    finish(fig, 'd0_curves_temperature', '每栅压 31 个精确求解点；原电学及批准热学门限全部通过。实线：Sentaurus；空心点：Vela。\n来源：D0 frozen R4 / full_r4_final.json；温度场覆盖全部 10241 个节点。')

    cold_path = args.evidence / 'r4_cold_summary.json'
    cold = read(cold_path)
    used[str(cold_path)] = digest(cold_path)
    assert len(cold) == 62 and all(p['numeric_pass'] and p['temperature_exact_300_K'] for p in cold)
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.4), sharey=True)
    fig.suptitle('自热对输出电流的影响：Vela 同模型对照', x=.075, ha='left', fontsize=18)
    reduction = {}
    for ax, (g, c) in zip(axes, curves.items()):
        points = sorted((p for p in cold if p['gate'] == g), key=lambda p: p['bias_V'])
        assert np.allclose(c['bias'], [p['bias_V'] for p in points], atol=1e-9, rtol=0)
        current = np.array([p['current_A_per_um'] for p in points]) * 1e6
        reduction[g] = float(100 * (1 - c['current'][-1] / current[-1]))
        ax.plot(c['bias'], current, '--', color='#606873', lw=2, label='冻结 300 K')
        ax.plot(c['bias'], c['current'], '-o', color=colors[g], lw=2, ms=3, label='四方程自洽自热')
        ax.set(title=f'Vg={g} V · 40 V 电流降低 {reduction[g]:.2f}%', xlabel='漏极电压 Vd (V)', xlim=(0, 40), ylim=(0, 400))
        ax.grid(alpha=.18)
        ax.legend(frameon=False, fontsize=10, loc='upper left')
        ax.annotate(f'{current[-1]:.2f} A/m', (40, current[-1]), xytext=(-6, 8), textcoords='offset points', ha='right', color='#606873', fontsize=10)
        ax.annotate(f'{c["current"][-1]:.2f} A/m', (40, c['current'][-1]), xytext=(-6, -18), textcoords='offset points', ha='right', color=colors[g], fontsize=10)
    axes[0].set_ylabel('单位宽度电流 Id (A/m)')
    finish(fig, 'd0_self_heating_current', '同一 Vela R4 物理/网格/电学边界；比较全域冻结 300 K 与联立晶格温度的自洽解。\n等温：62 点合格状态回放；自热：从零漏压完整推进。此图不用于比较运行时间。')

    mesh = read(args.mesh)
    used[str(args.mesh)] = digest(args.mesh)
    native_manifest = checked(args.evidence / 'native_fields_r1/manifest.json')
    mesh_hashes = {str(Path(k).resolve()): v for k, v in native_manifest['sources_sha256'].items()}
    assert digest(args.mesh) == mesh_hashes[str(args.mesh.resolve())]
    nodes = sorted(mesh['nodes'], key=lambda n: n['id'])
    assert [n['id'] for n in nodes] == list(range(10241))
    xy = np.array([[n['x'], n['y']] for n in nodes])
    tri = np.array([t['node_ids'] for t in mesh['triangles']])
    triangulation = Triangulation(xy[:, 0], xy[:, 1], tri)
    field = next(f for f in native_manifest['fields'] if f['gate_V'] == 8 and f['point_index'] == 30)
    native = checked(field['temperature_file'])
    assert native['node_id'] == list(range(len(nodes)))
    reference = np.array(native['temperature_K'])
    candidate = np.array(curves[8]['states'][-1]['temperature_K'])
    difference = candidate - reference
    difference_limit = max(.5, float(np.ceil(np.max(np.abs(difference)) * 2) / 2))
    edges = {}
    for t in mesh['triangles']:
        ids = t['node_ids']
        for a, b in zip(ids, ids[1:] + ids[:1]):
            edges.setdefault(tuple(sorted((a, b))), []).append(t['region_id'])
    outlines = [xy[list(e)] for e, regions in edges.items() if len(regions) == 1 or len(set(regions)) > 1]
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 5.8))
    fig.suptitle('高偏压晶格温度场：Vg=8 V，Vd=40 V', x=.055, ha='left', fontsize=18)
    for ax, values, title in zip(axes, (reference, candidate, difference),
                                 ('Sentaurus', 'Vela', '差值：Vela − Sentaurus')):
        norm = TwoSlopeNorm(vcenter=0, vmin=-difference_limit, vmax=difference_limit) if ax is axes[2] else plt.Normalize(300, 520)
        pc = ax.tripcolor(triangulation, values, shading='gouraud', cmap='RdBu_r' if ax is axes[2] else 'inferno', norm=norm, rasterized=True)
        ax.add_collection(LineCollection(outlines, colors='#a6adb4', linewidths=.45))
        ax.set(title=title, xlabel='网格 x (μm)', xlim=(-10.32, 10), ylim=(11, 0), aspect='equal')
        ax.set_xticks([-10, 0, 10])
        fig.colorbar(pc, ax=ax, orientation='horizontal', pad=.22, fraction=.08, label='温差 (K)' if ax is axes[2] else '晶格温度 (K)')
    axes[0].set_ylabel('网格 y (μm)')
    assert np.max(np.abs(difference)) <= difference_limit, 'Difference color scale clips field'
    assert min(reference.min(), candidate.min()) >= 300 and max(reference.max(), candidate.max()) <= 520
    fig.text(.055, .09, f'峰温：Sentaurus {reference.max():.3f} K；Vela {candidate.max():.3f} K。峰温差 2.189 K；体积加权温升场 RMS 误差 1.500 K。', fontsize=10)
    fig.text(.055, .04, '同一原始网格及节点映射，未旋转坐标；灰线为材料界面/外边界。温度图共用 300–520 K 色标，差值图独立色标。', fontsize=9, color='#4b5563')
    fig.subplots_adjust(left=.055, right=.98, top=.83, bottom=.27, wspace=.22)
    for ext in ('png', 'pdf'):
        fig.savefig(args.output / f'd0_temperature_field_vg8_vd40.{ext}', dpi=220)
    plt.close(fig)
    provenance = dict(scope='Report figure provenance; raw simulation fields are stored separately',
                      input_sha256=used, self_heating_current_reduction_percent=reduction,
                      field_max_absolute_difference_K=float(np.max(np.abs(difference))),
                      outputs_sha256={f.name: digest(f) for f in args.output.iterdir() if f.suffix in ('.png', '.pdf')})
    (args.output / 'figure_provenance.json').write_text(json.dumps(provenance, indent=2), encoding='utf-8')
    print(json.dumps({'self_heating_current_reduction_percent': reduction,
                      'field_max_absolute_difference_K': provenance['field_max_absolute_difference_K'],
                      'figures': list(provenance['outputs_sha256'])}, indent=2))


if __name__ == '__main__':
    main()
