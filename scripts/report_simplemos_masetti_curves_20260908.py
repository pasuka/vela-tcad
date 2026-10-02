"""Report current errors and solve qualification separately; no interpolation."""
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import simplemos_masetti_curve_vela_20260908 as v

a, d, OUT, LOCAL, REPO = v.a, v.d, v.OUT, v.LOCAL, v.REPO
REPORT = REPO/'docs/validation/simplemos_masetti_fullcurve_validation_2026-09-08.md'


def main():
    a.verify(OUT/'comparison_evidence.json')
    rows = a.rows(OUT/'comparison.csv')
    summary = a.read(OUT/'summary.json')
    refs = a.rows(OUT/'native_points.csv')
    mapping_checks = []
    for c in a.read(OUT/'vela_contract.json')['cases']:
        geo, _ = v.V.m.previous.prior.support(c)
        for index in range(51):
            export = LOCAL/'native_exports'/c['case']/f'vg_{index:03d}'
            error = v.mapping.original.m73.coordinate_error(export, geo)
            assert math.isfinite(error) and error<=1e-12
            fn = v.mapping.original.m73.scalar(export/'fields/eQuasiFermiPotential_region0.csv')
            fp = v.mapping.original.m73.scalar(export/'fields/hQuasiFermiPotential_region0.csv')
            assert set(fn)==set(fp) and len(fn)==942
            mapping_checks.append(dict(case=c['case'], index=index, coordinate_error_um=error, semiconductor_nodes=len(fn)))
    a.write_csv(OUT/'mapping_checks.csv', mapping_checks)
    staged_evidence = []
    for path in (OUT/'staged').glob('*/evidence.json'):
        a.verify(path)
        case = path.parent.name
        for source in (LOCAL/'stage'/case).iterdir():
            if source.is_file():
                assert a.sha(source)==a.sha(LOCAL/'native_raw/bundle'/case/source.name)
        staged_evidence += [path, path.parent/'freeze.json']
    for device in ('n19', 'n23'):
        for vd in (.05, 1.):
            points = [r for r in rows if r['device']==device and float(r['vd'])==vd]
            assert len(points)==51
            assert sorted(float(r['vg']) for r in points)==v.n.GRID
            for r in points:
                for arm in ('continuation', 'native'):
                    error = 100*(float(r[arm+'_Id_A_per_um'])/float(r['sentaurus_Id_A_per_um'])-1)
                    assert abs(error-float(r[arm+'_error_percent']))<=1e-12
    attempts = [r for arm in ('continuation', 'native') for r in a.rows(OUT/(arm+'_attempts.csv'))]
    accepted = [r for r in attempts if r['qualified']=='True']
    failures = [r for r in attempts if r['qualified']!='True']
    curve_rows, bands, anchors = [], [], []
    old = a.rows(v.old.OUT/'current_extraction_fix/comparison.csv')
    for c in summary['curves']:
        subset = [r for r in rows if r['device']==c['device'] and float(r['vd'])==c['vd']]
        qualified = [r for r in subset if r['comparison_qualified']=='True']
        curve_rows.append(f"| {c['device']} | {c['vd']:g} | {c['qualified']}/51 | {c['error_min_percent']:+.8f}% | {c['error_max_percent']:+.8f}% | {c['worst_vg']:g} |" if qualified else f"| {c['device']} | {c['vd']:g} | 0/51 | — | — | — |")
        for lo, hi in ((0., .5), (.5, 1.)):
            selected = [r for r in qualified if lo<=float(r['vg'])<=hi and (lo==0 or float(r['vg'])>lo)]
            if selected:
                worst = max(selected, key=lambda r:abs(float(r['continuation_error_percent'])))
                bands.append(dict(device=c['device'], vd=c['vd'], lower_V=lo, upper_V=hi, qualified=len(selected),
                                  max_absolute_error_percent=abs(float(worst['continuation_error_percent'])), worst_vg=float(worst['vg'])))
        for vg in (.8, 1.):
            point = next(r for r in subset if float(r['vg'])==vg)
            previous = next(r for r in old if r['device']==c['device'] and float(r['vd'])==c['vd'] and float(r['vg'])==vg and r['arm']=='generated')
            previous_id = float(previous['Id_A_per_um'])
            previous_native = previous_id/(1+float(previous['native_relative_error']))
            anchors.append(dict(device=c['device'], vd=c['vd'], vg=vg,
                vela_relative_drift=float(point['continuation_Id_A_per_um'])/previous_id-1,
                native_relative_drift=float(point['sentaurus_Id_A_per_um'])/previous_native-1,
                comparison_qualified=point['comparison_qualified']=='True'))
    n23_table = []
    for vg in (0., .02, .1, .2, .3, .5, .8, 1.):
        points = [next(r for r in rows if r['device']=='n23' and float(r['vd'])==vd and float(r['vg'])==vg) for vd in (.05, 1.)]
        cells = [f"{float(r['continuation_error_percent']):+.8f}%" + ('' if r['comparison_qualified']=='True' else '（未合格诊断）') for r in points]
        n23_table.append(f'| {vg:g} | {cells[0]} | {cells[1]} |')
    deep = next(r for r in rows if r['device']=='n23' and float(r['vd'])==1. and float(r['vg'])==0.)
    higher = [r for r in rows if r['device']=='n23' and float(r['vg'])>=.3 and r['comparison_qualified']=='True']
    a.write_csv(OUT/'voltage_bands.csv', bands)
    a.write_csv(OUT/'anchor_replays.csv', anchors)
    v.csv_union(OUT/'failed_attempts.csv', failures or [dict(note='No failed DC attempts')])
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 7.8), sharex='col', layout='constrained')
    for column, device in enumerate(('n19', 'n23')):
        for vd, color in ((.05, '#1764ab'), (1., '#cc5833')):
            subset = sorted([r for r in rows if r['device']==device and float(r['vd'])==vd], key=lambda r:float(r['vg']))
            x = np.array([float(r['vg']) for r in subset])
            ref = np.array([float(r['sentaurus_Id_A_per_um']) for r in subset])
            current = np.array([float(r['continuation_Id_A_per_um']) for r in subset])
            error = np.array([float(r['continuation_error_percent']) for r in subset])
            good = np.array([r['comparison_qualified']=='True' for r in subset])
            axes[0,column].semilogy(x, np.abs(ref), color=color, label=f'Sentaurus, Vd={vd:g} V')
            axes[0,column].semilogy(x[good], np.abs(current[good]), 'o', ms=3, mfc='none', color=color, label=f'Vela qualified, Vd={vd:g} V')
            axes[1,column].plot(x, np.where(good, error, np.nan), '.-', color=color, label=f'Vd={vd:g} V')
            if (~good).any():
                axes[1,column].plot(x[~good], error[~good], 'x', color=color, label='Unqualified diagnostic')
        axes[0,column].set_title(device+' | 300 K, Masetti, automatic box')
        axes[0,column].set_ylabel('|Id| [A/um]')
        axes[1,column].set_ylabel('100 (Id Vela / Id Sentaurus - 1) [%]')
        axes[1,column].set_xlabel('Vg [V]')
        axes[1,column].axhline(0, color='.65', lw=.7)
        for ax in axes[:,column]:
            ax.grid(alpha=.2)
            ax.legend(fontsize=8)
    fig.suptitle('Matched-model Id-Vg comparison | 51 targets per curve')
    fig.savefig(LOCAL/'idvg_comparison.png', dpi=180)
    fig.savefig(LOCAL/'idvg_comparison.svg')
    plt.close(fig)
    q = [r for r in rows if r['comparison_qualified']=='True']
    maxfield = lambda key:max((float(r[key]) for r in q), default=math.nan)
    lines = [
        '# SimpleMOS Masetti 完整 0–1 V 曲线验证', '',
        '日期：2026-09-08。延续已通过八点验证的显式生产组合，补齐 n19/n23 × Vd=0.05/1 V，Vg=0–1 V、步长 0.02 V，每条 51 点。', '',
        f"原生参考 {summary['native_qualified']}/204 点合格；Vela 延续路径 {summary['continuation_qualified']}/204、原生状态独立初始化 {summary['independent_qualified']}/204；双初始化和参考共同合格 **{summary['comparison_qualified']}/204**。", '',
        '## 电流误差', '',
        '统一定义为 100×(Id_Vela/Id_Sentaurus−1)，正值表示 Vela 偏大。以下统计只包含双初始化及原生参考均合格点。没有把小扰动预测门槛当成绝对电流误差门槛，也没有新增或放宽电流百分比验收线。', '',
        '| NWell | Vd/V | 完整资格 | 最小有符号差 | 最大有符号差 | 最大绝对差的 Vg/V |',
        '| --- | ---: | ---: | ---: | ---: | ---: |', *curve_rows, '',
        '![Id–Vg 和逐点误差](../../build-release/simplemos_masetti_curves_20260908/idvg_comparison.png)', '',
        '图中的叉号如存在，仅表示未合格点的诊断数值；不会进入合格误差统计。原生和 Vela 使用同一偏置逐点直接比较，没有插值。', '',
        '### n23 的电压依赖', '',
        '| Vg/V | Vd=0.05 V | Vd=1 V |', '| ---: | ---: | ---: |', *n23_table, '',
        f"n23、Vd=1 V、Vg=0 V：Sentaurus Id={float(deep['sentaurus_Id_A_per_um']):.12e} A/μm，Vela Id={float(deep['continuation_Id_A_per_um']):.12e} A/μm，绝对有符号差={float(deep['continuation_Id_A_per_um'])-float(deep['sentaurus_Id_A_per_um']):.9e} A/μm，相对差={float(deep['continuation_error_percent']):+.8f}%。该点双初始化资格为 {deep['dual_qualified']}。绝对量很小与相对差接近 0.83% 同时成立，应分别报告。", '',
        f"n23 的 Vg≥0.3 V 合格点中，两个 Vd 的最大绝对相对差为 {max((abs(float(r['continuation_error_percent'])) for r in higher),default=math.nan):.8f}%。本轮确认剩余较明显的相对差集中于高 Vd 的深关断端；尚未完成该漏电差的原因归属。", '',
        '## 严格资格与失败保留', '',
        f"实际 Vela DC 尝试 {len(attempts)} 次，合格 {len(accepted)} 次，失败记录 {len(failures)} 次。失败后最多在同偏置重载一次，原记录保留；只有合格状态才可作为下一个偏置的延续种子。", '',
        f"合格 DC 中最大逐行比值 {max(float(r['max_row_ratio']) for r in accepted):.9g}（门槛 1e-6），最大 KCL/Id {max(float(r['kcl_over_Id']) for r in accepted):.9g}（门槛 1e-8），最大端口提取/守恒功能量相对差 {max(float(r['port_relative']) for r in accepted):.9g}（门槛 1e-8）。", '',
        f"合格配对最大 ψ/φn/φp 差分别为 {maxfield('psi_max_V'):.9g}/{maxfield('phin_max_V'):.9g}/{maxfield('phip_max_V'):.9g} V（各 1e-6 V）；密度相对差 {maxfield('density_max_relative'):.9g}（1e-4），两条 Vela 路径的 Id 相对差 {maxfield('dual_Id_relative'):.9g}（1e-6）。", '',
        '全部 1814 个自由 Si 载流子行逐行检查，无密度/通量筛选或零尺度排除。全局连续性闭合沿用原守恒接触边通量对 SRH 源积分、1e-6 门槛和 1e-10 源下限；低于源下限不声明任意相对精度。', '',
        '## 模型、网格和运行条件', '',
        'n19/n23 的 NWell Boron 为 1e17/2e17 cm⁻³，1480/1482 节点、2742/2746 三角形，网格坐标 μm，密度 m⁻³，宽度 1 μm，电流 A/μm。300 K，Boltzmann/no-BGN、matched ni、总杂质 Masetti、掺杂 SRH；表面迁移率、高场饱和、Auger、雪崩及 DG 关闭。', '',
        '显式组合为 Delaunay box 转移、cell_material 介电装配、element_box 输运和迁移率、signed_transport 三项 Poisson 电荷体积。生产现代 SI 常数、原 SRH 体积、原标量线搜索保持不变。未启用手册常数、稳定范数比较或 ψ 余量累积诊断。', '',
        'Vela 使用现有 UCRT64 Release 生产 runner、Eigen SparseLU/COLAMD、4 次线性迭代修正；max_iter=200、reltol=1e-7、abstol=1e-12、原 stall floor、ψ/准费米更新上限 0.35/0.025 V、contact_basin。原生为 Sentaurus T-2022.03-SP2、ExtendedPrecision(128)、Method=Super、Digits=12、ErrRef(e/h)=1e-2、RhsMin=1e-20、Iterations=40。', '',
        '原生从既有 Masetti 0.8 V 保存态回闭合，下降到零后向上扫描。Vela 第一条路径从既有合格 0.8 V 状态向下、向上分别延续；第二条路径每个偏置独立使用原生势和准费米势、以 Vela 未改的 ni/Vt 重算一致初始密度，再完成自洽求解。原生密度导出另行保留，未混称为重算密度。', '',
        '虚拟机有 4 核。运行中以独立冻结的调度补充停止串行调度父进程，保留其正在运行的 n19 DC 子进程，并行启动尚未开始的两条 n23 曲线；仿真指令和输入逐字节不变。已完成原生曲线可提前取回并执行同一套独立初始化任务；提前取回的全部文件已与最终完整原始包逐字节核对。', '',
        '## 覆盖范围与复核', '',
        '这批结果覆盖当前简化模型的四条完整曲线。PhuMob/Lombardi/HFS 恢复、其他控制组合、全物理场一致性及新增 SRH 局部扰动仍是独立验证范围。先前四个 n23 点的已通过结论保留。', '',
        '后续可优先在 n23、Vd=1 V、Vg=0–0.2 V 比较绝对漏电及守恒生成复合源的分量，沿原有 SRH 体积校准证据区分候选；本轮没有据此修改模型，也不以当前差值直接认定 SRH 为唯一原因。', '',
        '新增验证脚本的 4 项单元检查通过，覆盖双路径必须合格、少数载流子场门槛、非有限数拒绝和扫描偏置覆盖。本轮未修改 C++ 求解器或默认参数，未重复运行此前已有历史失败的全套 CTest。', '',
        f"既有八点重放最大 Vela Id 相对漂移 {max(abs(r['vela_relative_drift']) for r in anchors):.9g}，新旧原生参考最大漂移 {max(abs(r['native_relative_drift']) for r in anchors):.9g}；详见锚点表。", '',
        '机器证据和原始输出保留在本地，未提交仿真输出；所有失败、输入哈希及远端原始包均保留。执行前冻结原生合同与 Vela 合同，结果另行封存，旧合同未改写。', '',
        '- [逐点电流及资格](../../reference_tcad/simplemos_sentaurus2022/masetti_curves_20260908/comparison.csv)',
        '- [分电压区间统计](../../reference_tcad/simplemos_sentaurus2022/masetti_curves_20260908/voltage_bands.csv)',
        '- [既有八点重放](../../reference_tcad/simplemos_sentaurus2022/masetti_curves_20260908/anchor_replays.csv)',
        '- [失败尝试](../../reference_tcad/simplemos_sentaurus2022/masetti_curves_20260908/failed_attempts.csv)',
        '- [原生准备与导出脚本](../../scripts/simplemos_masetti_curve_native_20260908.py)',
        '- [Vela 执行脚本](../../scripts/simplemos_masetti_curve_vela_20260908.py)',
        '- [汇总脚本](../../scripts/report_simplemos_masetti_curves_20260908.py)',
        '- [最终证据](../../reference_tcad/simplemos_sentaurus2022/masetti_curves_20260908/final_evidence.json)', '',
    ]
    REPORT.write_text('\n'.join(lines), encoding='utf-8')
    d.matrix.freeze(OUT/'final_evidence.json', [OUT/'comparison_evidence.json', REPORT, Path(__file__).resolve(),
        REPO/'tests/regression/test_simplemos_masetti_curves.py', OUT/'voltage_bands.csv', OUT/'anchor_replays.csv',
        OUT/'failed_attempts.csv', OUT/'mapping_checks.csv', OUT/'parallel_dispatch_freeze.json',
        OUT/'parallel_dispatch_addendum.json', LOCAL/'parallel_finish.sh',
        LOCAL/'idvg_comparison.png', LOCAL/'idvg_comparison.svg', LOCAL/'unit_tests.log'] + staged_evidence)
    print(REPORT)


if __name__ == '__main__':
    main()
