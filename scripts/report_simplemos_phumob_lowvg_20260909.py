"""Seal script QA and render the bounded low-Vg validation report."""
import argparse
import ast
import math
import re
import subprocess
from pathlib import Path

import validate_simplemos_phumob_lowvg_20260909 as s

a, d, OUT, LOCAL = s.a, s.d, s.OUT, s.LOCAL
SCRIPT_NAMES = ('validate_simplemos_phumob_lowvg_20260909.py', 'analyze_simplemos_phumob_lowvg_20260909.py',
    'review_simplemos_phumob_lowvg_20260909.py', 'audit_simplemos_phumob_lowvg_export_precision_20260909.py',
    'audit_simplemos_phumob_lowvg_srh_20260909.py', 'audit_simplemos_phumob_lowvg_srh_fields_20260909.py',
    'report_simplemos_phumob_lowvg_20260909.py')


def qa():
    scripts = [s.REPO / 'scripts' / n for n in SCRIPT_NAMES]
    for path in scripts:
        ast.parse(path.read_text(encoding='utf8'))
    tests = []
    root = LOCAL / 'qa'
    root.mkdir(exist_ok=False)
    for name, count in [('test_simplemos_bgn_restore.py', 5), ('test_simplemos_masetti_curves.py', 4)]:
        path = s.REPO / 'tests/regression' / name
        result = subprocess.run(['D:/msys64/ucrt64/bin/python.exe', '-X', 'utf8', str(path)],
                                cwd=s.REPO, capture_output=True, text=True)
        text = result.stdout + result.stderr
        (root / (name + '.log')).write_text(text, encoding='utf8')
        passed = result.returncode == 0 and f'Ran {count} tests' in text and '\nOK' in text
        tests.append(dict(test=name, count=count, exit_code=result.returncode, passed=passed))
        assert passed, text
    a.write(OUT / 'qa_summary.json', dict(python_syntax_checks=len(scripts), tests=tests,
        regression_tests_passed=sum(r['count'] for r in tests), production_source_changed=False,
        full_ctest_repeated=False, previous_full_ctest='772/788; same 16 historical failures. Previous stage only, not a new execution.'))
    a.write(OUT / 'srh_failure_ledger.json', dict(
        initial_failure='Initial SRH manifest lookup used SRHRecombination; actual field name is srhRecombination. Frozen initial inputs and script retained; no SRH result admitted from that run.',
        amendment='Separate srh_fields output uses the exact exported name. An initial output-prefix typo in its input-manifest path failed before freezing and was corrected. Physics, fields, units and thresholds unchanged.',
        original_manifest=str(OUT / 'srh_freeze.json'), corrected_manifest=str(OUT / 'srh_fields/srh_freeze.json')))
    d.matrix.freeze(OUT / 'qa_evidence.json', scripts + [OUT / 'qa_summary.json', OUT / 'srh_failure_ledger.json',
        OUT / 'srh_freeze.json', OUT / 'srh_fields/srh_evidence.json'] + list(root.glob('*.log')))
    print(a.read(OUT / 'qa_summary.json'), flush=True)


def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |', '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                     ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows])


def report():
    for path in ('review_evidence.json', 'qa_evidence.json', 'srh_fields/srh_evidence.json'):
        a.verify(OUT / path)
    summary = a.read(OUT / 'review_summary.json')
    post = a.read(OUT / 'post_summary.json')
    cal = a.read(OUT / 'calibration_summary.json')
    precision = a.read(OUT / 'export_precision_summary.json')
    low = a.rows(OUT / 'dc_comparison.csv')
    attempts = [r for r in a.rows(OUT / 'all_attempts.csv') if r['origin'] == 'low']
    failures = [r for r in attempts if r['qualified'] != 'True']
    currents = []
    for r in sorted(low, key=lambda x: (x['device'], float(x['vd']), float(x['vg']))):
        currents.append([r['device'], r['vg'], r['vd'], f"{float(r['native_Id_A_per_um']):.10e}",
            f"{float(r['legacy_versus_native_percent']):+.8f}%", f"{float(r['candidate_versus_native_percent']):+.8f}%",
            r['legacy_dual_qualified'], r['candidate_dual_qualified']])
    srh = a.rows(OUT / 'srh_fields/srh_summary.csv')
    srhtable = [[r['device'], r['vg'], r['vd'], f"{100*float(r['integrated_source_relative']):+.7f}%",
        f"{float(r['rate_contribution_charge_equivalent_A_per_um']):+.6e}",
        f"{float(r['volume_contribution_charge_equivalent_A_per_um']):+.6e}"] for r in srh]
    fields = a.rows(OUT / 'native_field_comparison.csv')
    fieldtable = []
    for point in sorted([r for r in low if r['device'] == 'n23'], key=lambda r: (float(r['vd']), float(r['vg']))):
        group = {r['field']: r for r in fields if r['key'] == point['key'] and r['stage'] == 'candidate'}
        fieldtable.append([point['vg'], point['vd']] + [f"{float(group[k]['max_absolute']):.6e}" for k in ('psi', 'phin', 'phip')] +
                          [f"{float(group[k]['max_relative']):.6e}" for k in ('electrons_m3', 'holes_m3')])
    srhnodes = a.rows(OUT / 'srh_fields/srh_nodes.csv')
    hot = sorted([r for r in srhnodes if r['key'] == 'm65_n23_vd_1p000000_endpoint_vg_000'],
                 key=lambda r: abs(float(r['volume_contribution_per_m_s'])), reverse=True)[:8]
    hottest = ', '.join(r['node_id'] for r in hot)
    failuretable = [[r['stage'], r['device'], r['vg'], r['vd'], r['arm'], r['attempt'],
                     r['iterations'], f"{float(r['max_row_ratio']):.6e}", r['failure']] for r in failures]
    body = f'''# SimpleMOS PhuMob 低栅压与剩余源项验证

日期：2026-09-09。分支：`codex/simplemos-sdevice-validation`。前置结果见 [PhuMob 单元平均八点验证](simplemos_phumob_box_candidate_validation_2026-09-09.md)。

**已补齐 Vg=0/0.2 V 的八个低栅压点及两种初始化，并完成原生单元迁移率、导出精度和 SRH 积分账本。** 与已有 Vg=0.8/1 V 八点合并后，候选双初始化合格 {summary['qualified_points']['candidate']}/16，legacy 合格 {summary['qualified_points']['legacy']}/16。n23、Vd=1 V、Vg=0 V 的候选 Id 差仍为 **−0.80813524%**，绝对差约 **−3.00069e-18 A/μm**，不能由高栅压误差小而放行整个曲线。

SRH 体积差的电荷等价值约为 +3.00624e-18 A/μm，与上述差的绝对量级接近；同体积速率差仅约 −6.15312e-22 A/μm。这是优先级很高的源项候选，**尚不是经同源扰动校准的漏端电流因果归因**。本轮未修改生产 C++、模型默认、物理常数、SRH 体积或接受条件。

## 1. 范围与输入资格

同一 SimpleMOS Si/SiO2 网格：n19 为 1480 节点、n23 为 1482 节点；每点独立检查 907 个自由 Si 节点、1814 个载流子行。Windows UCRT64 Release，UMFPACK；沿用上轮已冻结的 runner 和静态库。Sentaurus 使用 T-2022.03-SP2，原生独立目录为 `{s.REMOTE}`，16 次 DC（8 个 PhuMob、8 个同偏置 Masetti 控制）已结束、取回并全部通过版本、偏置和 KCL 检查。16 份 TDR 场导出完成。

两端 plain PhuMob、OldSlotboom、独立匹配的基础 ni、300 K、Boltzmann、掺杂相关 SRH；Enormal/HFS 仍关闭。Vela 两组只切换 `legacy` 与显式 `element_box_phumob`；配套输运、Poisson 电荷及介电几何保持原组合。输入采用 `unit_scaling`；电流 A/μm、电势 V，场表密度 m⁻³。原生迁移率比较使用 cm²/(V·s)，SRH 速率使用 cm⁻³·s⁻¹；转换在脚本中显式执行。

沿用原两条初始化路径和最多一次同偏置重载。逐行 1e-6、KCL/Id 1e-8、端口一致性 1e-8、双初始化最大势差 1e-6 V、密度相对差 1e-4、Id 相对差 1e-6 均未改变。全局源的相对条件继续保留原 1e-10 下限；低于下限不算独立的相对源闭合通过。

## 2. 低栅压电流与初始化

误差定义为 `(Id_Vela/Id_Sentaurus−1)×100%`。表中 Vela 电流选自合格的 Vela 初始化；两条路径的最终资格单独列出，未通过的比较不能按正式双初始化合格点使用。

{table(['NWell','Vg (V)','Vd (V)','原生 Id (A/μm)','legacy 误差','候选误差','legacy 双初始化','候选双初始化'], currents)}

新增低栅压实际尝试 {summary['low_attempts']} 次，其中 {summary['low_failed_attempts']} 次未通过，全部保留。首次或重载失败如下：

{table(['组','NWell','Vg','Vd','初始化','attempt','迭代','最大行比','失败原因'], failuretable)}

合并 16 点的最大绝对相对误差：legacy {summary['max_absolute_Id_error_percent']['legacy']:.8f}%，候选 {summary['max_absolute_Id_error_percent']['candidate']:.8f}%；{summary['improved_points']}/16 点的误差幅度减小。资格仍以逐点标记为准。完整值与高低 NWell 配对见 [16 点账本](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/sixteen_points.csv)、[配对](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/nwell_pairs.csv)、[尝试记录](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/all_attempts.csv)。这些是四个 Vg 的控制点，未重算 51 点完整曲线。

## 3. Jacobian、少数载流子与场差

八个合格 Vela 初始化候选状态恢复 SRH 后，各用 ψ/φn/φp 方向和四个幅度检查 Jv。较小三幅度的 {post['gated_checks']} 个正式分块全部通过，最大相对差 {post['max_relative']:.8e}，门槛 1e-4；弱电子/空穴交叉源块仍单列，不用全残差差分代替独立弱源校准。该导数检查独立于另一初始化是否通过。

候选选取态的最大载流子行比为 {summary['candidate_worst_row']:.8e}。所有高低点按载流子分别保存最差五行、局部密度、是否为少数载流子、原始残差、绝对边通量及 SRH 项；源下限激活 {summary['source_floor_active']}/{summary['source_components']}。见 [行排序](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/minority_row_ranking.csv)、[源下限与未截底比值](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/source_support.csv)、[Jv](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/post_summary.json)。

n23 候选在自由 Si 节点上的最大原生场差如下；密度列为无量纲最大相对差，不是百分数：

{table(['Vg','Vd','max |Δψ| (V)','max |Δφn| (V)','max |Δφp| (V)','max 相对 Δn','max 相对 Δp'], fieldtable)}

完整场最大值、RMS 和最差节点见 [场差账本](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/native_field_comparison.csv)。

## 4. G 截断与原生端口重放的边界

56000 个单元/载流子比较及 59360 个边系数比较完成。Masetti 单元公式与原生最大相对差 4.44e-16，八个 Masetti 原生电流控制均在原 1e-8 条件内重闭合。PhuMob 未截断单元最大差 2.22e-16；电子最大差约 1.38265e-9；空穴最大差 {cal['exact_phumob_cell_max_relative']:.8e}，2960 个含截断顶点的空穴单元仍未通过 1e-7。

使用上轮高 Vg 冻结的两个有效 G 下限作为诊断，新增低 Vg 单元最大差 {cal['heldout_phumob_cell_max_relative']:.8e}。本轮拟合参数数为零。该结果验证已定位的差异形式能预测新状态，仍不能确定原生内部极小值搜索或截断算法，也未将反推参数写入生产。

低 Vg 的 16 个原生准费米势 SG 端口重放均未通过 1e-6，最大相对差 {cal['native_replay_max_relative']:.8f}。失败同时出现于 Masetti 与 PhuMob，不能将原生直接输出 Id 的资格与导出势场重放资格混用。

对每个端口先按节点合并有符号准费米势敏感度，再估计 binary64 最近舍入的半 ULP 影响。估计范围为参考电流的 {precision['min_half_ulp_bound_over_Id']:.8g}～{precision['max_half_ulp_bound_over_Id']:.8g} 倍，16 个重放偏差均在该范围内。正负一个 ULP 的直接重放与线性预测最大差 {precision['max_one_ulp_linearization_relative']:.8e}。这支持导出舍入足以解释量级；它不是未知原生高精度值的恢复或因果证明。上述失败不能调宽门槛后计为通过，固定态端口替换量也不能直接当作 Id 归因。

证据：[单元分组](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/native_groups.csv)、[校准汇总](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/calibration_summary.json)、[端口重放](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/native_ports.csv)、[导出精度](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/export_precision.csv)。

## 5. SRH 速率和体积的独立账本

使用原生 runtime 单元顶点 box 测度形成 Si 源体积；Vela 保留生产 all-cell SRH 体积。生产行积分通过同一边探针的物理粒子通量/内部通量比例换回物理单位，再除以当前体积取得 SRH 速率。比较只覆盖同一批自由 Si 节点，分解采用恒等式 `R_V V_V − R_S V_S = (R_V−R_S)V_S + R_V(V_V−V_S)`；没有改变源体积或重算候选。

下表后两列为积分源差乘 q 后的电荷等价值，**不是已校准的漏端 ΔId**：

{table(['NWell','Vg','Vd','积分源相对差','同体积速率差 (A/μm)','体积差 (A/μm)'], srhtable)}

n23、Vd=1、Vg=0 的体积贡献集中节点排序为 {hottest}。源账本闭合到浮点舍入量级；生成源积分相对差约 −1.55715%，与漏端相对差并不相等。应先校准这些局部源到端口的响应，再考虑体积候选自洽修改。见 [逐节点源账本](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/srh_fields/srh_nodes.csv)、[源汇总](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/srh_fields/srh_summary.csv)。

初次 SRH 辅助审计的字段名大小写检查失败，原输入和脚本保留；修正版按实际 `srhRecombination` 字段名写入独立目录，未改变数据或单位。修正版开始前的输入路径错误也保留说明，见 [辅助审计失败记录](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/srh_failure_ledger.json)。

## 6. 验证与后续

新增脚本完成语法检查；既有 BGN 模型恢复 5 项和完整曲线驱动 4 项回归共 9 项通过。本轮没有修改 C++ 或重建二进制，未重复全量 CTest；前置阶段 772/788、16 项历史身份/缺失文件失败仍是上一轮记录，不能当成本轮新运行。输入、失败、原生归档、导出、数值账本和 QA 分别冻结，见 [复核汇总](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/review_summary.json)、[脚本 QA](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/qa_summary.json)。未提交或推送。

下一阶段优先在已合格状态上，以 n23 高 Vd 体积贡献热点及低 NWell 控制做局部 SRH 源/体积的双幅度正负响应校准；原生端应直接记录守恒电流或足够精度的源响应，避免依赖已失去小梯度分辨率的导出准费米势重放。只有同扰动响应、收敛、双初始化及端口守恒通过后，才开展有限体积替换。G 截断算法资格、弱源资格及初始化失败轨迹继续保留。Enormal、高场饱和和完整 0–1 V 曲线尚未恢复或放行。
'''
    path = s.REPO / 'docs/validation/simplemos_phumob_lowvg_validation_2026-09-09.md'
    assert not path.exists()
    path.write_text(body, encoding='utf8')
    d.matrix.freeze(OUT / 'report_evidence.json', [Path(__file__).resolve(), path, OUT / 'review_evidence.json', OUT / 'qa_evidence.json', OUT / 'srh_fields/srh_evidence.json'])
    print(path, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('qa', 'report'))
    globals()[parser.parse_args().action]()
