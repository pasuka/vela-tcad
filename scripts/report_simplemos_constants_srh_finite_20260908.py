"""Seal finite constants/SRH results without changing frozen experiments."""
import ast
import re
from pathlib import Path

import restart_simplemos_srh_finite_20260908 as restart

t = restart.t
c, a, d, v, p = t.c, t.a, t.d, t.v, t.p
OUT = t.OUT.parent
REPORT = p.REPO / 'docs/validation/simplemos_constants_and_srh_finite_validation_2026-09-08.md'


def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                      '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                     ['| ' + ' | '.join(str(x) for x in row) + ' |' for row in rows])


def main():
    manifests = [c.OUT / 'validation_evidence.json',
                 t.continuation.OUT / 'validation_evidence.json',
                 t.OUT / 'build_evidence.json', t.OUT / 'preflight_evidence.json',
                 t.OUT / 'response_evidence.json', t.OUT / 'validation_evidence.json',
                 t.OUT / 'input_tests_evidence.json', t.OUT / 'native_rows/evidence.json',
                 restart.OUT / 'validation_evidence.json']
    for path in manifests:
        a.verify(path)
    cases = a.read(t.OUT / 'contract.json')['cases']
    bykey = {case['key']: case for case in cases}
    comparison = a.rows(c.OUT / 'comparison.csv')
    direct_constants = a.rows(c.OUT / 'dc.csv')
    continuation = a.rows(t.continuation.OUT / 'dc.csv')
    constant_dual = a.rows(t.continuation.OUT / 'dual.csv')
    preflight = a.rows(t.OUT / 'preflight.csv')
    small = a.rows(t.OUT / 'small_dc.csv')
    response = a.rows(t.OUT / 'response.csv')
    original_finite = a.rows(t.OUT / 'finite_dc.csv')
    selected = a.rows(restart.OUT / 'selected_states.csv')
    dual = a.rows(restart.OUT / 'dual.csv')
    local = a.rows(restart.OUT / 'local_fields.csv')
    fields = a.rows(restart.OUT / 'fields.csv')
    failures = a.rows(restart.OUT / 'remaining_violations.csv')
    input_tests = a.rows(t.OUT / 'input_tests.csv')
    assert len(selected) == 24 and len(dual) == 12
    assert all(x['qualified'] == 'True' for x in preflight + small + response + constant_dual)
    assert len(input_tests) == 9 and all(x['qualified'] == 'True' for x in input_tests)
    assert all(x['qualified'] == 'True' for x in selected if x['axis'] == 'joint')
    assert all(x['qualified'] == 'True' for x in dual if x['axis'] == 'joint')

    # The original initial-batch dual files test distances only. Retain them,
    # but qualify final pairs using both DC qualifications as well as distances.
    ports, currents, baseline_local, baseline_fields = [], [], [], []
    for row in selected:
        case = bykey[row['key']]
        dest = p.REPO / row['path']
        st = a.read(dest / 'post/functional.status.json')
        current = float(row['current_A_per_um'])
        err = max(abs(current / st['current_A_per_um'] - 1),
                  abs(st['current_A_per_um'] / st['contact_current_extractor_A_per_um'] - 1))
        ports.append(dict(key=row['key'], axis=row['axis'], label=row['label'],
                          relative=err, functional_consistent=err <= 1e-8,
                          DC_qualified=row['qualified'] == 'True'))
        baseline = next(x for x in comparison if x['key'] == row['key'] and x['label'] == 'manual_constants')
        currents.append(dict(key=row['key'], axis=row['axis'], label=row['label'],
                             current_A_per_um=current, native_relative_error=current / case['native_Id_A_per_um'] - 1,
                             delta_over_constant_baseline=current / float(baseline['current_A_per_um']) - 1,
                             qualified=row['qualified'] == 'True'))
    assert all(x['functional_consistent'] for x in ports)
    for case in cases:
        nodes = [case['mapped_nodes'][tag]['node'] for tag in ('792', '1057')] + [1000, 1009]
        met, loc = c.field_metrics(case, 'manual_constants', c.LOCAL / case['key'] / 'manual_constants', nodes=nodes)
        baseline_fields += met
        baseline_local += loc

    traces = []
    for failed in failures:
        path = restart.LOCAL / failed['key'] / failed['axis'] / failed['label'] / 'updates.csv'
        rows = [x for x in a.rows(path) if x['node_id'] == failed['node_id'] and x['carrier'] == failed['carrier']]
        last = rows[-1]
        traces.append(dict(key=failed['key'], node_id=failed['node_id'], carrier=failed['carrier'],
                           final_row_ratio=failed['ratio'], final_Poisson_block=failed['psi_block'],
                           **{k: last[k] for k in ('iteration', 'line_search_accepted', 'line_search_attempts',
                               'selected_damping', 'best_rejected_candidate', 'residual', 'raw_linear_step_V',
                               'capped_step_V', 'applied_step_V', 'raw_linear_residual',
                               'selected_trial_residual', 'qf_update_limit_V')}))
    for name, rows in [('selected_ports', ports), ('selected_currents', currents),
                       ('constant_baseline_local', baseline_local), ('constant_baseline_fields', baseline_fields),
                       ('remaining_update_trace', traces)]:
        a.write_csv(OUT / (name + '.csv'), rows)

    summary = dict(
        date='2026-09-08', scope='Vg=0.8 V, n19/n23 x Vd=0.05/1 V, original 1814 free-carrier-row gates',
        production_defaults_changed=False, new_native_simulations=False,
        constant_direct_DC=len(direct_constants),
        constant_direct_qualified=sum(x['qualified'] == 'True' for x in direct_constants),
        constant_continuation_DC=len(continuation),
        constant_continuation_qualified=sum(x['qualified'] == 'True' for x in continuation),
        constant_strict_dual=sum(x['qualified'] == 'True' for x in constant_dual),
        constant_max_absolute_Id_relative_error=max(abs(float(x['relative_error'])) for x in comparison if x['label'] == 'manual_constants'),
        SRH_preflight_passed=len(preflight), SRH_small_DC_passed=len(small), SRH_response_passed=len(response),
        SRH_max_response_prediction_relative=max(float(x['prediction_relative']) for x in response),
        SRH_initial_finite_qualified=sum(x['qualified'] == 'True' for x in original_finite),
        SRH_initial_finite_total=len(original_finite),
        SRH_restart=a.read(restart.OUT / 'summary.json'),
        SRH_selected_qualified=sum(x['qualified'] == 'True' for x in selected),
        SRH_strict_dual_qualified=sum(x['qualified'] == 'True' for x in dual),
        SRH_joint_qualified=sum(x['qualified'] == 'True' for x in selected if x['axis'] == 'joint'),
        SRH_joint_strict_dual_qualified=sum(x['qualified'] == 'True' for x in dual if x['axis'] == 'joint'),
        SRH_selected_ports_consistent=len(ports),
        native_SRH_fixed_export_replay=True, native_same_source_DC_response_calibrated=False,
        full_field_equivalence_established=False, all_candidate_paths_qualified=False,
        sixteen_conditions_or_full_curves_run=False,
        original_failed_DC_artifacts_retained=2 + 12,
        input_guard_cases=len(input_tests))
    a.write(OUT / 'final_summary.json', summary)

    constant_table = []
    joint_table = []
    field_table = []
    local_table = []
    for case in cases:
        get = lambda label: next(x for x in comparison if x['key'] == case['key'] and x['label'] == label)
        constant_table.append([case['device'], case['vd'], f"{100*float(get('baseline')['relative_error']):+.8f}%",
                               f"{float(get('manual_constants')['relative_error']):+.8e}"])
        joint = next(x for x in currents if x['key'] == case['key'] and x['axis'] == 'joint' and x['label'] == 'finite')
        jd = next(x for x in dual if x['key'] == case['key'] and x['axis'] == 'joint')
        joint_table.append([case['device'], case['vd'], f"{joint['native_relative_error']:+.6e}",
                            f"{max(float(x['max_row_ratio']) for x in selected if x['key']==case['key'] and x['axis']=='joint'):.6e}",
                            jd['qualified']])
        for field in ('phin', 'phip'):
            before = next(x for x in baseline_fields if x['key'] == case['key'] and x['field'] == field)
            after = next(x for x in fields if x['key'] == case['key'] and x['label'] == 'joint/finite' and x['field'] == field)
            field_table.append([case['device'], case['vd'], field,
                                f"{1e3*float(before['max_absolute']):.7g}",
                                f"{1e3*float(after['max_absolute']):.7g}", after['worst_node']])
        if case['vd'] != 1.:
            continue
        for tag, field in [('792', 'phin'), ('792', 'phip'), ('1057', 'phip')]:
            node = case['mapped_nodes'][tag]['node']
            before = next(x for x in baseline_local if x['key'] == case['key'] and x['node_id'] == node)
            after = next(x for x in local if x['key'] == case['key'] and x['label'] == 'joint/finite' and int(x['node_id']) == node)
            local_table.append([case['device'], node, field, f"{1e3*float(before[field+'_delta_V']):+.7f}",
                                f"{1e3*float(after[field+'_delta_V']):+.7f}"])
    # Report generation is downstream of the frozen numerical decisions.
    body = f"""# SimpleMOS 常数约定与局部 SRH 体积的有限自洽验证

日期：2026-09-08。四个 Vg=0.8 V 控制点完成了常数组合的有限自洽对照，最大 Id 相对差降至 **{summary['constant_max_absolute_Id_relative_error']:.8e}**；经固定三步延续，两条初始化路径全部通过原严格门槛。SRH 两节点联合有限修改的 **8/8 端点、4/4 双初始化**通过，局部准费米势明显改善，但全场仍有毫伏级差异。

三个 SRH 候选合计 24 个端点中最终 22 个合格、10/12 双初始化合格；两个单节点分支仍因空穴行未闭合而失败。保留全部初次失败记录，没有降低门槛或修改生产常数、SRH 默认策略。本轮不扩大到 16 工况或完整 0–1 V 曲线。

## 1. 工况、独立性与接受条件

延续[上一轮独立小响应校准](simplemos_constants_srh_response_validation_2026-09-07.md)。n19/n23 对应 NWell Boron 1e17/2e17 cm⁻³，Vg=0.8 V，Vd=0.05/1 V。网格分别 1480/1482 节点、2742/2746 单元；各有 942 个 Si 节点、907 个自由 Si 节点，检查全部 1814 个自由载流子行。几何坐标 μm、状态密度 m⁻³、宽度 1 μm，Id 单位 A/μm。

模型为 300 K、Boltzmann/no-BGN、ni=1.0750038488844236e10 cm⁻³、总杂质 Masetti、掺杂 SRH（τn=1e-5 s、τp=3e-6 s、Nref=1e16 cm⁻³、γ=1）；HFS、表面迁移率、Auger、雪崩、DG 关闭。沿用已验证的自动 Delaunay box、element-box 输运/迁移率、分材料介电系数及三项独立 signed Si Poisson 电荷体积。

Windows MSYS2 UCRT64、C++20 Release，DC 实际使用 Eigen SparseLU/COLAMD；HDF5/TDR、UMFPACK、SPQR 编译可用，本轮未切换 DC 后端。独立切向使用 SciPy SuperLU/COLAMD，加四次 long-double 残差修正。常数程序沿用全部 47 个 core 编译单元统一常数头的隔离构建；SRH 程序只替换其中装配器对象，继续链接其余同约定对象，不与旧生产 core 常数混用。

所有 DC 保持 max_iterations=200、relative_tolerance=1e-7、absolute_tolerance=1e-12、stall floor=1e-9、ψ 更新上限 0.35 V、准费米势上限 0.025 V、contact_basin 初始化、原标量线搜索、四次线性修正。全部载流子行要求 |R|/scale≤1e-6，不排除少数载流子；KCL/Id≤1e-8。全局源闭合继续使用 1e-6 比值及原 1e-10 源下限，低于下限不能宣称相应相对精度。双初始化要求势差≤1e-6 V、密度相对差≤1e-4、Id 相对差≤1e-6，且**两端各自严格合格**。

## 2. 常数有限对照：剩余电流差得到解释

诊断运行同时采用已有 Sentaurus 手册/原生回放给出的 q=1.602192e-19 C、kb=1.380662e-23 J/K、ε₀=8.8542e-12 F/m，未拟合数值。现代 SI 默认仍为 1.602176634e-19、1.380649e-23、8.8541878128e-12。此处只改变常数约定，原 SRH 体积保持不变。

{table(['NWell', 'Vd/V', '现代 SI 基态 Id 误差', '手册常数有限实算 Id 相对差（无量纲）'], constant_table)}

最后一列为 (Id_Vela/Id_native−1)，不是百分数。两个 Vd 下高 NWell 绝对差、低 NWell 控制及高低配对差均改善；说明这四点在已修正几何/介电/迁移率之后的剩余 Id 差主要来自常数约定。四点 ψ 最大绝对差为 {max(float(x['max_absolute']) for x in baseline_fields if x['field']=='psi'):.3g} V，接近双精度和导出精度下限；这不表示器件物理具有该绝对精度，也不能外推到未经计算的偏置。

直接运行 12 次中 10 次合格，原生来源初值的 n19/Vd=0.05、n23/Vd=1 两分支失败。之后从各自原始原生来源初值做冻结 α=0、0.5、1 三步常数延续，12/12 合格，四个终点与直接合格的 Vela 来源路径形成 4/4 严格双初始化。原两次失败及原始状态保留。[直接 DC](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/constants/dc.csv)、[比较](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/constants/comparison.csv)、[延续双初始化](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/constants/continuation/dual.csv)。

原始 `constants/dual.csv` 的 qualified 字段只表示状态距离检查，不能替代两端 DC 合格。最终采用延续路径的严格结果；未回写旧汇总。

## 3. SRH 候选：先独立小扰动，再有限替换

在上述合格手册常数状态上，用几何直接给出的比值定义 V_i(α)=V_original_i×[1+α(V_signedSi_i/V_original_i−1)]。n23 节点 792/1057 对应 n19 的空间匹配节点 791/1056；位置分别约 (0.1707979, 0.3)、(0.005045938, 0.1127039) μm，最大跨网格匹配距离 9.23e-10 μm，小于冻结的 1e-8 μm。

第一节点体积比为 1.17958462195；第二节点 n19/n23 分别为 0.217276367556/0.217274908657。执行第一节点、第二节点和联合三支路，α=±0.001、±0.0005 后才到 α=1。只改变 SRH 源积分及对应六项状态导数；三项 Poisson 体积、输运几何、固定状态通量均保持原定义。[空间与体积表](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/spatial_scope.csv)。

12/12 预检查通过：零扰动残差逐位一致，独立源公式与固定状态残差差分最大相对差 {max(float(x['source_relative']) for x in preflight):.6g}（门槛 1e-5）；完整线性切向残差最大 {max(float(x['linear_relative']) for x in preflight):.6g}（1e-8）；弱交叉块与预期体积缩放最大相对差 2.22e-16（1e-12）。同时通过固定输运不变性、强块 Jv 与端口一致性。

52/52 小扰动 DC、24/24 响应检查通过。主响应预先冻结为全部 907 个自由 Si 节点的 (ψ,φn,φp) 向量，独立完整 J 解的预测最大相对差 {summary['SRH_max_response_prediction_relative']:.6g}，双幅度比例差最大 {max(float(x['two_amplitude_relative']) for x in response):.6g}，均低于 1e-3；最大偶/奇比 {max(float(x['even_over_odd']) for x in response):.6g}，最小响应/零漂移比 {min(float(x['signal_over_zero_drift']) for x in response):.6g}。未以不敏感的 Id 代替少数载流子场校准。[预检查](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/preflight.csv)、[响应](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/response.csv)。

## 4. 有限 SRH 结果与资格

最初 24 次有限 DC 的低 Vd 12 次通过，高 Vd 12 次失败。冻结一次重启，逐支使用原失败状态继续同一有限参数，数值配置与门槛不变，仅增加只读更新诊断。12 次重启中 10 次通过，合并初次合格状态后为 22/24 严格端点、10/12 严格双初始化；全部 24 个选定端点的装配器、端口提取器和探针电流一致，但端口一致不抵销 DC 失败。

联合支路两种初始化的 8/8 端点和 4/4 双初始化均通过：

{table(['NWell', 'Vd/V', '联合后 Id 相对差（无量纲）', '两初值最大载流子行比值', '严格双初始化'], joint_table)}

仅就已修改位置，高 Vd 的有符号准费米势差如下。基线为手册常数、原 SRH 体积；联合结果取选定的合格 Vela 来源支路，另一初值已通过上述严格检查。

{table(['NWell', '本地节点', '场', '基线差/mV', '联合后差/mV'], local_table)}

n23 节点 792 电子密度相对差从约 −13.674% 降至 +0.3941%，1057 空穴密度相对差从 +16.289% 降至 +4.1860%。联合后这两点局部 SRH 体积率与原生的相对差分别约 +2.30e-10、−4.53e-9；体积率一致不意味着积分源与全场已全局一致。n23 节点 1000/1009 的 φp 差仍约 19.239/20.802 μV，变化很小。

全部 907 个自由 Si 节点的最大准费米势绝对差仍为：

{table(['NWell', 'Vd/V', '场', '基线最大差/mV', '联合后最大差/mV', '联合后最差节点'], field_table)}

因此，局部两体积能解释被选热点的大部分偏差，尚未消除邻近节点或其他 SRH 支撑范围的差异。联合状态高 Vd 的 n23 最大电子/空穴密度相对差仍约 11.84%/6.79%；只看电势和端口电流会遗漏这些少数载流子误差。[最终选定状态](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/restart/selected_states.csv)、[严格双初始化](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/restart/dual.csv)、[局部场](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/restart/local_fields.csv)、[全场误差](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/restart/fields.csv)。

初次 `srh/dual.csv` 同样只是状态距离指标；最终严格判定使用 `srh/restart/dual.csv`，失败态在场表中带 qualified_state=False，未混入上面的合格联合结果。

## 5. 原生固定状态 SRH 行回放与剩余失败

另用已导出的原生电势、准费米势、密度、SRH 体积率、box 几何及原生单元加权 Masetti 独立重建局部 SG 流出项，分别使用准费米势和密度两种公式，共保留 32 个行记录。下述源归一化指标统一使用 max(|净流出|, |原生 SRH×signed Si 体积|) 为分母。这是原生固定导出状态回放，**不是原生同源正负自洽响应校准**，未启动新 sdevice。

n23/Vd=1 的节点 792 电子/空穴行用 signed Si 体积后，源归一化残差约 2.40e-15/8.58e-16；沿用原 Vela SRH 体积则约 0.152244。节点 1057 的空穴行分别约 2.42e-12 与 3.60246，低 NWell 对应少数载流子行也支持同一几何解释。节点 1057 的多数电子行有严重大通量消减：绝对边通量和约 4.68e12 particles/(cm·s)，而原生积分源约 −0.112 particles/(cm·s)，准费米势导出回放净通量约 −4.83。该行无法在源尺度上反推体积；其异常/负推断体积保留在表中，不当作物理体积结论。[32 行回放](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/native_rows/row_replay.csv)。

最终仍失败的是 Vd=1、仅改第二节点、原生来源初始化的两个分支：

{table(['NWell', '空穴节点', '最终行比值（门槛 1e-6）', '原始=截断更新/V', '最后 Poisson 块残差'], [[bykey[x['key']]['device'], x['node_id'], f"{float(x['final_row_ratio']):.9g}", f"{float(x['raw_linear_step_V']):+.6e}", f"{float(x['final_Poisson_block']):.6e}"] for x in traces])}

两点最后均为 carrier_row_convergence_line_search_rejected，13 次尝试未接受；原始更新没有触及 0.025 V 上限。对应行 raw J·dx+R 约 −6.31e-41/+1.10e-39（求解器缩放单位），最优但被拒候选仍稍微降低该空穴行残差。diagnostic 的 applied_step_V 在 best_rejected_candidate=1 时描述候选更新，并未提交到最终状态。这表明更新截断不是这里的限制；下一步应独立核对接近 Poisson 求值精度下限时的全局 merit/线搜索决定。单行改善及本次双精度线性诊断仍不足以证明任意完整步可接受。[剩余失败](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/restart/remaining_violations.csv)、[原始/截断/候选更新](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/remaining_update_trace.csv)。

## 6. 完成范围与后续边界

本轮有限常数对照、SRH 分项/联合对照、双初始化、原生固定状态源行审计已执行。九项诊断输入拒绝测试通过（网格数、缺项、接触节点、重复、非正比值、多余记录、幅度越界、缺少幅度、旧/新接口冲突），各阶段哈希复核通过。未更改生产默认，未重跑生产全量 CTest；上一轮生产测试状况仍见[几何生产报告](simplemos_generated_box_mobility_validation_2026-09-07.md)，本报告不将其重新计为今日测试。

下一阶段先对上述两个失败态做固定原始方向的 trial residual/merit 独立精度扫描，区分 Poisson 消减与少数载流子真实下降，不改接受阈值；并在合格联合状态上审计 n23 的 794/1056、n19 的 793/794 及低 Vd 1086/1087 所在邻域的 SRH 支撑。只有邻域源定义与同扰动校准支持后，再冻结更完整的 SRH 体积候选；原生同源自洽响应资格尚未取得。常数兼容策略如进入生产，应做显式配置并保留 SI 默认。Vg=1、16 工况和全曲线仍需另行逐步验证。

脚本：[常数有限对照](../../scripts/validate_simplemos_constants_finite_20260908.py)、[常数延续](../../scripts/continue_simplemos_constants_finite_20260908.py)、[SRH 隔离构建](../../scripts/build_simplemos_srh_finite_20260908.py)、[SRH 校准与有限对照](../../scripts/validate_simplemos_srh_finite_20260908.py)、[一次重启](../../scripts/restart_simplemos_srh_finite_20260908.py)、[原生行审计](../../scripts/audit_simplemos_native_srh_rows_20260908.py)、[输入校验](../../scripts/test_simplemos_srh_finite_inputs_20260908.py)、[报告生成](../../scripts/report_simplemos_constants_srh_finite_20260908.py)。[最终汇总](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/final_summary.json)与[最终证据](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/final_evidence.json)关联全部冻结合同、原失败记录及本报告。
"""
    assert not REPORT.exists()
    REPORT.write_text(body, encoding='utf-8')
    scripts = [p.REPO / 'scripts' / (stem + '_20260908.py') for stem in (
        'validate_simplemos_constants_finite', 'continue_simplemos_constants_finite',
        'build_simplemos_srh_finite', 'validate_simplemos_srh_finite',
        'restart_simplemos_srh_finite', 'audit_simplemos_native_srh_rows',
        'test_simplemos_srh_finite_inputs', 'report_simplemos_constants_srh_finite')]
    for path in scripts:
        ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
    links = re.findall(r'\]\(([^)]+)\)', body)
    for link in links:
        if link.endswith('/final_evidence.json'):
            continue
        assert (REPORT.parent / link).resolve().exists(), link
    a.write(OUT / 'report_checks.json', dict(script_syntax_passed=len(scripts), local_links_checked=len(links),
                                            stage_hash_manifests_verified=len(manifests),
                                            strictly_qualified_joint_endpoints=8, strictly_qualified_joint_pairs=4,
                                            full_CTest_rerun=False))
    d.matrix.freeze(OUT / 'final_evidence.json', manifests + scripts + [REPORT] + list(OUT.glob('*.csv')) +
                    [OUT / 'final_summary.json', OUT / 'report_checks.json'])
    a.verify(OUT / 'final_evidence.json')
    print(summary, flush=True)
    print('Report and final evidence:', a.rel(REPORT), flush=True)


if __name__ == '__main__':
    main()
