"""Seal production migration, unchanged-cohort results and roundtrip limitation."""
import ast
import re
from pathlib import Path
import validate_simplemos_production_weighted_20260910 as v


def main():
    a,d=v.a,v.d;out=v.OUT
    manifests=[out/'comparison_evidence.json',out/'jvp_evidence.json',out/'source/production_evidence.json',
               out/'failed_restart_isolation/evidence.json',out/'failed_roundtrip/evidence.json']
    for p in manifests:a.verify(p)
    dc=a.read(out/'comparison_summary.json');jvp=a.read(out/'jvp_summary.json');source=a.read(out/'source/production_summary.json')
    trip=a.read(out/'failed_roundtrip/summary.json');points=a.rows(out/'comparison.csv')
    failed=[r for r in points if r['comparison_qualified']!='True'];assert len(failed)==1
    f=failed[0];root=v.LOCAL.parent
    tests=(root/'ctest_weighted.log').read_text(encoding='utf-8-sig')
    assert '16 tests failed out of 795' in tests
    def failures(text):return set(re.findall(r'\d+ - (\S+) \(Failed\)',text))
    before=failures((root/'ctest.log').read_text(encoding='utf-8-sig'));after=failures(tests)
    assert len(after)==16 and before==after
    script_names=['migrate_simplemos_numerics_20260910','validate_simplemos_production_precision_20260910',
        'validate_simplemos_production_unit_block_20260910','validate_simplemos_production_weighted_20260910',
        'validate_simplemos_production_source_20260910','validate_simplemos_production_weighted_source_20260910',
        'isolate_simplemos_production_failed_restart_20260910','finish_simplemos_failed_restart_isolation_20260910',
        'check_simplemos_failed_state_roundtrip_20260910','report_simplemos_production_precision_20260910']
    scripts=[v.REPO/'scripts'/(n+'.py') for n in script_names]
    for p in scripts:ast.parse(p.read_text(encoding='utf-8'),filename=str(p))
    table=[]
    for vg in (0.,.2,.8,1.):
        rr=[next(r for r in points if r['device']=='n23' and float(r['vg'])==vg and float(r['vd'])==vd) for vd in (.05,1.)]
        table.append(f"| {vg:g} | {float(rr[0]['vela_error_percent']):+.8f}% | {float(rr[1]['vela_error_percent']):+.8f}% | {'两点均通过' if all(r['comparison_qualified']=='True' for r in rr) else '高 Vd 双初始化未通过'} |")
    report=v.REPO/'docs/validation/simplemos_production_precision_and_sixteen_point_2026-09-10.md'
    prefix='../../reference_tcad/simplemos_sentaurus2022/phumob_numerics_production_20260910/weighted_merit/'
    text=f'''# SimpleMOS 数值修复生产移植与十六点复验

日期：2026-09-10。分支：`codex/simplemos-sdevice-validation`。

## 结论

低 Vd 局部源资格缺口已闭合，正式实现复现了合格诊断结果：20/20 DC、8/8 全场响应、16/16 原生端口响应通过，20 份状态及端口电流与诊断实现逐位一致。但使用原独立种子、同一实现和原重载协议的 16 点双初始化仍为 **15/16**，不能放行完整曲线及 Enormal/HFS。

未通过点由此前的 n23、Vd=1、Vg=0.2 转为 **n23、Vd=1、Vg=0** 的原生初始化。此次没有将新旧两版的合格点合并成 16/16。额外重载虽能通过，但不计入正式协议；零迭代诊断确认保存态往返存在残差表示差。

前置独立精度证据见 [低 Vd 源扰动精度报告](simplemos_lowvd_source_precision_2026-09-10.md)。本轮没有新运行 Sentaurus，复用了已授权并取回的同源原生结果；未修改常数或 SRH 体积，未提交或推送。

## 配置与生产修改

沿用 n19/n23 原匹配网格、原生 box 修改及 Si 输运几何，n19 为 1480 节点、n23 为 1482 节点；采用 signed Si 的三项 Poisson 电荷体积及逐单元介电系数。物理为 300 K、Boltzmann、OldSlotboom 与匹配基准 ni、plain PhuMob（`element_box_phumob`）、掺杂相关 SRH，Enormal/HFS 关闭。电流单位 A/μm，势单位 V，Vela 密度单位 m⁻³。实际 DC 使用 Eigen SparseLU 和四次线性迭代修正；UCRT64 Release 构建中的 UMFPACK/SPQR 可用性不代表此次使用它们。

新增三个显式选项，默认均保持原路径：

- `poisson_residual_precision: binary128`：由实际 packed 坐标与准费米参考值计算 Boltzmann 电荷和 Poisson 残差，最后舍入至 double。状态存储、原几何/材料系数、连续性与 Jacobian 不因该选项改变；排除 Fermi、量子、热发射接触和反馈替代分支。
- `stable_merit_comparison: true`：保持原 L2 或块范数的权重及自动归一化系数，用补偿求和判断平方范数差；不确定时使用精确整数或有理数符号。仅支持原 `merit` 线搜索且全局闭合 merit 关闭。只读残差/闭合探针关闭无用的线搜索选项，仍启用相同 Poisson 求值。
- `exact_dirichlet_updates: true`：在线性修正后令显式接触 identity 行的更新严格等于负残差，再进入原截断和线搜索。自由行更新不变；不支持非零载流子正则化。

另修正准费米保存态重启的运算顺序：先相减旧/新参考值，再相加存储增量，保留小于参考值舍入间隔的增量。此修复不改 CSV 格式；**不代表所有 packed 坐标的保存/加载一致性已解决**。接口说明见 [配置文档](../config_schema.md)。

## 原接受协议下的十六点结果

[输入及门槛]({prefix}validation_contract.json)、[全部尝试]({prefix}attempts.csv)、[双初始化比较]({prefix}comparison.csv)、[汇总]({prefix}comparison_summary.json)。

共 32 次首次尝试，6 次原协议允许的同偏压重载，总计 {dc['attempts']} 次，其中 {dc['failed_attempts']} 次失败完整保留。最终选取状态 {dc['qualified_selected']}/32 合格，双初始化 {dc['qualified']}/16 合格。门槛仍为全部自由 Si 载流子逐行 1e-6、KCL/Id 1e-8、端口一致性 1e-8，双初始化势差 1e-6 V、密度相对差 1e-4、Id 相对差 1e-6。

未通过点的原生初始化重载后，空穴节点 800 最大行比 {dc['max_row_ratio']:.9g}，仍有 91 个违规行。两条初始化路径的最大空穴准费米势差 {float(f['phip_max_V']):.9g} V，最大密度相对差 {float(f['density_max_relative']):.9g}，均超门槛；Id 相对差仅 {float(f['dual_Id_relative']):.6g}，不能替代行与场资格。

n23 的有符号电流误差为 `100*(Id_Vela/Id_Sentaurus-1)`，下表取合格 Vela 初始化；高 Vd、Vg=0 只具有单路径资格：

| Vg / V | Vd=0.05 V | Vd=1 V | 双初始化 |
|---|---:|---:|---|
{chr(10).join(table)}

全部 16 个合格 Vela 初始化状态的 [Jv 分块检查]({prefix}jvp_blocks.csv) 中，{jvp['checks']}/{jvp['checks']} 项正式检查通过，最大相对差 {jvp['max_relative']:.9g}，门槛 1e-4。弱 SRH 交叉块继续单独标记，不能用全残差差分覆盖其资格。64 个全局 SRH 分量均未激活旧源下限，未增加独立相对源闭合资格。

## 同源响应与生产一致性

[源预检]({prefix}source/preflight.csv)、[生产源结果]({prefix}source/production_summary.json)、[与诊断实现逐位比较]({prefix}source/diagnostic_identity.csv)、[原生端口对照]({prefix}source/native_ports.csv)。

原四个 Vg=0 工况、原八节点固定成对源、零源与 ±0.001/±0.0005 两档扰动不变。新编译的诊断覆盖层只加入固定源和只读完整 J 导出，Newton 及 Poisson 使用正式库。4/4 新鲜完整 J/源插入预检通过；20/20 DC、8/8 全场响应、16/16 原生漏/衬底端口响应通过。最大端口导数相对差 {source['max_native_relative']:.9g}，即 {100*source['max_native_relative']:.9g}%，低于原 0.1% 门槛。

20 份输出状态和端口电流与合格诊断实现逐位相同，证明此批次的移植一致性。通用旧汇总保留的 `native_response_qualified:false` 是历史未定义的通用资格字段，本轮具体端口资格见 `qualified_native_ports=16`；没有据此放行有限 SRH 体积替换或全部独立初始态。

## 剩余保存态往返问题

[四项同态隔离]({prefix}failed_restart_isolation/results.csv)、[零迭代往返]({prefix}failed_roundtrip/roundtrip.csv)、[范数变化]({prefix}failed_roundtrip/summary.json)。

冻结正式未通过态后，另外重载一次的 combined、关闭稳定 merit、关闭接触投影三项均在 1 次 Newton 更新后通过，最大行比 1.535602706e-7；只将 Poisson 改回 double 的隔离在 7 次更新后通过，最大行比 1.945244646e-10。这些额外诊断全部不计入正式双初始化资格，也不能证明关闭某项即可从原独立种子通过。

进一步只做保存—加载、**零次 Newton 更新**：ψ、两种准费米势、密度的导出值全部不变；但原结束残差 {trip['prior_final_residual']:.12g} 变为重载初始残差 {trip['first_reload_initial_residual']:.12g}，为原来的 {trip['norm_ratio']:.9g} 倍。第二次零迭代重载稳定。这说明导出的物理场不足以复现该状态的内部残差表示，值得核对 packed ψ/φn/φp、参考值与缩放往返；本轮尚未把差异唯一分配至某一种内部坐标。

## 验证与保留记录

Release 全目标构建完成。数值测试 12 个用例、1222 项断言通过，微小准费米增量重启测试 1 个用例、6 项断言通过。完整 CTest **779/795**，16 项失败与本轮前一次 778/794 的历史文件身份/缺失文件失败集合完全相同，无新增失败；原始输出见 [CTest 日志](../../build-release/phumob_numerics_production_20260910/ctest_weighted.log)。{len(scripts)} 个本轮生产验证脚本完成语法检查。

初始生产保护只接受显式 L2，32 个输入被配置拒绝，未进入 DC；随后单位块保护在非单位自动尺度初始态被拒绝，并发现只读探针不应启用线搜索保护。这两批原始记录及源/库快照保留，最终独立批次支持原块权重，不改初始态或门槛。新增测试曾缺 JSON 头文件，已修复；一个隔离的日志文件超过 Windows 路径长度，保留原输出并在短路径仅重跑该控制。

## 后续门槛

优先独立审计并修复保存/加载的内部坐标与残差往返一致性，再从同一 32 个原独立种子复验，保留一次重载上限及本轮失败记录。该任务通过后才能宣称统一实现 16/16；随后运行四条完整 PhuMob 0–1 V 曲线，再开展 Enormal 和高场饱和的逐项本构/导数/自洽校准。本轮未启动这些曲线或新物理模型仿真。
'''
    # Reports are new immutable artifacts; do not overwrite a prior conclusion.
    assert not report.exists();report.write_text(text,encoding='utf-8')
    for target in re.findall(r'\]\(([^)]+)\)',text):assert (report.parent/target).resolve().exists(),target
    summary=dict(production_source=source,cohort=dc,jvp=jvp,roundtrip=trip,ctest_total=795,ctest_passed=779,
        historical_failures=sorted(after),new_failures=0,scripts_parsed=len(scripts),full_curves_launched=False,
        Enormal_HFS_launched=False,native_new_runs=0,report=str(report.relative_to(v.REPO)))
    a.write(out/'completion_summary.json',summary)
    d.matrix.freeze(out/'completion_evidence.json',manifests+[out/'completion_summary.json',report,
        root/'ctest_weighted.log',root/'numerics_test_weighted.log',root/'warm_restart_test.log',root/'build_weighted.log']+scripts)
    print(summary,flush=True)

if __name__=='__main__':main()
