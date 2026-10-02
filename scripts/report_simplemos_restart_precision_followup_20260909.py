"""Summarize both precision axes and preserve every failed experiment."""
import ast
import math
import re
from pathlib import Path
import validate_simplemos_restart_source_followup_20260909 as w
import audit_simplemos_e1089_linesearch_20260909 as v
import validate_simplemos_e1089_stable_merit_20260909 as m
s=w.s;OUT=w.OUT
REPORT=s.REPO/'docs/validation/simplemos_restart_and_e1089_precision_followup_2026-09-09.md'
STATUS=s.REPO/'docs/validation/simplemos_branch_status.md'

def link(path,title):return f'[{title}](../../{path.relative_to(s.REPO).as_posix()})'

def run():
    assert not REPORT.exists()
    manifests=[OUT/'freeze.json',OUT/'identity_evidence.json',OUT/'preflight_evidence.json',OUT/'dc_evidence.json',OUT/'evidence.json',OUT/'zero_plateau/evidence.json',OUT/'settled_source/freeze.json',OUT/'settled_source/preflight_evidence.json',OUT/'settled_source/dc_evidence.json',OUT/'settled_source/evidence.json',v.OUT/'build_evidence.json',v.OUT/'freeze.json',v.OUT/'replay_evidence.json',v.OUT/'precision/freeze.json',v.OUT/'precision/evidence.json',m.OUT/'build_evidence.json',m.OUT/'freeze.json',m.OUT/'dc_evidence.json']
    manifests.append(m.OUT/'disabled_identity/evidence.json')
    for path in manifests:s.a.verify(path)
    original=s.a.rows(OUT/'response.csv');dc=s.a.rows(OUT/'dc.csv');settled=s.a.rows(OUT/'settled_source/dc.csv');stable=s.a.rows(m.OUT/'dc.csv');precision=s.a.rows(v.OUT/'precision/merit.csv');checks=s.a.rows(v.OUT/'precision/checks.csv');flux=s.a.rows(v.OUT/'precision/independent_flux_checks.csv')
    repaired=next(r for r in stable if r['key']=='failed_e1089');control=next(r for r in stable if r['key']=='qualified_control');source=[r for r in stable if r['kind']=='source']
    repaired_status=s.a.read(Path(repaired['dest'])/'config.status.json');before=s.a.read(s.LOCAL/'restart_coordinates/dc/native_reload/config.status.json')
    ratio= float(repaired['max_row_ratio']);passed=repaired['qualified']=='True'
    summary=dict(original_source_DC_qualified=sum(r['qualified']=='True' for r in dc),original_source_DC=len(dc),original_field_responses_qualified=sum(r['qualified']=='True' for r in original),original_field_responses=len(original),zero_plateau_restarts=4,zero_plateau_exactly_unchanged=all(float(r['drift_V'])==0 for r in s.a.rows(OUT/'zero_plateau/results.csv')),settled_DC_qualified=sum(r['qualified']=='True' for r in settled),settled_DC=len(settled),stable_source_DC_qualified=sum(r['qualified']=='True' for r in source),stable_source_DC=len(source),e1089_repair_qualified=passed,e1089_repair_max_row=ratio,e1089_repair_iterations=int(repaired['iterations']),e1089_repair_failure=repaired['failure'],e1089_original_max_row=before['carrier_row_convergence']['max_ratio'],e1089_original_KCL_over_Id=1.901076798419972e-9,e1089_repair_KCL_over_Id=float(repaired['kcl_over_Id']),precision_snapshots=len(checks),actual_rejected_but_exactly_descending=sum(r['would_decrease']=='True' and r['label']!='base' for r in precision if r['mode']=='production'),independent_BGNSG_checks=len(flux),maximum_precision_agreement=max(float(r['precision_agreement']) for r in checks),maximum_rounded_reconstruction=max(float(r['reconstruction']) for r in checks),tests='3 Catch2 cases, 57 assertions; 13 real trial signs include 8 exact decreases',production_changed=False,threshold_changed=False,finite_volume_replacement=False)
    s.a.write(OUT/'final_summary.json',summary)
    table=[]
    for key in sorted({r['key'] for r in original}):
        r=next(x for x in original if x['key']==key and x['amplitude']=='full');h=next(x for x in original if x['key']==key and x['amplitude']=='half');table.append(f"| {key} | {float(r['signal_over_drift']):.6g} | {float(h['signal_over_drift']):.6g} | {sum(x['qualified']=='True' for x in original if x['key']==key)}/2 |")
    outcome='通过全部当前 DC 资格' if passed else '仍未通过全部当前 DC 资格'
    body=f'''# SimpleMOS 重启漂移与电子 1089 收敛精度跟进

日期：2026-09-09。工作分支 `codex/simplemos-sdevice-validation`。接续 [原生局部成对源校准](simplemos_phumob_native_pair_source_validation_2026-09-09.md)。

**完成两个独立数值轴的定位及隔离验证。电子 1089 的精确范数差隔离重算{outcome}，最终最大逐行比 {ratio:.9e}；低 Vd 全场资格仍未放行。** 仅重启修复的原基态实验仍为 6/8 全场通过；稳定基态的零扰动可逐位重复，但正负扰动暴露线搜索问题。原失败全部保留。

## 配置和资格

Vg=0 的 n19/n23 × Vd=0.05/1 V 用于局部源；n23、Vd=1、Vg=0.2 用于电子 1089 失败态及同偏置合格控制。网格分别 1480/1482 节点，各 907 个自由 Si 节点，1814 个受验收载流子行。300 K、Boltzmann、OldSlotboom、plain PhuMob、匹配 ni、掺杂相关 SRH；保持已验证的几何、介电系数及 `element_box_phumob`，SRH 仍用原 all-cell 体积。

Windows MSYS2 UCRT64 Release C++20，实际 Eigen SparseLU，四次线性迭代修正。独立切向为 SciPy SuperLU 加四次 long-double 缺陷修正。沿用原逐行 1e-6、KCL/Id 1e-8、max_iter=200、abstol=1e-12、reltol=1e-7 及原残差停滞条件；全场双幅度、预测相对门槛 1e-3，偶/奇比 0.01，信号/零漂移至少 100。未修改接受门槛。全局弱源沿用原下限，不新增独立弱源相对资格。

本轮没有新增 Sentaurus 仿真，复用此前原生 20/20 DC、8/8 端口响应校准证据；物理常数、物理源定义、α=±0.001/±0.0005 和生产 C++ 均未改变。所有代码变化为新脚本及 ignored 构建目录中的隔离程序。

## 低 Vd 漂移：重启坐标与稳定基态分开检查

把上轮已验证的“先减参考值、再加增量” Newton 对象与原校准源装配器组合。第一次全部节点身份检查失败；审计证明变化仅为 Dirichlet 接触上的微小残留投影到零，全部自由 Si 状态和载流子残差逐位一致。没有扩大物理验收域或放宽阈值，原预检失败和逐项接触差另存。{link(OUT/'identity_domain_amendment.json','检查域说明')}、{link(OUT/'contact_projection.csv','接触差异')}。

原基态 20/20 DC 通过，场响应仍为 6/8：

| 工况 | 全幅信号/漂移 | 半幅信号/漂移 | 全场资格 |
| --- | --- | --- | --- |
{chr(10).join(table)}

n23 低 Vd 的零源漂移主要为电势：自由 Si ψ 向量差的 L2 范数 5.21825e-14 V，最大节点差 2.45e-15 V；φn/φp 向量差约 1.10940e-14/8.57958e-15 V。重启增量修复并未显著改变原信号/漂移。{link(OUT/'zero_drift_components.csv','漂移分量')}、{link(OUT/'response.csv','原基态响应')}。

从四个合格零源结果中，分别对两个低 Vd 控制继续做两次零源重启，共 4/4 严格 DC 通过，均 0 次迭代、自由 Si 场差为零。该结果说明第一次零源运行完成了基态微调，不能以此回写原实验的漂移。随后独立冻结稳定基态，保持原物理源值不变，重算完整 J、源插入和切向预检，两个控制均通过。

然而稳定基态的 10 次 DC 只有 {sum(r['qualified']=='True' for r in settled)}/10 通过（两个零控制）；8 个正负扰动均保留为失败，全场 0/4。仅启用精确范数差的隔离组也只有 {sum(r['qualified']=='True' for r in source)}/10 DC 通过。零漂移不代表扰动已收敛，失败状态的差分不能被视为有效响应；表中的历史 6/8 资格仍保持原义。{link(OUT/'zero_plateau/results.csv','两次零重启')}、{link(OUT/'settled_source/preflight.csv','稳定态切向预检')}、{link(OUT/'settled_source/dc.csv','稳定态失败批次')}、{link(m.OUT/'dc.csv','精确比较隔离批次')}。

## 电子 1089：实际回溯与独立高精度残差

从上轮 52 次重算后的最终失败态执行一次实际 Newton 步，保留生产行缩放、四次线性修正、更新截断与回溯。另对同偏置合格控制做对照。输出开关两组均逐位保持最终状态、退出状态和重复残差；失败态导出基态及全部 13 个回溯候选，合格控制导出基态及其实际候选。

独立 Decimal 60/100 位复核包含全部 Poisson 和两类自由 Si 载流子行，区分：生产舍入通量/源项的高精度求和、原双精度内核操作数、分离状态的精确物理坐标、未舍入的数学试探状态。BGN-SG 明确保留 ni 梯度漂移，并与单独表达的密度 SG 验证。{len(flux)} 个独立 SG 检查通过；14 个快照的 60/100 位最大归一化差 {summary['maximum_precision_agreement']:.3e}，生产舍入项重建最大差 {summary['maximum_rounded_reconstruction']:.3e}。这些误差按通量/源项规模归一化；不等价于小残差本身的相对误差。

电子 1089 基态逐行比在生产、内核高精度、分离状态下均约 1.04590e-6，确实超过 1e-6。原始和截断步完全相同，完整 J 的独立线性缺陷在该行仅约原残差的 4.90e-18；不能归因为截断或线性求解缺陷。

全步把该行降到 4.80e-12，但生产总范数增加约 7.79%，Poisson 块阻止接受。α=1/32 及更小的 8 个候选，生产残差平方和实际下降约 2.44e-36～1.93e-38，却被原双精度范数比较判为不下降。α=1/32 时，1443 个非零预期 ψ 更新全部被双精度状态加法吸收；此时载流子仍有可解析的下降。高精度“分离/数学状态”与生产状态必须分开解释，不能据其下降强行接收生产候选。{link(v.OUT/'precision/merit.csv','13次回溯与范数差')}、{link(v.OUT/'precision/coordinates.csv','状态更新吸收')}、{link(v.OUT/'precision/linear.csv','线性缺陷')}。

## 精确范数差的隔离自洽结果

仅将原未加权 L2 的严格下降判断替换为平方和差的可靠符号比较：先用 long-double 补偿求和及误差界，必要时用精确二进制整数平方和。继续保留原零到零例外、正载流子检查和回溯网格。测试包括大小量级混合、次正规数、非有限值及此次 13 个真实回溯，3 个 Catch2 用例、57 项断言通过。关闭该选择器时，电子 1089 和 n19 低 Vd 带源对照的已完成轨迹均逐位复现。

| 项目 | 原失败态 | 精确比较隔离重算 |
| --- | --- | --- |
| 最大载流子行比 | 1.045902376e-6 | {ratio:.9e} |
| 行超标数 | 1 | {repaired['row_violations']} |
| KCL/Id | 1.901076798e-9 | {float(repaired['kcl_over_Id']):.9e} |
| 接受迭代数 | 本次重载 0 | {repaired['iterations']} |
| 当前 DC 资格 | 失败 | {'通过' if passed else '失败'} |
| Id (A/μm) | {before['current_total_A_per_um']:.12e} | {float(repaired['current_A_per_um']):.12e} |

本轮最终退出原因：`{repaired_status.get('convergence_reason') or repaired_status.get('failure_reason')}`。同偏置合格控制 {'通过' if control['qualified']=='True' else '失败'}，最大行比 {float(control['max_row_ratio']):.3e}。低 Vd 的 8 个非零扰动最终失败态在精确/普通比较之间逐位相同，支持“仅改比较不足以解除该阻塞”的判断；n19 零控制多走了两个微小下降步，状态不计逐位相同。该单点结果没有覆盖完整双初始化和 16 工况资格，也不能替代低 Vd 失败的小扰动验证。不能将深关断电流差异归因为本次发现的范数比较问题。

## 后续边界与证据

下一步应针对稳定低 Vd 基态的真实带源回溯，独立检查 Poisson 残差求值及小电势更新的表示精度；使用同源、同幅度、n19/n23 配对和既有严格门槛做隔离对照。精确范数差的生产移植应在补齐相关回归及控制范围后另行开展。当前不进行有限 SRH 体积替换、Enormal/HFS 或完整曲线放行。n23、Vg=0、Vd=1 的现有电流误差仍为 −0.80813524%。

保留全部失败，包括第一次接触身份检查、稳定基态正负扰动，以及测试构建首次漏传 Eigen 的 include 路径（修正测试命令后通过；求解程序编译/链接已成功）。本轮未修改生产文件，未跑新的全量 CTest，未提交或推送。{link(OUT/'final_summary.json','机器汇总')}、{link(m.LOCAL/'test_result.log','本轮数值测试')}。
'''
    REPORT.write_text(body,encoding='utf-8')
    old=STATUS.read_text(encoding='utf-8');start=old.index('## 当前原生局部源校准与剩余资格');end=old.index('## 前置 PhuMob 低栅压控制',start)
    entry=f'''## 当前重启漂移与电子 1089 验证

最新报告为 [重启漂移与电子 1089 精度跟进]({REPORT.name})。重启坐标修复后，原四点局部源仍为 20/20 DC、6/8 全场通过；低 Vd 首次零源微调后，连续两次重启场差为零，但新稳定基态正负扰动失败，不能将零漂移视作全场放行。稳定基态普通比较及精确比较的 DC 分别为 {sum(r['qualified']=='True' for r in settled)}/10、{sum(r['qualified']=='True' for r in source)}/10。

实际电子 1089 回溯确认 8 个下降候选被普通双精度范数比较漏判。Decimal 60/100 位独立 BGN-SG/SRH/Poisson 检查通过，精确范数差的 3 个 Catch2 用例、57 项断言通过。隔离重算{outcome}，最终最大行比 {ratio:.9e}，迭代 {repaired['iterations']} 次；完整双初始化与 16 工况资格未因此放行。

前置 [原生同源校准](simplemos_phumob_native_pair_source_validation_2026-09-09.md) 20/20 DC、8/8 端口响应及最大原始响应差 0.001002% 保留。本轮实际 SparseLU、四次线性修正；未新增原生仿真、未改生产 C++、常数、门槛或 SRH 体积。下一步检查低 Vd 带源候选的 Poisson 残差及电势表示精度，再决定精确比较的生产移植和有限替换。深关断 Id 误差仍为 −0.80813524%。本轮未提交或推送。

'''
    STATUS.write_text(old[:start]+entry+old[end:],encoding='utf-8')
    scripts=[s.REPO/'scripts'/name for name in ('validate_simplemos_restart_source_followup_20260909.py','qualify_simplemos_restart_source_domain_20260909.py','check_simplemos_lowvd_zero_plateau_20260909.py','validate_simplemos_settled_lowvd_source_20260909.py','audit_simplemos_e1089_linesearch_20260909.py','analyze_simplemos_e1089_precision_20260909.py','validate_simplemos_e1089_stable_merit_20260909.py','report_simplemos_restart_precision_followup_20260909.py')]
    scripts.append(s.REPO/'scripts/check_simplemos_stable_merit_disabled_20260909.py')
    for path in scripts:ast.parse(path.read_text(encoding='utf-8'))
    links=re.findall(r'\]\(([^)]+)\)',body)
    for target in links:assert (REPORT.parent/target).resolve().exists(),target
    s.a.write(OUT/'report_qa.json',dict(scripts_syntax_checked=len(scripts),links_checked=len(links),manifests_verified=len(manifests),direct_hash_entries_verified=sum(len(s.a.read(p)['input_hashes']) for p in manifests)))
    s.d.matrix.freeze(OUT/'completion_evidence.json',manifests+scripts+[OUT/'zero_drift_components.csv',m.OUT/'source_identity.csv',OUT/'final_summary.json',OUT/'report_qa.json',REPORT,STATUS]);s.a.verify(OUT/'completion_evidence.json');print(summary,flush=True)

if __name__=='__main__':run()
