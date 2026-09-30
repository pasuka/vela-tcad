"""Seal local source calibration and restart diagnostics, including all failures."""
from pathlib import Path
import math
import py_compile
import calibrate_simplemos_phumob_local_source_20260909 as s

REPORT=s.REPO/'docs/validation/simplemos_phumob_local_source_and_restart_validation_2026-09-09.md'

def link(path,label):return f'[{label}](../../{Path(path).relative_to(s.REPO).as_posix()})'

def run():
    out=s.OUT;local=s.LOCAL;scaled=out/'scaled_source';restart=out/'restart_coordinates'
    manifests=[out/'build_evidence.json',out/'jacobian_build_evidence.json',out/'freeze.json',out/'preflight_evidence.json',out/'native/freeze.json',out/'failed_steps/evidence.json',out/'failed_trace/evidence.json',scaled/'build_evidence.json',scaled/'jacobian_build_evidence.json',scaled/'freeze.json',scaled/'jacobian_probe_freeze.json',scaled/'preflight_evidence.json',scaled/'dc_evidence.json',scaled/'evidence.json',scaled/'port_calibration/evidence.json',scaled/'input_guards/evidence.json',restart/'build_evidence.json',restart/'freeze.json',restart/'identity_evidence.json',restart/'dc_evidence.json']
    count=0
    for path in manifests:s.a.verify(path);count+=len(s.a.read(path)['input_hashes'])
    dc=s.a.rows(scaled/'dc.csv');response=s.a.rows(scaled/'response.csv');ports=s.a.rows(scaled/'port_calibration/ports.csv');identity=s.a.rows(restart/'identity.csv');repair=s.a.rows(restart/'dc.csv')
    cases=s.a.read(scaled/'contract.json')['cases'];table=[]
    for c in cases:
        p=next(r for r in ports if r['key']==c['key'] and r['amplitude']=='full');rr=[r for r in response if r['key']==c['key']]
        table.append(f"| {c['device']} | {c['vd']} | {float(p['DC_Id_derivative_A_per_um']):+.9e} | {max(float(r['prediction_relative']) for r in rr):.3e} | {min(float(r['signal_over_drift']) for r in rr):.5g} | {sum(r['qualified']=='True' for r in rr)}/2 |")
    environment={k:v for k,v in s.V.environment().items() if k.startswith('VELA_')}
    assert 'VELA_LINEAR_SOLVER' not in environment
    s.a.write(out/'execution_environment.json',dict(VELA_environment=environment,selected_local_linear_backend='sparselu',selection_basis='Default in current src/solver/LinearSolver.cpp; no VELA_LINEAR_SOLVER override in the environment used by drivers.',compiled_optional_backends=['UMFPACK','SPQR'],native_status='not_uploaded_not_run',source_default_changed=False))
    summary=dict(local_source_DC=len(dc),qualified_local_source_DC=sum(r['qualified']=='True' for r in dc),field_responses=len(response),qualified_field_responses=sum(r['qualified']=='True' for r in response),port_responses=len(ports),qualified_port_responses=sum(r['qualified']=='True' for r in ports),maximum_port_prediction_relative=max(float(r['prediction_relative']) for r in ports),restart_identity_corrected_passed=sum(r['qualified']=='True' for r in identity if r['mode']=='corrected'),restart_DC_qualified=sum(r['qualified']=='True' for r in repair),restart_DC=len(repair),native_uploaded=False,native_started=False,production_cpp_changed=False,finite_volume_replacement=False,threshold_changed=False,manifests_verified=len(manifests),direct_hash_entries_verified=count)
    s.a.write(out/'summary.json',summary)
    body=f'''# SimpleMOS PhuMob 局部源校准与重启精度验证

日期：2026-09-09。分支 `codex/simplemos-sdevice-validation`。接续 [低栅压验证](simplemos_phumob_lowvg_validation_2026-09-09.md)。

**局部源的 20/20 次 Vela DC 和 8/8 项漏端响应检查通过；全场响应为 6/8，整体校准仍未放行。** 同时确认正式重启分支的参考值相消问题：隔离修复在 7/7 状态上逐位保留场与残差，但两条原失败路径的重算仍未合格。生产 C++、物理模型、SRH 体积、默认值和接受门槛均未修改。

Sentaurus 的 20 工况输入包已经生成，上传在执行前被自动审批拒绝，尚未上传或启动仿真。已向用户明确请求本批文件和目的地的授权，保留既有授权的上下文。原生同源资格、有限体积替换以及 Enormal/HFS 扩展继续待办。

## 1. 配置、单位与实验边界

采用合格的 n19/n23、Vg=0、Vd=0.05/1 V 四个状态。网格分别 1480/1482 节点，每点 907 个自由 Si 节点、1814 个载流子行。300 K、Boltzmann、plain PhuMob、OldSlotboom、匹配基础 ni、掺杂相关 SRH；显式 `element_box_phumob`、element-box 输运、signed Si Poisson 电荷体积、cell-material 介电系数和 Delaunay box 修正。SRH 仍使用原 all-cell 体积。

本轮 Windows MSYS2 UCRT64 C++20 Release，实际线性后端为 **Eigen SparseLU**，`VELA_LINEAR_SOLVER` 未设置；静态库支持 UMFPACK/SPQR不表示本轮选择它们。上一轮文档的 UMFPACK 标注不能作为本轮后端记录。跨轮迭代数不作严格后端性能 A/B；本轮 legacy/corrected 重启对照使用同一 SparseLU。独立切向使用 SciPy SuperLU 与四次 long-double 缺陷修正。原生输入保留 T-2022.03-SP2、ExtendedPrecision(128)、Super、Digits=12。

仍使用逐行 1e-6、KCL/Id 1e-8；全局源相对条件保留旧 1e-10 下限，未激活的弱源不计独立相对闭合通过。场响应比较全部自由 Si 节点的 ψ/φn/φp 向量：预测差及双幅度一致性 1e-3，偶/奇比 0.01，信号/零漂移至少 100。没有为失败点放宽条件。

## 2. 明确定义的同源方向

在 n23 的 795、792、1114、794、791、1115、1203、1188 八个预选节点取源，n19 按物理坐标独立匹配。冻结每节点积分源 `S_i = R_V(base) × (V_native − V_Vela)`；它是与状态无关的电子/空穴成对源，两个连续性方程的额外状态偏导均为零。SRH 本构及其 Jacobian 不替换。

源在 Vela 的归一化连续性行单位与物理粒子通量之间显式换算；Sentaurus 的 PMI 速率为同一积分源除以原生 Si box 体积。物理积分源单位为 particles/(m·s)，原生速率 cm⁻³s⁻¹，Id 为 A/μm。每工况分别从同一合格基态运行 α=0、±0.001、±0.0005；α=1 仅用于固定状态插入检查，未作有限自洽替换。

初版遗漏归一化前后的 `C0×D0` 换算，四工况固定源插入误差均为 1，被预检拦截；没有运行其 DC。修正版另存 `scaled_source`，原失败完整保留。修正版 4/4 零身份、源插入、Jacobian 逐位不变及固定态边通量不变通过，源相对差最大 {max(float(r['source_relative']) for r in s.a.rows(scaled/'preflight.csv')):.3e}，独立线性闭合最大 {max(float(r['linear_relative']) for r in s.a.rows(scaled/'preflight.csv')):.3e}。9/9 非法输入保护通过。

另有两项工具准备失败保留：坐标字典转换错误在写入工况前停止；生产 runner 不支持旧 `parameter_jacobian` 入口，新增只读工具适配器后再导出完整 J。适配器沿用当前 solver API，没有替换 DC 实现。{link(scaled/'unit_amendment.json','单位修正说明')}、{link(scaled/'preflight.csv','固定源预检')}、{link(scaled/'input_guards/checks.csv','输入保护')}。

## 3. Vela 自洽及端口响应

20/20 DC 通过原逐行和 KCL 条件，最大逐行比 {max(float(r['max_row_ratio']) for r in dc):.3e}，最大 KCL/Id {max(float(r['kcl_over_Id']) for r in dc):.3e}。漏端的完整 J 切向、解析电流梯度与伴随权重三者独立比对，再与实际正负 DC 差分比较：8/8 通过，最大预测相对差 {summary['maximum_port_prediction_relative']:.3e}，最大切向/伴随差 {max(float(r['duality_relative']) for r in ports):.3e}。

| NWell | Vd (V) | dId/dα (A/μm) | 全场预测最大相对差 | 最小场信号/漂移 | 全场资格 |
| --- | --- | --- | --- | --- | --- |
{chr(10).join(table)}

n23 高 Vd 的正向响应约 +3.06656e-18 A/μm，方向和数量级支持八热点 SRH 体积假设；原始漏端差为 −3.00069e-18 A/μm。但此处测量的是冻结源的局部导数，尚无原生同源校准，且 n23 低 Vd 两档场响应的信号/漂移仅 68.64/34.32，低于 100。**不能据此声称有限体积替换已改善 Id，也不能把端口通过覆盖全场失败。** 当前 n23、Vg=0、Vd=1 的正式候选误差仍为 −0.80813524%。{link(scaled/'dc.csv','全部DC')}、{link(scaled/'response.csv','全场响应')}、{link(scaled/'port_calibration/ports.csv','独立端口校准')}。

## 4. 失败初始化与重启相消

旧第一失败态最差为空穴 983，逐行比 2.84227e-5；旧重载最终仅电子 967/1089 超标，最差 5.45097e-6。只读 `evaluateStep` 的完整试探步能降低这些载流子残差，却使 Poisson 范数上升；100 位对已导出的双精度残差平方求和仍判为不下降，因此仅把标量范数比较算得更精确并不足以接收该完整步。

注意：`evaluateStep` 使用直接默认线性求解，不包含正式 `solve()` 的行权重/四次线性修正，也不执行线搜索，因此其方向不能冒充实际 Newton 方向。另做两条各最多三步的生产轨迹诊断，原始和截断更新相同、出现回溯接受，但仍失败。该预算诊断不是第三次正式资格重试。{link(out/'failed_steps/summary.csv','只读步与范数')}、{link(out/'failed_trace/contract.json','实际轨迹范围')}。

进一步发现只读打包器已用 `(old_reference−new_reference)+increment`，正式 warm-start 分支仍用 `(old_reference+increment)−new_reference`。在偏置参考值附近，即使 long double 也会损失远小于参考值 ULP 的增量。隔离程序仅调整电子/空穴两处分组顺序，保留接触投影和一致性保护。

在两个失败态、同偏置合格态及四个 Vg=0 控制态做 14 次零迭代对照：原实现 7/7 未通过状态身份检查；最大增量相对改变约 7.376%，最大残差改变量达到原行尺度的约 99.276%。绝对准费米势变化只有约 1e-21～3e-20 V，密度几乎不变，因此只看绝对势误差会漏检。修正版 7/7 场值、增量、密度和导出载流子残差均逐位保留。该结论仅覆盖这组同偏置参考坐标，尚未覆盖更换偏置参考的通用回归。{link(restart/'identity.csv','零迭代身份对照')}。

修正版又按原 200 次上限分别从两个失败输出重算：

| 输入 | 接受迭代数 | 最差逐行比 | 违规行数 | 结果 |
| --- | --- | --- | --- | --- |
{chr(10).join(f"| {r['label']} | {r['iterations']} | {float(r['max_row_ratio']):.9e} | {r['row_violations']} | {r['failure']} |" for r in repair)}

第二路径在 52 次后剩电子 1089 一行，逐行比 1.045902376e-6，仍超过 1e-6；第 53 次尝试的 13 个回溯候选均未被接受。该行原始/截断电势更新相同，QF 更新约 −5.856e-25 V，线性缺陷已很小，仍有线搜索/残差求值问题。首条路径未能接受第一步。**重启修复已证明消除了输入相消，但未证明消除了收敛失败，暂不移植生产或改变双初始化资格。** {link(restart/'dc.csv','隔离修复DC')}。

## 5. 原生待执行与下一步

本地手册 T-2022.03 第 1214–1216、1230、1240、1266–1268 页支持使用简化 PMI 的 `pmi_float`、节点坐标和额外复合源。本批准备的 PMI 输出是常量，其状态导数为零，保留原 SRH；计划先四次零扰动一致性，再十六次正负响应，并核实八节点命中及 128 位算术。

待上传包为 `build-release/phumob_local_source_20260909/native/input.tgz`，目标 `sentaurus:/tmp/vela_simplemos_phumob_pair_source_20260909/input.tgz`。目录已创建，scp 在执行前被自动审批拒绝；PMI 尚未在虚拟机编译，20 次原生任务均未启动。原因是审批要求对具体网格/状态包和具体目的地授权；已提出明确批次确认，没有换通道绕过。{link(out/'native_transfer_status.json','传输状态')}、{link(out/'native/contract.json','原生合同')}。

下一步优先完成原生同源零控制与正负响应；本地需在重启相消修复的独立数值轴上重新检查低 Vd 零漂移，并定位电子 1089 失败态的实际回溯残差精度及 Poisson 小更新。只有全场、原生响应及严格收敛资格补齐后，才考虑状态相关 SRH 体积插值、双初始化有限替换和模型扩展。本轮没有提交或推送。

## 6. 可复核性

九个新增执行/检查脚本语法通过，源项非法输入 9/9 通过；仅构建隔离装配器、只读工具和重启坐标程序。未改生产 C++，未运行新的全量 CTest，不把上一轮结果计为本轮测试。报告前校验 {len(manifests)} 份清单、{count} 个直接输入哈希项；保留全部构建/配置/插入/响应/重启失败。{link(out/'summary.json','汇总')}、{link(out/'script_qa.json','脚本QA')}、{link(out/'execution_environment.json','本轮环境')}。
'''
    REPORT.write_text(body,encoding='utf-8')
    py_compile.compile(__file__,doraise=True)
    s.d.matrix.freeze(out/'final_evidence.json',manifests+[REPORT,Path(__file__).resolve(),out/'summary.json',out/'script_qa.json',out/'native_transfer_status.json',out/'execution_environment.json'])
    print(summary,flush=True)

if __name__=='__main__':run()
