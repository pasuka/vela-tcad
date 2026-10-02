"""Record the source-aware rejection audit and gated Poisson candidate results."""
from pathlib import Path
import ast
import validate_simplemos_poisson_half_source_20260910 as h
import validate_simplemos_combined_source_20260910 as combined
c=h.c;s=c.s;v=c.v
OUT=v.OUT

def run():
    manifests=[v.OUT/'replay_evidence.json',v.OUT/'precision/evidence.json',c.OUT/'kernel_check/evidence.json',c.OUT/'dc_evidence.json',h.OUT/'dc_evidence.json',c.OUT/'disabled_identity/evidence.json',c.OUT/'response/evidence.json',c.OUT/'response/jvp/evidence.json',c.OUT/'contact_step/dc_evidence.json',c.OUT/'contact_step/identity_rows/evidence.json',combined.OUT/'combined_evidence.json']
    for f in manifests:s.a.verify(f)
    rows=s.a.rows(c.OUT/'dc.csv');half=s.a.rows(h.OUT/'dc.csv');packed=[r for r in rows if r['mode']=='packed']+half
    checks=s.a.rows(v.OUT/'precision/checks.csv');linear=s.a.rows(v.OUT/'precision/linear.csv');identity=s.a.rows(v.OUT/'identity.csv')
    files=[s.REPO/'scripts'/name for name in ('audit_simplemos_lowvd_source_precision_20260910.py','analyze_simplemos_lowvd_source_precision_20260910.py','validate_simplemos_poisson_precision_20260910.py','check_simplemos_poisson_precision_20260910.py','validate_simplemos_poisson_half_source_20260910.py','calibrate_simplemos_packed_poisson_response_20260910.py','check_simplemos_packed_poisson_jvp_20260910.py','validate_simplemos_exact_contact_step_20260910.py','validate_simplemos_combined_source_20260910.py','report_simplemos_source_precision_20260910.py')]
    for f in files:ast.parse(f.read_text(encoding='utf-8'),filename=str(f))
    summary=dict(identity_replays=len(identity),qualified_identity=sum(r['identity']=='True' and r['repeat_identity']=='True' for r in identity),precision_snapshots=len(checks),max_precision_agreement=max(float(r['precision_agreement']) for r in checks),max_rounded_reconstruction=max(float(r['reconstruction']) for r in checks),max_linear_defect_over_residual=max(float(r['linear_defect_norm'])/float(r['base_norm']) for r in linear),kernel_DC=sum(r['mode']=='kernel' for r in rows),kernel_qualified=sum(r['mode']=='kernel' and r['qualified']=='True' for r in rows),packed_DC=len(packed),packed_qualified=sum(r['qualified']=='True' for r in packed),production_modified=False,uniform_16_point_qualification=False,full_phumob_curves=False)
    combo=s.a.read(combined.OUT/'combined_summary.json');summary['combined']=combo
    summary['poisson_direction_checks']=s.a.read(c.OUT/'response/jvp/summary.json')
    s.a.write(OUT/'stage_summary.json',summary)
    table=['| NWell | 源幅度 | 求值方式 | DC 资格 | 最大载流子行比 | KCL/Id | 迭代 |','| --- | --- | --- | --- | --- | --- | --- |']
    for r in rows+half:table.append(f"| {r['device']} | {r['alpha']} | {r['mode']} | {'通过' if r['qualified']=='True' else '失败'} | {float(r['max_row_ratio']):.9g} | {float(r['kcl_over_Id']):.6g} | {r['iterations']} |")
    combo_rows=s.a.rows(combined.OUT/'dc.csv');ct=['| 工况 | DC 通过 | 最大载流子行比 | 最大 KCL/Id | 迭代范围 |','| --- | --- | --- | --- | --- |']
    for key in dict.fromkeys(r['key'] for r in combo_rows):
        group=[r for r in combo_rows if r['key']==key]
        ct.append(f"| {key} | {sum(r['qualified']=='True' for r in group)}/5 | {max(float(r['max_row_ratio']) for r in group):.6g} | {max(float(r['kcl_over_Id']) for r in group):.6g} | {min(int(r['iterations']) for r in group)}–{max(int(r['iterations']) for r in group)} |")
    native_zero=combined.ORIGINAL_NATIVE.parent/'zero_points.csv';native=s.a.rows(native_zero);zero_rows=[]
    for row in combo_rows:
        if row['label']!='zero':continue
        target=next(r for r in native if r['key']==row['key']);current=float(row['current_A_per_um']);reference=float(target['Id_A_per_um'])
        zero_rows.append(dict(key=row['key'],vela_A_per_um=current,native_A_per_um=reference,signed_error_percent=(current/reference-1)*100,qualified=row['qualified']=='True'))
    s.a.write_csv(OUT/'zero_current_comparison.csv',zero_rows)
    n23_high=next(r for r in zero_rows if r['key']=='m65_n23_vd_1p000000_endpoint')
    n23_low=next(r for r in zero_rows if r['key']=='m65_n23_vd_0p050000_endpoint')
    report=s.REPO/'docs/validation/simplemos_lowvd_source_precision_2026-09-10.md'
    text=f'''# SimpleMOS 低 Vd 带源 Poisson 精度定位与隔离验证

日期：2026-09-10。工作分支：`codex/simplemos-sdevice-validation`。

本轮完成实际带源拒绝轨迹定位、独立高精度参考、两种 Poisson 求值候选及双幅度低 Vd 自洽对照。求和精度候选为 {summary['kernel_qualified']}/{summary['kernel_DC']} DC 通过；从保存坐标高精度重算物理量的候选为 {summary['packed_qualified']}/{summary['packed_DC']} 通过。再加入精确接触行更新，扩展四个控制工况，共 {combo['qualified_DC']}/{combo['DC']} DC、{combo['qualified_responses']}/{combo['responses']} 全场响应、{combo['qualified_native_ports']}/{combo['native_port_responses']} 原生端口响应通过。未将不同实现的合格状态拼接成统一 16/16，也未开始完整 PhuMob 曲线、Enormal 或高场饱和。

## 配置与冻结条件

复用 n19/n23 的既有二维网格（分别 1480/1482 节点），Vg=0 V、300 K，Boltzmann、OldSlotboom、plain PhuMob `element_box_phumob`、掺杂相关 SRH。低 Vd 隔离使用 0.05 V，联合控制另覆盖 Vd=1 V。保持已验证联合几何、单元介电系数及原 SRH 体积。电流单位 A/μm。

八个局部节点的等量电子/空穴源、正负幅度 ±0.001/±0.0005 与前轮相同；零源基态来自已冻结的第二次零漂移重启。源独立于状态，不新增 Jacobian 项。该源试验不等于有限 SRH 体积替换。

实际 Windows UCRT64 Release、Eigen SparseLU、四次线性修正，稳定范数差和重启坐标修复仍为诊断叠加实现。保留 200 次迭代、绝对/相对残差 1e-12/1e-7、既有残差平台 1e-9 退出支路、载流子逐行 1e-6、KCL/Id 1e-8、行缩放、步长上限和回溯网格。平台接受仍须通过原逐行资格，不能用绝对小残差覆盖载流子失败。未改变物理常数或接受门槛。

## 实际拒绝原因

四个正负全幅失败轨迹均逐位复现最终状态、退出状态与重复残差，共 56 个快照。对每组基态及实际 trial 0/5/12 共 16 个快照做全部方程行的 Decimal 60/100 位检查，明确在两个载流子方程中保留同一独立源项。精度一致性最大 {summary['max_precision_agreement']:.6g}，舍入项重建误差最大 {summary['max_rounded_reconstruction']:.6g}，均按通量/源项规模归一化。

原始步与截断步四组完全相同，独立完整 J 线性缺陷/残差范数最大 {summary['max_linear_defect_over_residual']:.6g}。四个实际全步均已将载流子行比降至约 2.98e-11～1.16e-10，但 Poisson 块使总范数增加而拒绝；不能归因为载流子不下降或更新截断。

n19/n23 正扰动全步最大电势更新约 9.03e-17/8.99e-17 V。最小回溯时自由节点电势更新全部被状态加法舍入吸收，剩下 39 个发生极小变化的电势节点均为 Dirichlet 接触。只提高范数比较精度不足以修复这类候选。

保持生产双精度内核操作数的高精度求值，四组全步仍不下降；从实际保存坐标高精度形成 ψ/n/p 的诊断使两个 n23 全步下降，但 n19 仍不下降。无舍入的数学试探状态四组均下降。后者不是可直接接受的生产状态，仅用于区分求值误差与更新表示误差。

## 两种独立求值候选

- `kernel`：用 binary128 重算 Poisson 边差、乘积与求和，ψ/n/p 仍使用现有双精度值。
- `packed`：额外从实际双精度保存坐标、准费米参考和增量，以 binary128 重算物理电势及 Boltzmann 浓度，再形成 Poisson 残差。

两者保持原双精度状态存储及原 Jacobian；迁移率、连续性通量、端口提取、源项和体积不变。采用 Boost `cpp_bin_float_quad` 的软件 binary128，仅作为有保护的隔离候选，不是全系统高精度求解。其 47,392 个 Poisson 节点检查通过独立 100 位参考，最大按项规模误差约 2.31e-33，小于冻结 1e-25。关闭选择器时，n19/n23 正扰动重放均逐位保持原轨迹。

{chr(10).join(table)}

所有失败和原接受预算保留。DC 通过仍不自动等于双幅度全场响应、初始化不变性或生产实现资格；后续必须依据完整对照和相应门槛决定是否移植。

Poisson 对 ψ/φn/φp 三个方向的三档步长检查，两个低 Vd 基态共 18/18 通过；最大相对差 1.39453e-9，小于 1e-6。该结论限于此次 Poisson 块，不能代替任意载流子交叉块审计。重新导出完整 J、固定源、输运通量的预检通过；源项没有新增导数。

## 接触约束与四点联合验证

四个实际失败基态各有 246 个严格单位矩阵边界行，其中 49 个线性求解更新带有约 1e-45 的约束舍入偏差。只把这些单位行的步长设为其残差的相反数，即精确求解原边界方程，自由行步长逐位不变；约束行线性缺陷变为零。没有修改边界偏置、自由状态、接受门槛或物理模型。

三个隔离对照全部通过：n19 负全幅从 97 次迭代降为 10 次；n23 零源从 171 次降为 3 次，正全幅从 200 次降为 3 次。两项 n23 电流逐位相同，n19 相对变化约 1.4e-16。原长迭代轨迹仍保留，不能以此覆盖初始失败。

随后固定该数值组合，在 n19/n23 × Vd=0.05/1 V 上重算完整 J、源插入、正负双幅度及零源。继续复用已经运行并验证的原生同源实验，未新增原生仿真。

{chr(10).join(ct)}

联合全场按原预测相对误差 1e-3、双幅度 1e-3、偶/奇比 1e-2、信号/漂移 ≥100 判断；端口按同扰动原生相对误差 1e-3 判断。原生漏端/衬底响应最大相对差 {combo['max_native_relative']:.8g}。注意：通用源分析器的 `summary.json` 保留旧字段 `native_response_qualified=false`；本轮独立原生端口结果应读取 `native_ports.csv` 与 `combined_summary.json` 的明确端口计数，不能把通用字段解释为这 16 项端口检查失败。

当前完成数值候选与源响应资格，尚未把这些诊断叠加实现移入正式求解器；下一阶段是生产移植、统一实现下的 16 个控制点与双初始化复验，通过后再扩展四条 PhuMob 完整曲线及恢复 Enormal/HFS。此次合格源响应不等于有限 SRH 体积替换已获验证。

## 证据与范围

- [输出身份](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/identity.csv)、[分块实际范数差](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/precision/actual_deltas.csv)、[高精度参考](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/precision/merit.csv)。
- [独立 C++ 内核验证](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/kernel_check/summary.json)、[全幅及零源 DC](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/dc.csv)、[半幅 DC](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/half_amplitude/dc.csv)。
- [Poisson 方向导数](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/response/jvp/checks.csv)、[低 Vd 全场响应](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/response/response.csv)、[接触行隔离](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/contact_step/dc.csv)。
- [四工况联合 DC](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/contact_step/full_source/dc.csv)、[联合全场响应](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/contact_step/full_source/response.csv)、[原生端口响应](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/contact_step/full_source/native_ports.csv)。
- [阶段汇总](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/stage_summary.json)。

本轮新增诊断脚本、诊断 C++ 头文件和报告；未修改生产 C++，未新增 Sentaurus 仿真，未提交或推送。联合零源 n23、Vg=0 的电流误差为 {n23_low['signed_error_percent']:+.8f}%（Vd=0.05 V）和 {n23_high['signed_error_percent']:+.8f}%（Vd=1 V）；后者 Vela 电流为 {n23_high['vela_A_per_um']:.16g} A/μm，原生为 {n23_high['native_A_per_um']:.16g} A/μm。数值修复没有消除既有深关断电流差。[当前零源电流对比](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/zero_current_comparison.csv)。完整 PhuMob 16 点、双初始化及四条曲线仍待统一实现后验证。
'''
    assert not report.exists(),report
    report.write_text(text,encoding='utf-8')
    s.d.matrix.freeze(OUT/'stage_evidence.json',manifests+files+[c.HEADER,native_zero,OUT/'zero_current_comparison.csv',OUT/'stage_summary.json',report])
    print(summary,flush=True)

if __name__=='__main__':run()
