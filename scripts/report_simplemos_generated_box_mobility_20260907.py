"""Seal geometry/mobility production results and retained failure provenance."""
import math
import xml.etree.ElementTree as ET
from pathlib import Path
import complete_simplemos_generated_box_mobility_20260907 as done

s=done.s;a=s.a;d=s.d;v=s.v;OUT=s.OUT;LOCAL=s.LOCAL
REPORT=s.p.REPO/'docs/validation/simplemos_generated_box_mobility_validation_2026-09-07.md'

def main():
    a.verify(OUT/'validation_evidence.json');dc=a.rows(OUT/'dc.csv');comp=a.rows(OUT/'comparison.csv')
    assert len(dc)==32 and all(x['qualified']=='True' for x in dc)
    statechecks=[];pairing=[]
    for c in s.prior.cases():
        for arm in ('legacy','native_cells','generated','generated_native_seed'):
            before=done.INITIAL_LOCAL/c['key']/arm;after=LOCAL/c['key']/arm
            identical=a.sha(before/'state.csv')==a.sha(after/'state.csv')
            Id=a.read(after/'config.status.json')['contact_currents_A_per_um']['drain'];func=a.read(after/'post/functional.status.json')['current_A_per_um']
            statechecks.append(dict(key=c['key'],arm=arm,first_attempt_state_bit_identity=identical,DC_port_vs_residual_functional=abs(Id/func-1),qualified=identical and abs(Id/func-1)<=1e-8))
    for vg in (.8,1.):
        for vd in (.05,1.):
            rows={(r['device'],r['arm']):r for r in comp if float(r['vg'])==vg and float(r['vd'])==vd}
            rel=lambda device,arm:float(rows[(device,arm)]['native_relative_error'])
            high=abs(rel('n23','generated'))<abs(rel('n23','legacy'));low=abs(rel('n19','generated'))<=abs(rel('n19','legacy'))
            paired=abs(math.log((1+rel('n23','generated'))/(1+rel('n19','generated'))))<=abs(math.log((1+rel('n23','legacy'))/(1+rel('n19','legacy'))))
            pairing.append(dict(vg=vg,vd=vd,high_absolute_improved=high,low_absolute_nonworse=low,paired_log_improved=paired,qualified=high and low and paired))
    assert all(x['qualified'] for x in statechecks+pairing)
    a.write_csv(OUT/'extraction_fix_identity.csv',statechecks);a.write_csv(OUT/'physical_gates.csv',pairing)
    suite=ET.parse(done.INITIAL_LOCAL/'full_ctest.xml').getroot();previous=a.read(s.prior.OUT/'ctest_audit.json');old={x['test']:x for x in previous['failures']};failures=[]
    for t in suite.findall('testcase'):
        if t.find('failure') is None:continue
        text=t.findtext('system-out','');name=t.attrib['name']
        failures.append(dict(test=name,category=old.get(name,{}).get('category','new_failure'),previously_failed=name in old,evidence=[x for x in text.splitlines() if any(w in x for w in ('FAIL:','AssertionError','FileNotFoundError'))]))
    total=int(suite.attrib['tests']);assert {x['test'] for x in failures}==set(old)
    audit=dict(total=total,passed=total-len(failures),failed=len(failures),failures=failures,new_failures=[],old_evidence_updated=False,
               after_extraction_fix='7 targeted CTests passed plus all 32 DC states bit-identical and terminal outputs independently agree. Full 780-test run preceded the two-line extraction fix.')
    a.write(OUT/'ctest_audit.json',audit)
    jvp=a.rows(OUT/'jvp.csv');strong=[x for x in jvp if x['unchanged_SRH_cross_block']=='False' and float(x['step_V'])<1e-5];weak=[x for x in jvp if x['arm']=='generated' and x['unchanged_SRH_cross_block']=='True']
    assert len(strong)==336 and len(weak)==48 and all(x['fixed_state_cross_block_bit_identical']=='True' for x in weak)
    dual=a.rows(OUT/'dual.csv');migration=a.rows(OUT/'state_equivalence.csv');identity=a.rows(OUT/'cell_identity.csv')
    summary=dict(a.read(OUT/'summary.json'),status='completed_opt_in_generated_cell_box_and_mobility',physical_gates=4,qualified_physical_gates=4,
        geometry_max_absolute=max(float(x['coefficient_max_absolute']) for x in identity),cell_mu_max_relative=max(float(x['cell_mu_max_relative']) for x in identity),
        strong_JVP_checks=len(strong),strong_JVP_max_relative=max(float(x['relative_error']) for x in strong),weak_bit_identity_checks=48,
        max_dual_phi_V=max(float(x[k]) for x in dual for k in ('psi_max_V','phin_max_V','phip_max_V')),max_dual_density_relative=max(float(x['density_max_relative']) for x in dual),
        max_migration_phi_V=max(float(x[k]) for x in migration for k in ('psi_max_V','phin_max_V','phip_max_V')),retained_failed_DC=24)
    a.write(OUT/'final_summary.json',summary)
    table=[]
    for c in sorted(s.prior.cases(),key=lambda x:(x['vg'],x['device'],x['vd'])):
        rows={r['arm']:r for r in comp if r['key']==c['key']}
        table.append(f"| {c['device']} | {c['vg']:g} | {c['vd']:g} | {100*float(rows['legacy']['native_relative_error']):+.6f}% | {100*float(rows['generated']['native_relative_error']):+.6f}% |")
    report=f'''# SimpleMOS 原生等效单元迁移率与自动 box 几何验证

日期：2026-09-07。自动生成方案完成显式生产接入，32/32 DC、8/8 双初始化、8/8 原生表格方案迁移、32/32 端口一致性与四组高低 NWell 配对门槛通过。新方案无需原生几何或迁移率输入表。该结论仅覆盖已验证二维网格、300 K constant/Masetti 子集；全局默认、SRH 体积、常数及接受门槛不变。

## 1. 单元与边形成方式

每个三角形的原始局部边系数为 g(T,ij)=cot(对角)/2。当同一区域的两相邻单元分别贡献负、正 g 且总和非负时，将负贡献置零，另一贡献改为两者之和。n19/n23 分别发生 78/76 次这样的转移。全边系数和全节点 signed box 体积守恒，但单元内 M 的总和可能偏离该单元坐标面积，因而不能用坐标面积作迁移率平均分母。

处理后 M(T,i)=Σj[g(T,ij)·l(ij)²]/4；先按每个顶点的 Nd+Na 计算 Masetti μi，再计算 μT=Σi(Mi·μi)/ΣiMi。输运只累计 Si：μedge=ΣT(gT·μT)/ΣTgT。端口、耦合残差、解析 Jacobian 缓存和探针使用同一函数。这个模型只依赖固定掺杂，所以迁移率对 ψ、φn、φp 的偏导为零；不支持的场、表面或载流子依赖模型明确拒绝。

这是依据本地 Sentaurus 2022 手册第 38 章 pp.1177–1185 的单元平均/box 描述、原生 Measure/Coefficients 与独立单元 μ 导出验证出的二维子集实现。不能将其扩大为 Sentaurus 一般非 Delaunay、加权或三维交叠算法。全局非 Delaunay、非流形、显著负贡献的区域边界会拒绝；近零分支仅使用 128 个 double epsilon 的舍入界限。

配置分别为 `mesh_geometry.cell_box_policy=delaunay_transfer`、`poisson_permittivity_policy=cell_material`、`solver.region_resolved_interface_assembly.transport_edge_geometry=element_box`、`solver.mobility.edge_averaging=element_box`；三项 Poisson 电荷另选 `signed_transport`。不以布尔总开关替代组合。Gummel 载流子路径尚不支持新平均并明确拒绝。详见[配置说明](../config_schema.md)、[几何实现](../../include/vela/mesh/DelaunayBox.h)、[共享公式](../../include/vela/equation/AssemblerUtils.h)。

## 2. 八点电流与资格

表中误差为 (Id/IdSentaurus−1)×100%。旧联合方案为本轮相同种子、相同严格条件下的生产回放。

| NWell | Vg/V | Vd/V | 旧联合方案 | 自动 box＋单元/边迁移率 |
| --- | ---: | ---: | ---: | ---: |
{chr(10).join(table)}

24 组固定状态检查、16 组全单元检查通过。自动 g 最大绝对差 {summary['geometry_max_absolute']:.6g}，M 最大绝对差 3.03577e-18 μm²，μT 最大相对差 {summary['cell_mu_max_relative']:.6g}。这些全单元等式排除了通过拟合端口电流构造 μ 的做法。

336 个非弱小步长 Jv 检查通过，最大相对差 {summary['strong_JVP_max_relative']:.6g}，门槛 1e-4。48 个自动/导入几何间的弱 SRH 交叉块解析值同态逐位一致，沿用此前独立弱块资格；不将此等同任意 Jacobian 的新完整性验证。

全部 DC 保留 1814 个自由 Si 载流子行，最差比值 {summary['maximum_row_ratio']:.9g}，门槛 1e-6，位于 n23/Vg=1/Vd=1/generated（节点 1000 电子行）；裕量约 17.5%。最大 KCL/Id={summary['maximum_KCL_over_Id']:.6g}，门槛 1e-8。双初始化最大势差 {summary['max_dual_phi_V']:.6g} V、密度相对差 {summary['max_dual_density_relative']:.6g}，门槛 1e-6 V/1e-4；Id 门槛 1e-6。相比此前原生表格合格态，最大 Id 相对迁移差 {summary['maximum_migration_Id_relative']:.6g}，门槛 1e-8。

结果见[32 次 DC](../../reference_tcad/simplemos_sentaurus2022/generated_box_mobility_20260907/current_extraction_fix/dc.csv)、[双初始化](../../reference_tcad/simplemos_sentaurus2022/generated_box_mobility_20260907/current_extraction_fix/dual.csv)、[配对门槛](../../reference_tcad/simplemos_sentaurus2022/generated_box_mobility_20260907/current_extraction_fix/physical_gates.csv)。毫伏级极少数载流子场差仍属于剩余问题，本轮迁移精度不等于这些场差已经解决。

## 3. 首批失败与修复

首批 32 DC 中 8 个旧方案通过，24 个新方案在求解并导出状态后，由最终 ContactCurrent 构造器拒绝。原因是 `runNewtonSolveFromState` 没有把 regionResolvedInterfaceAssembly 传给端口提取器。已补齐传递并用独立新目录重跑；32 个状态与首批逐位一致，所有 DC 端口与残差功能量一致。该问题属于生产配置传递遗漏，不能把首批仅有合格行比值当作整体验证通过。

首批状态、stderr、失败资格、源码/二进制前像和哈希均保留在[首次证据](../../reference_tcad/simplemos_sentaurus2022/generated_box_mobility_20260907/failed_attempt_evidence.json)。修复后的[独立一致性检查](../../reference_tcad/simplemos_sentaurus2022/generated_box_mobility_20260907/current_extraction_fix/extraction_fix_identity.csv)证明修复没有改变已算出的状态或接受门槛。

## 4. 测试和运行条件

新增 4 个 Catch2 用例、35 断言通过，覆盖全边/节点守恒、非支持几何拒绝、非线性 Masetti 单元权重和三块 Newton Jv。完整 Release 构建成功；全套 CTest 764/780 通过，16 项失败仍为此前 13 个历史源码哈希/替代链检查和 3 个缺失历史源文件引用，新增失败为零。完整测试发生在两行端口传递修复之前；修复之后 7 项相关 CTest 与全部 32 DC/端口复核通过。详见[测试审计](../../reference_tcad/simplemos_sentaurus2022/generated_box_mobility_20260907/current_extraction_fix/ctest_audit.json)，不声明全套测试已通过。

Windows MSYS2 UCRT64、C++20 Release，实际 DC 使用 Eigen SparseLU/COLAMD；HDF5/TDR、UMFPACK、SPQR 编译可用。n19/n23 为 1480/1482 节点、2742/2746 三角形；Boron 1e17/2e17 cm⁻³，网格 μm，状态密度 m⁻³，宽度 1 μm，Id 单位 A/μm。

300 K、Boltzmann/no-BGN、matched ni=1.0750038488844236e10 cm⁻³、总杂质 Masetti、掺杂 SRH，HFS/表面/Auger/雪崩/DG 关闭。max_iter=200、reltol=1e-7、abstol=1e-12、stall floor=1e-9、ψ/准费米更新上限 0.35/0.025 V、contact_basin、原标量线搜索、线性修正 4 次；原净源闭合门槛 1e-6、source_floor=1e-10，低于净源下限不声明任意相对精度。

本阶段使用既有原生导出，无新虚拟机传输或 sdevice 运行。常数与节点 792/1057 SRH 独立校准在单独证据目录执行，不属于本报告的已通过结论。16 工况和完整 0–1 V 曲线仍遵循原放行门槛。代码未提交/合并。

执行入口：[初始验证](../../scripts/validate_simplemos_generated_box_mobility_20260907.py)、[修复后重跑](../../scripts/complete_simplemos_generated_box_mobility_20260907.py)、[报告生成](../../scripts/report_simplemos_generated_box_mobility_20260907.py)。[最终证据](../../reference_tcad/simplemos_sentaurus2022/generated_box_mobility_20260907/current_extraction_fix/final_evidence.json)关联冻结合同、代码、配置、状态、测试与报告。
'''
    assert not REPORT.exists();REPORT.write_text(report,encoding='utf-8')
    files=[Path(__file__).resolve(),REPORT,OUT/'validation_evidence.json',OUT/'final_summary.json',OUT/'ctest_audit.json',OUT/'physical_gates.csv',OUT/'extraction_fix_identity.csv',done.INITIAL_LOCAL/'full_ctest.xml',done.INITIAL_LOCAL/'full_ctest.log',done.INITIAL_LOCAL/'unit_tests.log',done.INITIAL_LOCAL/'extraction_fix_tests.log']
    d.matrix.freeze(OUT/'final_evidence.json',files);print(summary,flush=True)

if __name__=='__main__':main()
