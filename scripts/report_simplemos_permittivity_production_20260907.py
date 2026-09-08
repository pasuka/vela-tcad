"""Seal production permittivity migration results without changing old evidence."""
import math
from pathlib import Path
import xml.etree.ElementTree as ET
import validate_simplemos_permittivity_production_20260907 as m

a=m.a;d=m.d;v=m.v;LOCAL=m.LOCAL;OUT=m.OUT
REPORT=m.p.REPO/'docs/validation/simplemos_permittivity_production_validation_2026-09-07.md'

def main():
    a.verify(OUT/'validation_evidence.json')
    dc=a.rows(OUT/'dc.csv');comp=a.rows(OUT/'comparison.csv');dual=a.rows(OUT/'dual_initialization.csv')
    migration=a.rows(OUT/'migration_equivalence.csv');ports=a.rows(OUT/'ports.csv');pre=a.rows(OUT/'preflight.csv');jvp=a.rows(OUT/'jvp.csv')
    assert len(dc)==48 and all(r['qualified']=='True' for r in dc)
    assert len(dual)==16 and all(r['qualified']=='True' for r in dual)
    assert all(r['qualified']=='True' for r in migration+ports+pre)
    strong=[r for r in jvp if r['unchanged_SRH_cross_block']=='False' and float(r['step_V'])<1e-5]
    weak=[r for r in jvp if r['unchanged_SRH_cross_block']=='True']
    assert len(strong)==224 and len(weak)==96 and all(r['fixed_state_cross_block_bit_identical']=='True' for r in weak)
    statechecks=[];localchecks=[];pairing=[]
    for c in m.cases():
        geo,mask=v.m.previous.prior.support(c)
        before=d.ordered(m.old.LOCAL/c['key']/'all_native_K/replacement/state.csv',geo.count)
        after=d.ordered(LOCAL/c['key']/'native_cells_native_mu/state.csv',geo.count)
        diff=v.m.previous.prior.delta_states(before,after,mask)
        passed=max(diff[k] for k in ('psi_max_V','phin_max_V','phip_max_V'))<=1e-6 and diff['density_max_relative']<=1e-4
        statechecks.append(dict(key=c['key'],**diff,qualified=passed))
        for node in (1000,1009):
            values={f+'_migration_delta_V':float(v.m.previous.prior.physical(after[node],f)-v.m.previous.prior.physical(before[node],f)) for f in ('psi','phin','phip')}
            localchecks.append(dict(key=c['key'],node=node,**values))
    for vg in (.8,1.):
        for vd in (.05,1.):
            select={(r['device'],r['arm']):r for r in comp if float(r['vg'])==vg and float(r['vd'])==vd}
            for arm in ('cell_local','native_cells','native_cells_native_mu'):
                get=lambda dev,ar:float(select[(dev,ar)]['current_A_per_um'])
                ref=get('n23',arm)/get('n19',arm)
                native=float(select[('n23',arm)]['native_current_A_per_um'])/float(select[('n19',arm)]['native_current_A_per_um'])
                baseline=get('n23','legacy')/get('n19','legacy')
                high=abs(float(select[('n23',arm)]['relative_error']))<abs(float(select[('n23','legacy')]['relative_error']))
                low=abs(float(select[('n19',arm)]['relative_error']))<=abs(float(select[('n19','legacy')]['relative_error']))
                pair=abs(math.log(ref/native))<=abs(math.log(baseline/native))
                pairing.append(dict(vg=vg,vd=vd,arm=arm,high_absolute_improved=high,low_absolute_nonworse=low,paired_log_improved=pair,qualified=high and low and pair))
    a.write_csv(OUT/'state_migration.csv',statechecks);a.write_csv(OUT/'local_state_migration.csv',localchecks);a.write_csv(OUT/'physical_gates.csv',pairing)
    assert all(r['qualified'] for r in statechecks)
    xml=LOCAL/'full_ctest.xml';suite=ET.parse(xml).getroot()
    prior=a.read(v.OUT/'ctest_audit.json');previous={r['test']:r for r in prior['failures']};failures=[]
    for test in suite.findall('testcase'):
        if test.find('failure') is None:continue
        name=test.attrib['name'];output=test.findtext('system-out','')
        lines=[x for x in output.splitlines() if any(word in x for word in ('FAIL:', 'AssertionError','FileNotFoundError'))]
        category='missing_historical_source' if 'FileNotFoundError' in output else ('historical_source_hash' if 'frozen' in output or name=='simplemos_m46_full_matrix_requalification' else 'unclassified')
        failures.append(dict(test=name,category=category,previously_failed=name in previous,evidence=lines))
    total=int(suite.attrib['tests']);failed=int(suite.attrib['failures']);assert len(failures)==failed
    audit=dict(total=total,passed=total-failed,failed=failed,failures=failures,new_failures=sorted({r['test'] for r in failures}-set(previous)),
               previous_total=prior['total'],previous_failed=prior['failed'],old_hashes_updated=False)
    a.write(OUT/'ctest_audit.json',audit)
    assert not audit['new_failures'] and {r['test'] for r in failures}==set(previous)
    native_errors=[abs(float(r['relative_error']))*100 for r in comp if r['device']=='n23' and r['arm']=='native_cells_native_mu']
    worst=max(dc,key=lambda r:float(r['max_row_ratio']))
    summary=dict(status='completed_opt_in_production_migration',DC=48,qualified_DC=48,dual_initialization_groups=16,qualified_dual=16,qualified_ports=48,
        preflight_checks=16,source_relative_max=max(float(r['poisson_relative']) for r in pre),strong_small_step_JVP_checks=224,strong_JVP_relative_max=max(float(r['relative_error']) for r in strong),
        weak_cross_checks=96,weak_cross_evidence='Analytic bit identity at identical states; not a new FD completeness claim.',
        legacy_Id_drift_max=max(float(r['legacy_current_drift']) for r in migration),migration_Id_drift_max=max(float(r['native_cell_migration_current_drift']) for r in migration),
        max_row_ratio=float(worst['max_row_ratio']),max_row_case=worst['key']+'/'+worst['arm'],max_KCL_over_Id=max(float(r['kcl_over_Id']) for r in dc),
        max_dual_phi_V=max(float(r[k]) for r in dual for k in ('psi_max_V','phin_max_V','phip_max_V')),max_dual_density_relative=max(float(r['density_max_relative']) for r in dual),
        max_dual_Id_relative=max(float(r['Id_relative']) for r in dual),max_migration_phi_V=max(float(r[k]) for r in statechecks for k in ('psi_max_V','phin_max_V','phip_max_V')),
        native_mobility_high_NWell_absolute_error_percent=[min(native_errors),max(native_errors)],unit_tests=6,unit_assertions=45,
        ctest_passed=total-failed,ctest_total=total,ctest_failed=failed,ctest_new_failures=0,
        production_code_changed=True,default_changed=False,acceptance_changed=False,native_mobility_promoted_to_production=False,full_curve_started=False,remote_simulations=0)
    a.write(OUT/'summary.json',summary)
    table=[]
    for c in sorted(m.cases(),key=lambda c:(c['vg'],c['device'],c['vd'])):
        rows={r['arm']:r for r in comp if r['key']==c['key']}
        table.append(f"| {c['device']} | {c['vg']:g} | {c['vd']:g} | "+' | '.join(f"{100*float(rows[arm]['relative_error']):+.6f}%" for arm in ('legacy','cell_local','native_cells','native_cells_native_mu'))+' |')
    text=f'''# SimpleMOS 逐材料介电装配生产接入与八点验证

日期：2026-09-07。结论：显式逐单元材料介电装配已接入正式 Poisson、Gummel、耦合 Newton 与耦合接触反应电荷路径；48/48 次 DC、16/16 组双初始化、48/48 组端口一致性检查通过原门槛。新实现复现此前已校准的全部原生 K 结果。

**单独介电修改仍不满足物理放行门槛。** 采用 Vela 自身迁移率时，低 NWell 的四点绝对误差均增大，高 NWell 的 Vg=1 V 两点也变差。加上此前独立校准的原生有效迁移率输入，才能达到高 NWell {min(native_errors):.8f}%–{max(native_errors):.8f}% 的绝对电流误差。因此新策略保持显式选择，现有默认和物理接受标准不变；本轮没有将原生迁移率替换写入正式算法。

## 1. 实现范围

新增 `mesh_geometry.poisson_permittivity_policy`：默认 `legacy_average` 保留旧行为；`cell_material` 使用 Gij=Σcell εcell(d/l)cell。逐单元几何先与自身材料 ε 相乘，再累加到公共边。该系数由一个共享函数供给三个装配器，耦合残差与解析 Jacobian 使用同一缓存，接触反应电荷读取同一残差路径。

`cell_material` 无附加数据时使用 BoxGeometryBuilder 实际生成的单元贡献，包含原有负 cotangent 的回退/截断选择。可独立提供 `poisson_cell_edge_coefficients`，每个 Tri3 单元一条记录，明确 cell_id、原序 node_ids 和局部边 (0,1)/(1,2)/(2,0) 的三个无量纲 d/l。输入不包含 ε，材料数据库仍决定 ε。

输入检查拒绝单元数量/节点顺序错误、非有限或负系数、退化单元及冲突策略；不将错误编号或不支持的几何静默平均。n19/n23 分别核对 2742/2746 单元、8226/8238 个系数。已导出的原生系数足够，本轮没有上传或启动新 sdevice 仿真。

原生 box 修改仍由明确输入提供，尚未成为 Vela 自动生成的通用规则。输运几何、三项 Poisson 电荷体积、SRH 体积、常数、SG 公式和接受门槛保持各自原定义。本次未改变 PN2D BV 模板的原子策略或全局默认。

源码：[共享系数](../../include/vela/equation/AssemblerUtils.h)、[几何构建](../../src/mesh/BoxGeometryBuilder.cpp)、[配置解析](../../src/simulation/ConfigParsing.cpp)。接口见[配置说明](../config_schema.md)。

## 2. 隔离对照与电流误差

下表为 (Id/IdSentaurus−1)×100%。前三列均使用 Vela 正式迁移率；最后一列的原生有效迁移率保持为独立诊断输入，其残差、Jacobian 迁移率缓存和端口三处使用同一因子。它用于检查介电修复迁移的一致性，不代表原生迁移率实现已生产化。

| NWell | Vg/V | Vd/V | 原联合方案 | 局部单元几何 | 原生单元几何 | 原生单元几何＋原生有效迁移率 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
{chr(10).join(table)}

完整 DC 分组：旧策略回放 8 次，局部单元几何 8 次，原生单元几何＋Vela 迁移率双初始化 16 次，原生单元几何＋原生有效迁移率双初始化 16 次。完整记录见[电流表](../../reference_tcad/simplemos_sentaurus2022/permittivity_production_20260907/comparison.csv)与[高低配对门槛](../../reference_tcad/simplemos_sentaurus2022/permittivity_production_20260907/physical_gates.csv)。

仅介电改变时，高低 NWell 误差配对与绝对误差没有全部通过，说明原有误差存在抵消；不能依据介电离散公式正确就将单项修改作为整机精度放行结论。原生 K 与原生有效迁移率的组合在四组配对上均通过高 NWell 改善、低 NWell 不恶化和配对 log 误差改善条件。

## 3. 数值资格与迁移精度

- 固定状态 16/16 组通过。新旧原生 K 装配的 Poisson 差异，相对实际介电扰动源最大 {summary['source_relative_max']:.6g}，门槛 1e-8；旧策略回放的残差保持逐位一致。全部载流子残差和输运边通量保持逐位一致。
- 224 个非弱块的小步长 Jv 检查全部通过，最大相对差 {summary['strong_JVP_relative_max']:.6g}，门槛 1e-4。96 个弱 SRH 交叉块检查保持同态解析值逐位一致，沿用原先独立弱块资格；此次不把弱块差分噪声或单方向通过扩大为任意 Jacobian 完整性结论。
- 48/48 DC 通过全部 1814 个自由载流子行及原全局闭合检查。最差逐行比值 {summary['max_row_ratio']:.9g}，位于 `{summary['max_row_case']}`，仍低于 1e-6，仅约 {(1-summary['max_row_ratio']/1e-6)*100:.2f}% 裕量；最大 KCL/Id 为 {summary['max_KCL_over_Id']:.6g}，门槛 1e-8。
- 16/16 组双初始化通过。最大势差 {summary['max_dual_phi_V']:.6g} V、密度相对差 {summary['max_dual_density_relative']:.6g}、Id 相对差 {summary['max_dual_Id_relative']:.6g}；原门槛依次为 1e-6 V、1e-4、1e-6。
- 新旧原生 K 方案的最大 Id 相对漂移 {summary['migration_Id_drift_max']:.6g}；全自由 Si 节点 ψ/φn/φp 最大迁移差 {summary['max_migration_phi_V']:.6g} V。节点 1000/1009 的迁移差单列于[局部状态表](../../reference_tcad/simplemos_sentaurus2022/permittivity_production_20260907/local_state_migration.csv)。这验证生产接入复现旧合格状态，不代表那些状态剩余的原生场误差已经消失。

## 4. 测试与环境

Release 完整构建成功。新增 6 个 Catch2 用例、45 个断言全部通过：不等界面贡献、材料/节点顺序、均匀材料与长度单位、负 cotangent 分支、法向位移和串联电容解析解、三条装配路径与反应电荷、残差/Jacobian 差分、几何重建与错误输入。见[测试源码](../../tests/test_permittivity_assembly.cpp)。

全套 CTest **{total-failed}/{total} 通过**，仍为原有同一组 {failed} 项失败：13 项历史源码哈希/替代链检查、3 项引用缺失 `tests/test_mos_mixed_material.cpp`。新增失败为 0；旧哈希及失败记录未改写。详见[逐项审计](../../reference_tcad/simplemos_sentaurus2022/permittivity_production_20260907/ctest_audit.json)。因此不声称全套 CTest 已通过。

Windows MSYS2 UCRT64、C++20 Release；实际 DC 使用 Eigen SparseLU/COLAMD。构建确认 HDF5/TDR、UMFPACK、SPQR 可用，DC 未切换后端。n19/n23：Boron 1e17/2e17 cm⁻³，1480/1482 节点，942 个 Si 节点、907 个自由 Si 节点；网格 μm，状态密度 m⁻³，宽度 1 μm，电流 A/μm。

300 K、Boltzmann/no-BGN、matched ni=1.0750038488844236e10 cm⁻³，Masetti 总杂质与掺杂相关 SRH；HFS、表面迁移率、Auger、雪崩、DG 关闭。max_iter=200、reltol=1e-7、abstol=1e-12、stall floor=1e-9，ψ/准费米更新上限 0.35/0.025 V，contact_basin、原标量线搜索、线性迭代修正 4 次。全局连续性 tolerance=1e-6、source_floor=1e-10；低于净源下限的结果不声明任意相对精度。

## 5. 后续边界与证据

本轮完成逐材料介电装配的显式生产接口与八点迁移验证。下一步应解释并实现原生有效迁移率的单元/边形成方式，验证原生 box 修改的自动生成范围；随后沿已有合格状态对 q/ε0/热电压约定和节点 792/1057 的 SRH 源体积分开展独立同扰动校准。既有毫伏级极少数载流子场差仍需保留检查。整个组合与生产回归接续完成前，维持 16 工况和完整 0–1 V 曲线的原放行门槛。

验证脚本：[执行入口](../../scripts/validate_simplemos_permittivity_production_20260907.py)、[报告入口](../../scripts/report_simplemos_permittivity_production_20260907.py)。阶段顺序为 build、prepare、preflight、run、analyze；prepare 后合同与代码哈希冻结，既有输出拒绝覆盖，重跑需新目录与新合同。生产变更前源码和二进制前像已独立保存，不重写旧证据。完整摘要与关联见[最终证据](../../reference_tcad/simplemos_sentaurus2022/permittivity_production_20260907/final_evidence.json)。本轮未提交或合并代码。
'''
    assert not REPORT.exists();REPORT.write_text(text,encoding='utf-8',newline='\n')
    paths=[Path(__file__).resolve(),REPORT,OUT/'validation_evidence.json',OUT/'summary.json',OUT/'ctest_audit.json',OUT/'state_migration.csv',OUT/'local_state_migration.csv',OUT/'physical_gates.csv',
           LOCAL/'full_build.log',LOCAL/'unit_tests.log',LOCAL/'full_ctest.log',xml,m.p.REPO/'tests/test_permittivity_assembly.cpp',m.p.REPO/'CMakeLists.txt',m.p.REPO/'docs/config_schema.md']
    a.write(OUT/'final_evidence.json',dict(status=summary['status'],summary=summary,input_hashes={a.rel(x):a.sha(x) for x in paths}))
    print(summary,flush=True)

if __name__=='__main__':main()
