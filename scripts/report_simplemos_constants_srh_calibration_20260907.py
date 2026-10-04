"""Report qualified small responses without promoting finite replacements."""
from pathlib import Path
import math
import validate_simplemos_constants_srh_calibration_20260907 as t

a=t.a;d=t.d;OUT=t.OUT;LOCAL=t.LOCAL
REPORT=t.p.REPO/'docs/validation/simplemos_constants_srh_response_validation_2026-09-07.md'

def main():
    a.verify(OUT/'validation_evidence.json');response=a.rows(OUT/'response.csv');dc=a.rows(OUT/'dc.csv');pred=a.rows(OUT/'predictions.csv');pre=a.rows(OUT/'preflight.csv')
    assert len(dc)==68 and all(x['qualified']=='True' for x in dc)
    assert len(response)==32 and all(x['qualified']=='True' for x in response)
    manual=t.p.REPO/'build-release/m79_research/sdevice_ug_local_2022.txt';text=manual.read_text(encoding='utf-8')
    assert all(x in text for x in ('double kB = 1.380662e-23;', 'const double eps0 = 8.8542e-14;', 'const double e0 = 1.602192e-19;'))
    q=1.602192e-19/1.602176634e-19-1;eps=8.8542e-12/8.8541878128e-12-1;thermal=(1.380662e-23/1.602192e-19)/(1.380649e-23/1.602176634e-19)-1
    screening=[];volumes=[];ports=[];constants=[]
    for c in t.cases():
        prs={x['axis']:float(x['unit_prediction']) for x in pred if x['key']==c['key']}
        Id=a.read(t.s.LOCAL/c['key']/'generated/config.status.json')['contact_currents_A_per_um']['drain'];native=c['native_Id_A_per_um']
        dq=prs['q_fixed_Vt']*q;de=prs['eps0']*eps;dv=prs['thermal_voltage']*thermal
        screening.append(dict(key=c['key'],native_Id_A_per_um=native,qualified_Id_A_per_um=Id,actual_relative_error=Id/native-1,
            predicted_q_delta_A_per_um=dq,predicted_eps0_delta_A_per_um=de,predicted_thermal_delta_A_per_um=dv,
            predicted_remaining_relative_error=(Id+dq+de+dv)/native-1,interpretation='First-order screening from separately calibrated derivatives. No finite replacement solve.'))
        for job in next(x for x in a.read(OUT/'contract.json')['cases'] if x['key']==c['key'])['fixed_jobs']:
            if Path(job['path']).stem!='functional':continue
            r=a.read(Path(job['path']).with_suffix('.status.json'));err=abs(r['current_A_per_um']/r['contact_current_extractor_A_per_um']-1)
            ports.append(dict(key=c['key'],axis=job['axis'],alpha=job['alpha'],port_relative=err,qualified=err<=1e-8))
            constants.append(dict(key=c['key'],axis=job['axis'],alpha=job['alpha'],**r['diagnostic_constants']))
        if c['device']=='n23':
            geo,xy,el,vol,K,parts,*_=t.s.prior.old.l.g.geometry(c['device'])
            for node in (792,1057):volumes.append(dict(key=c['key'],node=node,signed_Si_volume_m2=float(vol['Si'][node]),original_SRH_volume_m2=float(geo.volumes['all_cell'][node]),signed_to_original_ratio=float(vol['Si'][node]/geo.volumes['all_cell'][node]),finite_change_tested=False))
    assert all(x['qualified'] for x in ports)
    a.write_csv(OUT/'constant_compatibility_screening.csv',screening);a.write_csv(OUT/'node_volume_scope.csv',volumes);a.write_csv(OUT/'fixed_port_consistency.csv',ports);a.write_csv(OUT/'actual_constants.csv',constants)
    summary=dict(a.read(OUT/'summary.json'),maximum_prediction_relative_error=max(float(x['prediction_relative_error']) for x in response),
        maximum_two_amplitude_relative=max(float(x['two_amplitude_relative']) for x in response),maximum_even_over_odd=max(float(x['even_over_odd']) for x in response),minimum_signal_over_zero_drift=min(float(x['signal_over_zero_drift']) for x in response),
        maximum_row_ratio=max(float(x['max_row_ratio']) for x in dc),maximum_KCL_over_Id=max(float(x['kcl_over_Id']) for x in dc),maximum_fixed_source_relative=max(float(x['parameter_source_relative']) for x in pre),maximum_linear_residual_relative=max(float(x['linear_residual_relative']) for x in pre),
        coherent_core_translation_units=47,fixed_port_checks=len(ports),compatibility_q_relative=q,compatibility_eps0_relative=eps,compatibility_Vt_relative=thermal)
    a.write(OUT/'final_summary.json',summary)
    stable=[]
    for r in screening:
        c=next(c for c in t.cases() if c['key']==r['key']);stable.append(f"| {c['device']} | {c['vd']:g} | {100*r['actual_relative_error']:+.8f}% | {100*(r['predicted_remaining_relative_error']-r['actual_relative_error']):+.8f}% | {100*r['predicted_remaining_relative_error']:+.3e}% |")
    srh=[]
    for r in response:
        if not r['axis'].startswith('srh') or r['amplitude']!='full':continue
        c=next(c for c in t.cases() if c['key']==r['key']);srh.append(f"| {c['vd']:g} | {r['node']} | {r['observable']} | {float(r['prediction']):+.9g} | {float(r['actual_odd']):+.9g} | {float(r['prediction_relative_error']):.4g} |")
    report=f'''# SimpleMOS 常数约定与局部 SRH 源的独立响应校准

日期：2026-09-07。Vela 四个 Vg=0.8 V 控制点的 q、ε₀、热电压，以及高 NWell 两个 Vd 下节点 792/1057 的局部 SRH 源方向，完成固定状态源/切向与正负双幅度自洽对照：16/16 方向预检查、68/68 DC、32/32 响应检查通过原门槛。最大响应预测相对差 {summary['maximum_prediction_relative_error']:.7g}，门槛 1e-3。

本轮校准的是 Vela 内部同一参数变化的响应。未更改生产常数或 SRH 体积，未对 Sentaurus 常数、原生 SRH 源实现进行新扰动，也未执行有限替换。已有原生导出和手册足够定义本轮独立参数方向，无新虚拟机文件传输或 sdevice 仿真。

## 1. 隔离与参数定义

沿[自动 box/单元迁移率已合格状态](simplemos_generated_box_mobility_validation_2026-09-07.md)继续；生产几何、掺杂、材料、SRH 模型和边界配置不变。单独重新编译全部 47 个 vela_core 编译单元和诊断 runner，所有 q/kb/ε₀ 使用同一头文件；不将改过常数的少数对象与旧 core 库混合链接。诊断默认零扰动的四点残差及端口与生产程序逐位一致。

三个常数方向分别为：q 与 kb 同乘 (1+α) 以固定 Vt，ε₀ 单独乘 (1+α)，以及固定 q、T 和材料输入仅令 kb 乘 (1+α) 以改变热电压。q 方向覆盖全部电荷、端口前因子、缩放和边界计算，不是仅改变 Poisson 中的一个 q。常数扰动幅度为 ±1e-4、±5e-5；参数源由固定物理状态的残差与端口差分得到，以 1e-5/5e-6 两步复核。

SRH 方向仅把选定自由 Si 节点的原 SRH 积分体积乘 (1+α)，同时改变残差、六项状态导数和求解方程源诊断。Poisson 电子、空穴和掺杂电荷已有独立体积数组，保持不变。冻结模型关闭雪崩与载流子对角下限，因而这个局部 vol 修改没有引入其他源。幅度 ±0.001/±0.0005。

SRH 的源导数独立取自该节点求解方程源项，同时加到电子和空穴行；与实际残差差分复核。主要响应量预先冻结为 φn(792)、φp(1057)，Id 同时记录。不能因端口对极少数载流子几乎不敏感，就用小 Id 差判定其场状态已经正确。

完整冻结合同见[contract.json](../../reference_tcad/simplemos_sentaurus2022/constants_srh_calibration_20260907/contract.json)，实际运行常数见[actual_constants.csv](../../reference_tcad/simplemos_sentaurus2022/constants_srh_calibration_20260907/actual_constants.csv)。

## 2. 独立线性与自洽资格

使用基态完整 Jacobian，独立稀疏 LU 解 J·dx=−∂R/∂α，再以 long double 重算线性残差并修正四次。端口响应为显式端口导数−λ·源，局部势响应由完整 dx 转换为 V；不以局部对角近似替代耦合。最大线性残差相对值 {summary['maximum_linear_residual_relative']:.6g}，门槛 1e-8；参数源的最大双步/独立公式误差 {summary['maximum_fixed_source_relative']:.6g}，门槛 1e-5。全部 {len(ports)} 个固定状态端口与残差功能量一致。

68 次 DC 保留全部 1814 个自由 Si 载流子行和原接受门槛。最差逐行比值 {summary['maximum_row_ratio']:.7g}，门槛 1e-6；最大 KCL/Id={summary['maximum_KCL_over_Id']:.6g}，门槛 1e-8。32 项响应检查全部通过：最大双幅度比例差 {summary['maximum_two_amplitude_relative']:.7g}（门槛 1e-3），最大偶/奇响应比 {summary['maximum_even_over_odd']:.7g}（门槛 0.01），最小响应/零回放漂移 {summary['minimum_signal_over_zero_drift']:.6g}（门槛 100）。零回放电流仍检查原 1e-5 dex 限制。

响应/漂移比只度量本轮回放可见的漂移；基态漂移为零时使用 1e-300 分母保护，极大的数值不代表无限精度或已知绝对误差界限。低 Vd 的微小 SRH 场响应仍通过预测、双幅度和符号三个独立检查。

数据：[DC](../../reference_tcad/simplemos_sentaurus2022/constants_srh_calibration_20260907/dc.csv)、[响应](../../reference_tcad/simplemos_sentaurus2022/constants_srh_calibration_20260907/response.csv)、[预检查](../../reference_tcad/simplemos_sentaurus2022/constants_srh_calibration_20260907/preflight.csv)。

## 3. SRH 局部场响应

下表为 n23、Vg=0.8 V、α=0.001 时的预测变化及正负 DC 中央奇响应，单位 V。

| Vd/V | 节点 | 响应 | 预测/V | 自洽奇响应/V | 预测相对差 |
| ---: | ---: | --- | ---: | ---: | ---: |
{chr(10).join(srh)}

高 Vd 下，SRH 源体积的千分之一变化对应节点 792 的约 −23.458 μV 和节点 1057 的约 +3.438 μV 响应；相同节点在低 Vd 下的响应远小。这验证了这两个独立局部源在当前合格状态上的场耦合，未证明原生 SRH 体积应如何替换。

原生 signed Si 与当前 SRH 体积之比分别约 1.17958、0.217275，已独立列于[体积范围](../../reference_tcad/simplemos_sentaurus2022/constants_srh_calibration_20260907/node_volume_scope.csv)。对应有限变化远超本轮小扰动，不能直接把线性预测当作修复后的场误差。完整局部 ψ/φn/φp、密度、SRH 积分源见[局部表](../../reference_tcad/simplemos_sentaurus2022/constants_srh_calibration_20260907/local_fields.csv)。

## 4. 剩余电流差的常数约定解释

本地 Sentaurus 2022 手册内置 PMI 示例明确列出 q=1.602192e-19 C、kb=1.380662e-23 J/K、ε₀=8.8542e-14 F/cm。此前原生 Poisson/SG 固定状态回放已独立检验这些规定值，本轮没有拟合常数。Vela 的现行 SI 值 q=1.602176634e-19、kb=1.380649e-23、ε₀=8.8541878128e-12 F/m 保持不变。

映射到本轮三个方向，规定差分别为 q: {q:.10g}、ε₀: {eps:.10g}、Vt: {thermal:.10g}。将已校准的一阶导数与这些已知差相乘，可形成下表筛查。最后一列是预测，**不是有限替换后的自洽结果**。

| NWell | Vd/V | 当前实算 Id 误差 | 常数约定预测改变量 | 一阶预测剩余误差 |
| --- | ---: | ---: | ---: | ---: |
{chr(10).join(stable)}

这使常数约定成为解释 Vg=0.8 V 剩余电流差的有力候选；仍需同一完整常数约定下的有限替换与双初始化，才能把预测提升为实算结论。现行 q/kb 是 SI 精确定义值，复现原生旧常数属于兼容性对照，不应据此改写全局物理常数默认值。

## 5. 完成范围与下一步

本轮完成生产几何/迁移率自动形成与四点常数、两节点 SRH 的独立小响应校准。接下来先做手册规定常数组合的有限自洽对照和双初始化，再对节点 792/1057 的 SRH 源体积做分项、联合有限对照；保持两个 Vd 和低 NWell 控制。SRH 的原生同源定义及响应仍需独立核对。Vg=1 V 的常数方向、16 工况和完整 0–1 V 曲线未在此扩展。

运行环境与几何生产报告一致：Windows MSYS2 UCRT64/C++20 Release，实际 DC 为 Eigen SparseLU/COLAMD，独立切向为 SciPy SuperLU/COLAMD 加 long-double 残差修正；HDF5/TDR、UMFPACK、SPQR 编译可用但本轮未切换 DC 后端。网格 μm、状态密度 m⁻³、宽度 1 μm、Id A/μm；300 K、Boltzmann/no-BGN、matched ni、总杂质 Masetti、掺杂 SRH，HFS/表面/Auger/雪崩/DG 关闭。

本轮没有新增生产常数或 SRH 改动，验证由零扰动逐位一致、独立源/J 切向、固定端口与 68 次严格 DC 完成。生产代码测试状况见几何报告：764/780 CTest 通过，16 项既有失败保留。脚本：[隔离编译](../../scripts/build_simplemos_constants_srh_calibration_20260907.py)、[校准执行](../../scripts/validate_simplemos_constants_srh_calibration_20260907.py)、[报告生成](../../scripts/report_simplemos_constants_srh_calibration_20260907.py)。[最终证据](../../reference_tcad/simplemos_sentaurus2022/constants_srh_calibration_20260907/final_evidence.json)关联合同、实际常数、源、状态、资格及本报告。
'''
    assert not REPORT.exists();REPORT.write_text(report,encoding='utf-8')
    d.matrix.freeze(OUT/'final_evidence.json',[Path(__file__).resolve(),REPORT,manual,OUT/'validation_evidence.json',OUT/'final_summary.json']+list(OUT.glob('*.csv')))
    print(summary,flush=True)

if __name__=='__main__':main()
