"""Summarize n23 evidence without mixing model versions or inventing a new sweep."""
import csv
import hashlib
import json
import math
from pathlib import Path

REPO=Path(__file__).resolve().parents[1]
BASE=REPO/'reference_tcad/simplemos_sentaurus2022'
OUT=REPO/'build-release/simplemos_commit_20260908/n23'
REPORT=REPO/'docs/validation/simplemos_n23_idvg_status_2026-09-08.md'


def rows(path):
    with path.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))


def table(head,data):
    return '\n'.join(['| '+' | '.join(head)+' |','| '+' | '.join(['---']*len(head))+' |']+['| '+' | '.join(map(str,row))+' |' for row in data])


def main():
    oldpath=BASE/'fullfield_validation_20260906/native_precision/comparison/points.csv'
    prodpath=BASE/'generated_box_mobility_20260907/current_extraction_fix/comparison.csv'
    contractpath=BASE/'generated_box_mobility_20260907/contract.json'
    constpath=BASE/'constants_srh_finite_20260908/constants/comparison.csv'
    old=sorted([r for r in rows(oldpath) if r['device']=='n23'],key=lambda r:(float(r['vd']),float(r['vg'])))
    current=sorted([r for r in rows(prodpath) if r['device']=='n23' and r['arm']=='generated'],key=lambda r:(float(r['vg']),float(r['vd'])))
    contract=json.loads(contractpath.read_text(encoding='utf-8'));cases={c['key']:c for c in contract['cases']}
    constants=[r for r in rows(constpath) if r['device']=='n23' and r['label']=='manual_constants']
    assert len(old)==42 and len(current)==4 and len(constants)==2
    assert all(r['qualified']=='True' for r in old+current+constants)
    for vd in (.05,1.):
        assert [float(r['vg']) for r in old if float(r['vd'])==vd]==[i/20 for i in range(21)]
    production=[]
    for r in current:
        native=float(cases[r['key']]['native_Id_A_per_um']);vela=float(r['Id_A_per_um']);error=100*(vela/native-1)
        assert math.isclose(error,100*float(r['native_relative_error']),abs_tol=1e-11)
        production.append(dict(vg=float(r['vg']),vd=float(r['vd']),native_A_per_um=native,vela_A_per_um=vela,delta_A_per_um=vela-native,error_percent=error))
    for r in old:
        error=100*(float(r['vela_A_per_um'])/float(r['native_A_per_um'])-1)
        assert math.isclose(error,float(r['relative_error_percent']),abs_tol=1e-10)
    OUT.mkdir(parents=True,exist_ok=True)
    with (OUT/'production_points.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(production[0]));w.writeheader();w.writerows(production)
    for filename,data in [('historical_full_model.csv',old),('manual_constants_points.csv',constants)]:
        with (OUT/filename).open('w',encoding='utf-8',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    maxima=[]
    for vd in (.05,1.):
        subset=[r for r in old if float(r['vd'])==vd];worst=max(subset,key=lambda r:abs(float(r['relative_error_percent'])))
        maxima.append([vd,f"{min(float(r['relative_error_percent']) for r in subset):+.6f}%",f"{float(worst['relative_error_percent']):+.6f}%",worst['vg']])
    history=[]
    for vg in [i/20 for i in range(21)]:
        history.append([f'{vg:.2f}']+[f"{float(next(r for r in old if float(r['vg'])==vg and float(r['vd'])==vd)['relative_error_percent']):+.6f}%" for vd in (.05,1.)])
    inputs=[oldpath,prodpath,contractpath,constpath]
    manifest=dict(inputs={x.relative_to(REPO).as_posix():hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs},
        formula='100*(Id_Vela/Id_Sentaurus-1)',current_points=4,historical_points=42,new_simulations=0,
        distinction='Current Masetti subset has only Vg=.8 and 1; historical 0-1 sweep uses PhuMob/Lombardi/HFS and earlier acceptance scope.')
    (OUT/'inputs.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    body=f'''# n23 Id–Vg 差异：当前点验证与历史全曲线

日期：2026-09-08。当前实现的显式自动 box/Masetti 组合，在已验证的 Vg=0.8/1.0 V、Vd=0.05/1.0 V 四点，相对 Sentaurus 的 Id 差为 **+0.005397%～+0.006752%**。新版整个 0–1 V 扫描尚未执行，不能把四个点外推为完整曲线。

统一误差定义为 **100×(Id_Vela/Id_Sentaurus−1)**；正值表示 Vela 偏大。n23 为 NWell Boron 2e17 cm⁻³，1482 节点、2746 个三角形，300 K，宽度 1 μm，电流单位 A/μm。Vela 使用 Windows UCRT64 Release、Eigen SparseLU/COLAMD。

## 当前生产实现：四个已验证点

这些点采用 Boltzmann/no-BGN、matched ni、总杂质 Masetti、掺杂 SRH，关闭表面迁移率及高场饱和。显式选择 Delaunay box 转移、逐材料介电装配、element-box 输运/迁移率、独立 signed Si Poisson 电荷体积；现代 SI 常数和 SRH 默认体积不变。这是已实现且显式启用的组合，不是全局默认。

{table(['Vg/V','Vd/V','Sentaurus Id / A·μm⁻¹','Vela Id / A·μm⁻¹','ΔId / A·μm⁻¹','相对差'],[[r['vg'],r['vd'],f"{r['native_A_per_um']:.12e}",f"{r['vela_A_per_um']:.12e}",f"{r['delta_A_per_um']:+.6e}",f"{r['error_percent']:+.8f}%"] for r in production])}

该批原始资格为 32/32 DC、8/8 双初始化和 32/32 端口一致性通过；自由 Si 的全部 1814 个载流子行按 1e-6 门槛检查。参见[自动 box 与迁移率验证](simplemos_generated_box_mobility_validation_2026-09-07.md)。本次只整理已有数据，没有新增 DC。

## 手册常数的隔离诊断

在同一简化模型的 Vg=0.8 V 点，同时采用手册 q/kb/ε₀ 后：

{table(['Vg/V','Vd/V','Id 相对差（无量纲）','换算成百分数'],[[r['vg'],r['vd'],f"{float(r['relative_error']):+.9e}",f"{100*float(r['relative_error']):+.9e}%"] for r in constants])}

原生来源路径的最终资格以冻结延续的严格双初始化记录为准，见[常数与局部 SRH 有限验证](simplemos_constants_and_srh_finite_validation_2026-09-08.md)。生产常数没有改成手册数值；这些接近导出精度的点差不能解释为整个器件或完整曲线具有同样精度。

最近报告中的 **0.104723%** 是邻域 SRH 小扰动的物理场响应与 Jacobian 预测之间的差，**不是 Id–Vg 相对 Sentaurus 的电流误差**。该响应偏差几乎可由剩余空穴残差投影解释，但原响应校准仍为失败，见[稳定范数差及 SRH 校准](simplemos_stable_merit_step_and_srh_calibration_2026-09-08.md)。

## 历史完整 0–1 V 曲线

最近一次覆盖整个 0–1 V 的 n23 数据是 2026-09-06 的 42 点（两个 Vd、每条 21 点）。它采用 **PhuMob + Lombardi + 向量准费米梯度高场饱和**，以及当时的几何/体积实现。它与上述 Masetti 简化模型的原生参考电流也不同，不能把两组百分数直接当成同一物理模型的纯算法前后改善。

历史资格包含载流子行筛选及带 floor 的全局闭合；当时的 qualified 不等于后来全部 1814 行的严格资格，少数载流子场问题仍保留。n23 参考采用加严后的 Sentaurus 128-bit 保存态重闭合。详见[历史扫描及资格边界](simplemos_jacobian_and_fullfield_followup_2026-09-06.md)。

{table(['Vd/V','最小差','最大差','最大差对应 Vg/V'],maxima)}

{table(['Vg/V','Vd=0.05 V 相对差','Vd=1 V 相对差'],history)}

因此，历史完整模型在 0.5–1 V 同样有较大差异；当前简化模型的新版完整曲线差异仍待重算。这两个结论应同时保留。

## 数据与复核

本表由 [汇总脚本](../../scripts/summarize_simplemos_n23_idvg_20260908.py)重新从以下本地证据读取，并复核百分数公式、两条历史曲线各 21 个唯一 Vg、四个当前点的原生电流映射及资格字段：

{chr(10).join('- `'+x.relative_to(REPO).as_posix()+'`，SHA256 `'+manifest['inputs'][x.relative_to(REPO).as_posix()]+'`。' for x in inputs)}

原始仿真和机器生成比较数据按仓库规则保留在本地，本轮不纳入 Git 提交。上面的完整数字表随文档提交；本次读取的输入哈希、汇总 CSV 位于 `build-release/simplemos_commit_20260908/n23/`。诊断脚本重放仍需要相应本地证据和忽略目录中的网格/状态/程序，单独检出代码不能恢复全部历史实验。
'''
    REPORT.write_text(body,encoding='utf-8')
    print(json.dumps(dict(production=production,historical_extrema=maxima,report=str(REPORT)),indent=2),flush=True)


if __name__=='__main__':main()
