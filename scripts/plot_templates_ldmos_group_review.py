"""Create source-backed static figures for the 12–17 September group review.

No solver is executed. Run from this worktree with UCRT64 Python and -X utf8.
All plotted data, provenance hashes and PNG/PDF exports accompany report.md.
"""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.collections import LineCollection
from matplotlib.colors import SymLogNorm, TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch
from matplotlib.tri import Triangulation
import numpy as np
import h5py

from analyze_templates_ldmos_stage4_d5 import read_curve
from extract_templates_ldmos_d0_temperature import values, geometry_hash

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'docs/reports/ldmos_review_2026-09-12_17'
B = ROOT/'reference_staging/templates_ldmos_symbolic_vm_20260917'
I = ROOT/'reference_staging/templates_ldmos_r11_isothermal_20260917'
N = ROOT/'reference_staging/templates_ldmos_d0_electrothermal_20260912/native_fields_r1'
USED = {}
EXPORTS = []


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def track(p, expected=None):
    p = Path(p).resolve()
    digest = sha(p)
    if expected is not None:
        assert digest == expected, f'Changed source: {p}'
    USED[str(p)] = digest
    return p


def read(p, expected=None):
    return json.loads(track(p, expected).read_text(encoding='utf-8'))


def table(name, rows):
    p = OUT/'data'/f'{name}.csv'
    with p.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def finish(fig, stem, note):
    fig.text(.06, .025, note, fontsize=9, color='#526170', va='bottom')
    fig.subplots_adjust(left=.07, right=.96, top=.84, bottom=.19, wspace=.31, hspace=.42)
    for ext in ('png', 'pdf'):
        p = OUT/'figures'/f'{stem}.{ext}'
        fig.savefig(p, dpi=200)
        EXPORTS.append(p)
    plt.close(fig)
    print('Exported', stem, flush=True)


def title(fig, text):
    fig.suptitle(text, x=.06, ha='left', y=.98, fontsize=18, fontweight='bold')


def main():
    (OUT/'figures').mkdir(parents=True, exist_ok=True)
    (OUT/'data').mkdir(exist_ok=True)
    font = Path('C:/Windows/Fonts/msyh.ttc')
    font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({'font.family': font_manager.FontProperties(fname=str(font)).get_name(),
        'font.size': 11, 'axes.unicode_minus': False, 'axes.spines.top': False,
        'axes.spines.right': False, 'pdf.fonttype': 42, 'axes.axisbelow': True,
        'text.color': '#20364a', 'axes.labelcolor': '#20364a', 'savefig.facecolor': 'white'})
    colors = {4: '#2369a5', 8: '#bc6023'}
    thermal = read(B/'electrothermal_audit.json')
    local = read(B/'joint_fields/summary.json')
    digest = read(B/'result_digest.json')
    iso = read(I/'completion_readout.json')
    assert thermal['status'] == local['status'] == iso['status'] == 'pass'
    assert len(thermal['thermal']) == len(local['points']) == 62
    assert all(t['qualification_pass'] for t in thermal['thermal'])
    hashes = {str(Path(p).resolve()): h for p,h in {**thermal['sources_sha256'], **local['sources_sha256']}.items()}
    def checked(p):
        p = Path(p).resolve()
        return read(p, hashes[str(p)])
    remote, copied = thermal['candidate_path_map']
    def mapped(p):
        assert p.startswith(remote+'/')
        return Path(copied)/p[len(remote)+1:]
    mesh_path = B/'copied_vm/cases/data/mesh.json'
    mesh = checked(mesh_path)
    xy = np.array([[p['x'],p['y']] for p in sorted(mesh['nodes'],key=lambda x:x['id'])])
    assert len(xy) == 10241
    cells = np.array([t['node_ids'] for t in mesh['triangles']])
    regions = np.array([t['region_id'] for t in mesh['triangles']])
    tri = Triangulation(xy[:,0],xy[:,1],cells)
    si_tri = Triangulation(xy[:,0],xy[:,1],cells,mask=regions != 0)
    edge_regions = {}
    for ids,r in zip(cells,regions):
        for a,b in zip(ids,np.roll(ids,-1)):
            edge_regions.setdefault(tuple(sorted((a,b))),[]).append(r)
    outline = [xy[list(e)] for e,r in edge_regions.items() if len(r)==1 or len(set(r))>1]
    def outline_on(ax):
        ax.add_collection(LineCollection(outline,colors='#64748b',linewidths=.45))
        ax.set(xlabel='原网格 x (μm)',ylabel='原网格 y (μm)')
        ax.set_aspect('equal')
    curves = {}
    curve_rows = []
    states40 = {}
    for g in (4,8):
        ledger = checked(B/f'audit_view/candidate_vg{g}/ledger.json')
        pts = ledger['exact_points']
        assert ledger['status']=='complete' and len(pts)==31
        cfg = checked(mapped(ledger['runs'][0]['directory'])/'input.json')
        native = read_curve(track(N/'normalized'/f'IdVd_Vg{1 if g==4 else 2}_n4_des_drain_curve.csv'))
        bias = np.array([p['bias_V'] for p in pts])
        assert np.allclose(bias,[p[0] for p in native],rtol=0,atol=1e-9)
        current=[]
        for p in pts:
            s=checked(mapped(p['result']))
            current.append(next(c['total_outflow_A_per_m'] for c in s['contacts'] if c['contact']=='drain'))
        states40[g]=s
        t=sorted([p for p in thermal['thermal'] if p['gate']==g],key=lambda p:p['bias_V'])
        c=dict(bias=bias,current=np.array(current),reference=np.array([p[1] for p in native])*1e6,
            peak=np.array([p['candidate_peak_K'] for p in t]),native_peak=np.array([p['reference_peak_K'] for p in t]))
        c['error']=100*np.abs(c['current'][1:]/c['reference'][1:]-1)
        assert np.isclose(c['error'].max(),thermal['metrics'][f'Vg{g}']['relative_error_percent']['max'],atol=1e-9)
        curves[g]=c
        for i,v in enumerate(bias):
            curve_rows.append(dict(gate_V=g,drain_V=v,vela_A_per_m=current[i],sentaurus_A_per_m=c['reference'][i],
                vela_peak_K=c['peak'][i],sentaurus_peak_K=c['native_peak'][i],relative_error_percent=None if i==0 else c['error'][i-1]))
    table('d0_curves',curve_rows)
    # Exact physical geometry, not a schematic reconstruction.
    fig,axs=plt.subplots(1,2,figsize=(13,6.5));title(fig,'原始 LDMOS 网格、材料与掺杂分布')
    axs[0].tripcolor(tri,facecolors=regions,cmap=matplotlib.colors.ListedColormap(['#b8cee0','#eed7a5','#eed7a5']),vmin=0,vmax=2,shading='flat')
    axs[0].triplot(tri,lw=.12,color='#526170',alpha=.38)
    for contact in mesh['contacts']:
        ids=contact['node_ids'];pos=xy[ids]
        axs[0].plot(pos[:,0],pos[:,1],'.',ms=2,label=contact['name'])
    axs[0].legend(fontsize=8,loc='upper left',ncol=2)
    net=(np.array(cfg['donors_m3'])-np.array(cfg['acceptors_m3']))/1e6
    lim=max(abs(net.min()),abs(net.max()))
    pc=axs[1].tripcolor(si_tri,net,cmap='RdBu_r',norm=SymLogNorm(linthresh=1e15,vmin=-lim,vmax=lim),shading='gouraud')
    cb=fig.colorbar(pc,ax=axs[1],fraction=.05,pad=.03,
        ticks=[-1e20,-1e18,-1e16,0,1e16,1e18,1e20])
    cb.ax.set_title(r'cm$^{-3}$',fontsize=10,pad=10)
    for ax in axs:outline_on(ax)
    axs[0].set_title('10,241 节点 · 19,782 三角形');axs[1].set_title('硅区净掺杂；对称对数色标')
    finish(fig,'01_mesh','原始坐标未旋转；蓝灰色为硅，浅黄色为氧化层。保留输入中的接触附近高掺杂，不将峰值当作沟道均匀掺杂。')
    fig,axs=plt.subplots(1,2,figsize=(13,5.8));title(fig,'R11 电热完整曲线：电流与自热温升均与参考接近')
    for g,c in curves.items():
        for ax,v,r in [(axs[0],c['current'],c['reference']),(axs[1],c['peak']-300,c['native_peak']-300)]:
            ax.plot(c['bias'],r,color=colors[g],lw=2,label=f'Sentaurus · Vg={g} V')
            ax.plot(c['bias'],v,'o',mfc='white',ms=4,color=colors[g],label=f'Vela R11 · Vg={g} V')
            ax.set(xlim=(0,40),xlabel='漏极电压 Vd (V)',ylim=(0,None));ax.grid(alpha=.18)
    axs[0].set(ylabel='单位宽度 Id (A/m)',title='D0 漏极电流');axs[1].set(ylabel='Tmax − 300 K (K)',title='全域峰值温升')
    axs[0].set_ylim(0,1.1*max(max(c['current'].max(),c['reference'].max()) for c in curves.values()))
    axs[1].set_ylim(0,1.1*max(max(c['peak'].max(),c['native_peak'].max())-300 for c in curves.values()))
    for c in curves.values():
        assert max(c['current'].max(),c['reference'].max()) < axs[0].get_ylim()[1]
        assert max(c['peak'].max(),c['native_peak'].max())-300 < axs[1].get_ylim()[1]
    axs[0].legend(fontsize=8.5,frameon=False)
    finish(fig,'02_d0_curves','每栅压 31 个精确点；线为参考、空心点为 R11，无拟合。重合程度较高，误差另见图 03。')
    fig,axs=plt.subplots(1,2,figsize=(13,5.5));title(fig,'R11 全曲线误差：电流小于 0.025%，峰温小于 0.042 K')
    for g,c in curves.items():
        axs[0].plot(c['bias'][1:],c['error'],'-o',ms=3,color=colors[g],label=f'Vg={g} V')
        axs[1].plot(c['bias'],abs(c['peak']-c['native_peak']),'-o',ms=3,color=colors[g],label=f'Vg={g} V')
    for ax in axs:ax.grid(alpha=.18);ax.set(xlabel='Vd (V)',xlim=(0,40),ylim=(0,None));ax.legend(frameon=False)
    axs[0].set(ylabel='漏极电流相对误差 (%)');axs[1].set(ylabel='峰值温度绝对误差 (K)')
    finish(fig,'03_d0_error','电流相对误差排除零电流点；零漏压单独检查平衡。曲线门限采用原合同，未新增“逐点 0.025%”验收规则。')
    # Isothermal curves use independent stage-matched native data.
    fig,axs=plt.subplots(1,2,figsize=(13,5.6));title(fig,'最终 R11 源码：D5 与 D4 等温双栅压完整曲线通过')
    iso_rows=[]
    for ax,stage in zip(axs,['D5','D4']):
        for g in (4,8):
            cp=I/f'{stage.lower()}_vg{g}/score/curve.csv'
            rp=ROOT/'reference_staging/templates_ldmos_stage4_idvd_ablation_20260903/analysis_r2/normalized'/('D5-no-IALMob' if stage=='D5' else 'D4-classical')/f'IdVd_Vg{1 if g==4 else 2}_n4_des_drain_curve.csv'
            v=read_curve(track(cp));r=read_curve(track(rp));assert len(v)==len(r)==31
            assert np.allclose([p[0] for p in v],[p[0] for p in r],rtol=0,atol=1e-9)
            ax.plot([p[0] for p in r],[p[1]*1e6 for p in r],color=colors[g],lw=2,label=f'Sentaurus · {g} V')
            ax.plot([p[0] for p in v],[p[1]*1e6 for p in v],'o',ms=3.5,mfc='white',color=colors[g],label=f'Vela · {g} V')
            for (x,y),(_,z) in zip(v,r):iso_rows.append(dict(stage=stage,gate_V=g,drain_V=x,vela_A_per_m=y*1e6,sentaurus_A_per_m=z*1e6))
        ax.set(title=stage+('：无 IALMob' if stage=='D5' else '：含 IALMob'),xlabel='Vd (V)',ylabel='Id (A/m)',xlim=(0,40),ylim=(0,None));ax.grid(alpha=.18);ax.legend(fontsize=8,frameon=False)
    table('isothermal_idvd',iso_rows)
    finish(fig,'04_isothermal_idvd','300 K 等温；各 62 点。与 D0 的物理配置和控制器不同，不把 D4→D0 的差值单独归因于自热。')
    q=read(I/'g3/qualification.json');v=read_curve(track(q['candidate']));r=read_curve(track(q['reference']))
    assert len(v)==len(r)==31
    table('g3_idvg',[dict(gate_V=x,vela_A_per_m=y*1e6,sentaurus_A_per_m=z*1e6) for (x,y),(_,z) in zip(v,r)])
    fig,axs=plt.subplots(1,2,figsize=(13,5.5));title(fig,'G3 转移特性：阈值误差 10.75 mV，最大跨导误差 1.61%')
    for ax in axs:
        ax.plot([p[0] for p in r],[p[1]*1e6 for p in r],lw=2,color='#2369a5',label='Sentaurus')
        ax.plot([p[0] for p in v],[p[1]*1e6 for p in v],'o',ms=4,mfc='white',color='#bc6023',label='Vela / R11 源码')
        ax.set(xlabel='Vg (V)',ylabel='Id (A/m)',xlim=(0,5));ax.grid(alpha=.18);ax.legend(frameon=False)
    axs[0].set_yscale('log');axs[0].set_title('对数坐标：亚阈值区');axs[1].set_title('线性坐标：强反型区');axs[1].set_ylim(bottom=0)
    finish(fig,'05_g3','Vd=0.1 V，Vg=0–5 V，31 点；G3 为无 IALMob 等温控制。亚阈值按原对数误差门限评估，不能套用 D0 的电流误差。')
    nm=read(N/'manifest.json');field=next(f for f in nm['fields'] if f['gate_V']==8 and f['point_index']==30)
    nt=read(field['temperature_file'],field['temperature_sha256'])
    ref=np.array(nt['temperature_K']);val=np.array(states40[8]['temperature_K']);err=val-ref
    assert nt['node_id']==list(range(len(xy)))
    fig,axs=plt.subplots(1,3,figsize=(14.5,5.6));title(fig,'Vg=8 V / Vd=40 V：温度场及同节点误差')
    limit=float(np.max(np.abs(err)))
    for ax,arr,label in zip(axs,[ref,val,err],['Sentaurus (K)','Vela R11 (K)','Vela − Sentaurus (K)']):
        iserr=ax==axs[2]
        pc=ax.tripcolor(tri,arr,shading='gouraud',cmap='RdBu_r' if iserr else 'inferno',norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit) if iserr else matplotlib.colors.Normalize(300,520))
        outline_on(ax);ax.set_title(label);fig.colorbar(pc,ax=ax,fraction=.055,pad=.025)
    finish(fig,'06_temperature_fields',f'两幅温度图共用 300–520 K 色标；差值未截断，最大同节点绝对温差 {limit:.5f} K。原始全域网格，无空间重采样。')
    fig,axs=plt.subplots(1,2,figsize=(13,5.6));title(fig,'关键局部量：全偏压密度与带边通过已批准门限')
    local_rows=[]
    for g in (4,8):
        points=[p for p in local['points'] if p['gate_V']==g]
        for carrier,ls,label in [('electrons_m3','-','电子'),('holes_m3','--','空穴')]:
            yy=[p['carriers'][carrier]['silicon_oxide_interface']['weighted_relative_rms']*100 for p in points]
            axs[0].plot([p['bias_V'] for p in points],yy,ls,color=colors[g],label=f'{g} V · 界面{label}')
        for band,ls,label in [('conduction_band_eV','-','导带'),('valence_band_eV','--','价带')]:
            axs[1].plot([p['bias_V'] for p in points],[p['bands'][band]['maximum_absolute_error_eV']*1000 for p in points],ls,color=colors[g],label=f'{g} V · {label}')
        for p in points:
            local_rows.append(dict(gate_V=g,drain_V=p['bias_V'],**{f'{c}_{group}_rms_percent':p['carriers'][c][group]['weighted_relative_rms']*100 for c in ['electrons_m3','holes_m3'] for group in ['all_silicon','silicon_oxide_interface']},**{k+'_max_meV':v['maximum_absolute_error_eV']*1000 for k,v in p['bands'].items()}))
    table('local_field_metrics',local_rows)
    for ax,lim,ylabel in [(axs[0],2,'界面密度体积加权相对 RMS (%)'),(axs[1],30,'全硅带边最大绝对差 (meV)')]:
        ax.axhline(lim,color='#667085',ls=':',label='批准门限');ax.set(xlabel='Vd (V)',ylabel=ylabel,xlim=(0,40),ylim=(0,lim*1.12));ax.grid(alpha=.18);ax.legend(fontsize=8,frameon=False,ncol=2)
    finish(fig,'07_local_metrics','密度仅统计参考值 > 各自全硅峰值 × 10^(-12) 的正体积节点；2% 为分组加权 RMS，并非逐节点误差上限。')
    # Paired local carrier fields from the verified native silicon ordering.
    mapping=next(Path(p) for p in nm['sources_sha256'] if Path(p).name=='LatticeTemperature_region0.csv')
    mapping=track(mapping,nm['sources_sha256'][str(mapping)])
    with mapping.open(encoding='utf-8') as f: ids=np.array([int(x['node_id']) for x in csv.DictReader(f)])
    tdr=ROOT/'reference_staging/templates_ldmos_d0_electrothermal_20260912/native_full_r1/raw/field_vg8_0030_des.tdr'
    track(tdr,field['tdr_sha256'])
    loc=checked(B/'joint_fields/vg8_30/output.json')['results'];assert np.array_equal(ids,[x['id'] for x in loc])
    with h5py.File(tdr,'r') as f:
        assert geometry_hash(f)==nm['geometry_sha256']
        native={k:values(f,name,0)*1e6 for k,name in [('electrons_m3','eDensity'),('holes_m3','hDensity')]}
    fig,axs=plt.subplots(1,2,figsize=(12.5,5.8));title(fig,'高压载流子场逐节点对照：同一硅区、同一偏压')
    density_rows=[]
    for ax,k,label in zip(axs,['electrons_m3','holes_m3'],['电子 n','空穴 p']):
        r=native[k];v=np.array([p[k]['value'] for p in loc]);active=(r>r.max()*1e-12)&(np.array(cfg['silicon_area_m2'])[ids]>0)
        xx=r[active]/1e6;yy=v[active]/1e6
        ax.scatter(xx,yy,s=5,alpha=.4,color='#2369a5',rasterized=True)
        lo=min(xx.min(),yy.min());hi=max(xx.max(),yy.max());ax.plot([lo,hi],[lo,hi],'--',lw=1,color='#bc6023',label='逐节点相等线')
        ax.set(xscale='log',yscale='log',xlabel=r'Sentaurus (cm$^{-3}$)',ylabel=r'Vela R11 (cm$^{-3}$)',title=f'{label} · {active.sum()} 个筛选节点');ax.grid(alpha=.14);ax.legend(fontsize=9,frameon=False)
        for node,a,b,used in zip(ids,r,v,active):density_rows.append(dict(node_id=int(node),carrier=k,reference_m3=a,candidate_m3=b,in_rms_gate=bool(used)))
    table('density_vg8_vd40',density_rows)
    finish(fig,'08_carrier_parity','Vg=8 V、Vd=40 V；显示所有满足原密度筛选的硅节点，不抽样。对数散点用于看分布，精度结论以图 07 的加权指标为准。')
    fig,axs=plt.subplots(1,2,figsize=(13.5,5.8));title(fig,'同 VM 完整配对：R11 墙钟合格，CPU 成本仍高于原生')
    pairs=[(4,0),(4,1),(8,0),(8,1)];xx=np.arange(4);width=.24
    for j,(variant,label,color) in enumerate([('baseline','R10','#8898ab'),('high','R11','#2369a5'),('native','Sentaurus','#bc6023')]):
        rows=[next(r for r in digest['runs'] if r['gate_V']==g and r['repeat']==rep and r['variant']==variant) for g,rep in pairs]
        for ax,key in [(axs[0],'external_wall_seconds'),(axs[1],'child_cpu_seconds')]:
            bars=ax.bar(xx+(j-1)*width,[r[key] for r in rows],width,label=label,color=color)
            ax.bar_label(bars,fmt='%.0f',fontsize=8,padding=3)
    for ax,label in zip(axs,['外部墙钟 (s)','进程 CPU (s)']):
        ax.set_xticks(xx,[f'Vg={g} V\n第{rep+1}轮' for g,rep in pairs]);ax.set(ylabel=label,ylim=(0,540));ax.grid(axis='y',alpha=.16);ax.legend(frameon=False,fontsize=9,ncol=3,loc='upper left')
    table('paired_performance',[{k:r[k] for k in ['gate_V','repeat','variant','external_wall_seconds','child_cpu_seconds']} for r in digest['runs']])
    finish(fig,'09_performance','每组独立完整配对，含初始化与输出；墙钟包含原生启动等待。不能据此宣称 Vela 求解内核更快，也不拼接最快轮次。')
    fig,axs=plt.subplots(1,2,figsize=(13,5.5));title(fig,'性能收益来源：接触长尾修复减少更新，显式偏导减少装配成本')
    counts=np.array([[366,321],[344,302],[344,303],[215,200]])
    for j,g in enumerate((4,8)):
        b=axs[0].bar(np.arange(4)+(j-.5)*.34,counts[:,j],.34,color=colors[g],label=f'Vg={g} V');axs[0].bar_label(b,padding=3,fontsize=9)
    axs[0].set_xticks(np.arange(4),['R7','R9/R10','R11','Sentaurus']);axs[0].set(ylabel='漏压阶段 Newton 更新次数',ylim=(0,430));axs[0].legend(frameon=False)
    for j,variant in enumerate(['baseline','high']):
        rows=[next(r for r in digest['runs'] if r['gate_V']==g and r['repeat']==rep and r['variant']==variant) for g,rep in pairs]
        b=axs[1].bar(xx+(j-.5)*.32,[r['performance']['assembly_seconds'] for r in rows],.32,color=['#8898ab','#2369a5'][j],label=['R10','R11'][j]);axs[1].bar_label(b,fmt='%.0f',padding=3,fontsize=8)
    axs[1].set_xticks(xx,[f'{g} V / 轮{rep+1}' for g,rep in pairs]);axs[1].set(ylabel='装配累计时间 (s，含初始化)',ylim=(0,335));axs[1].legend(frameon=False)
    for ax in axs:ax.grid(axis='y',alpha=.16)
    track(ROOT/'docs/validation/templates_ldmos_contact_sweep_2026-09-16.md')
    table('newton_counts',[dict(version=version,gate_V=g,updates=int(counts[i,j])) for i,version in enumerate(['R7','R9/R10','R11','Sentaurus']) for j,g in enumerate((4,8))])
    finish(fig,'10_newton_assembly','左：各冻结阶段证据，排除栅压预偏置，原生计正因子更新；不是单次共同实验。右：同轮 R10/R11 对照，装配含初始化。')
    # Rendered algorithm diagram; Mermaid source is also retained in report.md.
    fig,ax=plt.subplots(figsize=(11,10));title(fig,'生产电热求解流程：合格状态推进与失败回退')
    ax.set_xlim(0,10);ax.set_ylim(0,12);ax.axis('off')
    boxes=[(5,11,'冻结网格、物理参数、门限与程序哈希'),(5,9.6,'复用不可变几何准备\n中性 300 K → Poisson 预偏置 → 四方程平衡'),(5,8.0,'选择下个偏压：自适应步长 + 精确点截短\n仅用合格历史构造预测初值'),(5,6.4,'装配残差与耦合 Jacobian\n复用符号分析 → 数值分解 → Newton 方向'),(5,4.8,'变量限幅、线搜索、准费米重表示\n近稳态有限接触代数一致性检查'),(5,3.2,'原块 / 逐行 / KCL / 热平衡门限通过？'),(2.4,1.5,'通过：保存合格状态与检查点\n更新历史；继续至 40 V'),(7.6,1.5,'未通过：继续迭代\n或回退合格状态并减步')]
    for x,y,s in boxes:
        w=8 if y>2 else 4.4
        ax.add_patch(FancyBboxPatch((x-w/2,y-.48),w,.96,boxstyle='round,pad=0.1',facecolor='#eaf1f8',edgecolor='#507797',lw=1.2))
        ax.text(x,y,s,ha='center',va='center',fontsize=11)
    for ya,yb in [(11,9.6),(9.6,8),(8,6.4),(6.4,4.8),(4.8,3.2)]:ax.annotate('',(5,yb+.58),(5,ya-.58),arrowprops=dict(arrowstyle='->',color='#507797',lw=1.5))
    for x in [2.4,7.6]:ax.annotate('',(x,2.08),(5,2.62),arrowprops=dict(arrowstyle='->',color='#507797',lw=1.5))
    finish(fig,'11_solver_flow','流程是当前生产求解器的概念级说明；不声称复刻 Sentaurus 未公开的内部实现。终态验收门限独立于加速策略。')
    # Small tables used verbatim by the meeting report.
    summary=dict(period='2026-09-12 to 2026-09-17',d0_metrics=thermal['metrics'],
        thermal_40V={str(g):next(t for t in thermal['thermal'] if t['gate']==g and abs(t['bias_V']-40)<1e-9) for g in (4,8)},
        d0_total_points=62,isothermal_total_points=sum(x['points'] for x in iso['curves']),
        temperature_max_same_node_error_K=limit,mesh_nodes=len(xy),mesh_triangles=len(cells),
        performance_pairs=digest['pairs'],isothermal_curves=iso['curves'])
    (OUT/'data/summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
    track(Path(__file__))
    (OUT/'data/figure_provenance.json').write_text(json.dumps(dict(sources_sha256=USED,
        outputs_sha256={str(p.relative_to(OUT)):sha(p) for p in EXPORTS},
        scope='Frozen real evidence only; no new simulation, interpolation of bias points, smoothing or fitted reference adjustment.'),indent=2,ensure_ascii=False),encoding='utf-8')
    print('PASS: source hashes, point alignment, curve metrics and field mapping; 11 figures exported')


if __name__=='__main__':
    main()
