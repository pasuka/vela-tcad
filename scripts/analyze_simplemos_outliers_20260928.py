"""Analyze frozen results and native PLT files; does not run any solver."""
import csv, hashlib, json, math, shutil, statistics
from pathlib import Path
from sentaurus_import import parse_quoted_list, parse_values_block

ROOT = Path(__file__).resolve().parents[1]
B = ROOT / 'build/outlier_analysis_20260928'
O = B / 'analysis'

def rows(p):
    with p.open(newline='', encoding='utf-8') as f: return list(csv.DictReader(f))

def table(name, data):
    keys = list(dict.fromkeys(k for r in data for k in r))
    with (O/name).open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(data)

def number(r, k): return float(r[k])
def key(r): return r['case'], int(r['index'])
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    O.mkdir(exist_ok=True)
    points=rows(B/'all_points.csv'); outliers=rows(B/'outliers.csv')
    assert len(points)==816 and len({key(r) for r in points})==816
    assert sum(abs(number(r,'Id_error_percent'))>2 for r in points)==6
    for r in points:
        assert r['qualified']=='True'
        assert math.isclose((number(r,'current_A_per_um')/number(r,'native_Id_A_per_um')-1)*100,number(r,'Id_error_percent'),abs_tol=1e-11)
        assert number(r,'max_row_ratio')<=1e-6 and number(r,'kcl_over_Id')<=1e-8
    native={}; nt=[]; provenance=[]
    raw=ROOT/'build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics'
    for case in sorted({r['case'] for r in points}):
        dev=case.split('_')[0]
        p=raw/f'sentaurus_raw/sentaurus_bundle/{dev}/IdVg_{case}_des.plt'
        text=p.read_text(); names=parse_quoted_list(text,'datasets')
        data=[dict(zip(names,x)) for x in parse_values_block(text,len(names))]
        for r in [x for x in points if x['case']==case]:
            match=[x for x in data if abs(x['gate OuterVoltage']-number(r,'vg'))<1e-10]
            assert len(match)==1, (case,r['vg'],len(match))
            x=match[0]; assert x['drain TotalCurrent']==number(r,'native_Id_A_per_um')
            assert abs(x['drain OuterVoltage']-number(r,'vd'))<1e-10
            kcl=math.fsum(x[c+' TotalCurrent'] for c in ['source','drain','substrate','gate'])
            native[key(r)]=dict(x,native_kcl_A_per_um=kcl,native_kcl_over_Id=kcl/abs(x['drain TotalCurrent']))
            for c in ['source','drain','substrate','gate']:
                assert abs(x[c+' TotalCurrent']-x[c+' eCurrent']-x[c+' hCurrent']-x[c+' DisplacementCurrent'])<1e-12*max(abs(x[c+' TotalCurrent']),abs(x[c+' eCurrent']),1e-30)
                nt.append(dict(case=case,index=r['index'],vg=r['vg'],contact=c,
                    electron_A_per_um=x[c+' eCurrent'],hole_A_per_um=x[c+' hCurrent'],
                    displacement_A_per_um=x[c+' DisplacementCurrent'],total_A_per_um=x[c+' TotalCurrent']))
        for src,kind in [(p,'plt'),(raw/f'sentaurus_bundle/{dev}/{case}_des.cmd','deck')]:
            target=O/'native'/src.name; target.parent.mkdir(exist_ok=True);shutil.copy2(src,target)
            provenance.append(dict(source=str(src),file=str(target.relative_to(B)),sha256=sha(src),kind=kind))
    table('native_terminals_all816.csv',nt)
    table('native_kcl_all816.csv',[dict(case=r['case'],index=r['index'],vg=r['vg'],Id_error_percent=r['Id_error_percent'],native_kcl_A_per_um=native[key(r)]['native_kcl_A_per_um'],native_kcl_over_Id=native[key(r)]['native_kcl_over_Id']) for r in points])
    vela=rows(B/'vela_terminals.csv'); vt={(r['case'],int(r['index']),r['contact']):r for r in vela}; comp=[]
    for r in vela:
        x=native[key(r)]; c=r['contact']; row=dict(r)
        for car, label in [('electron','eCurrent'),('hole','hCurrent'),('total','TotalCurrent')]:
            row['native_'+car+'_A_per_um']=x[c+' '+label]
            row['delta_'+car+'_A_per_um']=number(r,car+'_A_per_um')-x[c+' '+label]
        comp.append(row)
    table('terminal_comparison_selected.csv',comp)
    for r in points:
        if (r['case'],int(r['index']),'drain') in vt:
            assert math.isclose(number(vt[r['case'],int(r['index']),'drain'],'total_A_per_um'),number(r,'current_A_per_um'),rel_tol=1e-12)
    compmap={(r['case'],int(r['index']),r['contact']):r for r in comp}
    source={key(r):r for r in rows(B/'srh_source_integrals.csv')}
    params=json.loads((ROOT/'reference_tcad/simplemos_sentaurus2022/simplemos_m8_original_physics_contract_v1.json').read_text())
    devices={r['id']:r for r in params['devices']}
    summary=[]
    for r in outliers:
        case,idx=key(r); x=native[case,idx]; c=compmap[case,idx,'substrate']; s=dict(r)
        s.update({k:v for k,v in devices[r['device']].items() if k not in ['id','nominal']})
        s.update(native_kcl_over_Id=x['native_kcl_over_Id'],native_kcl_percent=x['native_kcl_over_Id']*100,
                 native_substrate_electron_A_per_um=x['substrate eCurrent'],native_substrate_hole_A_per_um=x['substrate hCurrent'],
                 vela_substrate_electron_A_per_um=c['electron_A_per_um'],vela_substrate_hole_A_per_um=c['hole_A_per_um'],
                 substrate_total_delta_A_per_um=c['delta_total_A_per_um'],vela_free_Si_SRH_A_per_um=source[case,idx]['free_Si_SRH_A_per_um'])
        for offset,name in [(-1,'previous'),(1,'next')]:
            neighbor=next(t for t in points if key(t)==(case,idx+offset))
            s[name+'_error_percent']=neighbor['Id_error_percent']
        s['high_vd_error_percent']=next(t['Id_error_percent'] for t in points if key(t)==(r['device']+'_vd_1',idx))
        s['log10_Id_ratio']=math.log10(number(r,'current_A_per_um')/number(r,'native_Id_A_per_um'))
        summary.append(s)
    table('outlier_summary.csv',summary)
    # Audit carrier-resolved source balance; source sign follows collected port convention.
    balance=[]
    for k,s in source.items():
        ele=math.fsum(number(v,'electron_A_per_um') for v in vela if key(v)==k)
        hole=math.fsum(number(v,'hole_A_per_um') for v in vela if key(v)==k)
        srh=number(s,'free_Si_SRH_A_per_um')
        balance.append(dict(case=k[0],index=k[1],electron_terminal_sum_A_per_um=ele,hole_terminal_sum_A_per_um=hole,free_Si_SRH_A_per_um=srh,electron_plus_SRH_A_per_um=ele+srh,hole_minus_SRH_A_per_um=hole-srh))
    table('vela_carrier_source_balance.csv',balance)
    fields=rows(B/'field_statistics.csv'); table('outlier_field_statistics.csv',[r for r in fields if key(r) in {key(s) for s in outliers}])
    stats=dict(points=816,outliers=6,exact_native_matches=816,selected_points=len(source),
        max_abs_high_vd_error_percent=max(abs(number(r,'Id_error_percent')) for r in points if number(r,'vd')==1),
        max_abs_nonoutlier_error_percent=max(abs(number(r,'Id_error_percent')) for r in points if abs(number(r,'Id_error_percent'))<=2),
        native_kcl_over_1e_8_count=sum(abs(x['native_kcl_over_Id'])>1e-8 for x in native.values()),
        native_kcl_over_1percent_count=sum(abs(x['native_kcl_over_Id'])>.01 for x in native.values()),
        max_row_ratio=max(number(r,'max_row_ratio') for r in points),max_kcl_over_Id=max(number(r,'kcl_over_Id') for r in points),
        checks=['816 unique points; all qualified; six >2%','816 exact unique bias PLT matches; zero CSV current difference','native e+h+displacement=total','32 Vela edge-derived drain currents match audit Id'],
        native_spatial_difference_available=False,solver_invoked=False)
    (O/'validation.json').write_text(json.dumps(stats,indent=2)+'\n')
    (O/'native_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    # Six independent axes: zoom emphasizes isolated discrepancies without hiding full curves.
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(12,6.8),layout='constrained')
    for ax,r in zip(axes.flat,outliers):
        data=sorted([x for x in points if x['case']==r['case'] and number(x,'vg')<=.3],key=lambda x:number(x,'vg'))
        ax.plot([number(x,'vg') for x in data],[number(x,'Id_error_percent') for x in data],'o-',label='Id error')
        ax.axhline(2,color='tab:red',linestyle='--',label='+2% gate');ax.axhline(-2,color='tab:red',linestyle='--')
        ax.set(title=f"{r['device']}, Vd = 0.05 V",xlabel='Vg (V)',ylabel='Id error (%)');ax.grid(alpha=.25)
    fig.suptitle('Frozen 816-point run: low-Vd outliers (independent y scales)')
    fig.savefig(O/'outlier_error_zoom.png',dpi=170);plt.close(fig)
    print(json.dumps(stats,indent=2));print('\nSIX-POINT EVIDENCE')
    for r in summary:
        print(r['device'],*[f'{float(r[k]):.7g}' for k in ['vg','Id_error_percent','previous_error_percent','next_error_percent','high_vd_error_percent','native_kcl_percent','native_substrate_electron_A_per_um','vela_substrate_electron_A_per_um','native_substrate_hole_A_per_um','vela_substrate_hole_A_per_um']])

if __name__=='__main__': main()
