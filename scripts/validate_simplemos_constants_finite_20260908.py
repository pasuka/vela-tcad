"""Finite manual-constant compatibility control, four points/two seeds.

No global constant or production default changes. Reuse the coherent 47-TU
diagnostic binary qualified on the three independent parameter directions.
"""
import argparse,copy,math
from decimal import Decimal
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import validate_simplemos_constants_srh_calibration_20260907 as old

a=old.a;d=old.d;v=old.v;p=old.p;s=old.s
LOCAL=p.REPO/'build-release/simplemos_constants_srh_finite_20260908/constants'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/constants'
RUNNER=old.RUNNER
NATIVE=dict(q=1.602192e-19,kb=1.380662e-23,eps0=8.8542e-12)
BASE=dict(q=1.602176634e-19,kb=1.380649e-23,eps0=8.8541878128e-12)

def env(native=True):
    e=v.environment()
    if native:
        for key in NATIVE:e['VELA_CANDIDATE_'+key.upper()+'_RELATIVE']=format(NATIVE[key]/BASE[key]-1,'.17g')
    return e
def execute(path,native=True):return v.execute(path,RUNNER,env(native))
def cases():return old.cases()

def prepare():
    a.verify(old.OUT/'final_evidence.json');a.verify(old.OUT/'validation_evidence.json');a.verify(old.OUT/'build_evidence.json')
    assert all(x['qualified']=='True' for x in a.rows(old.OUT/'response.csv'))
    cs=cases();files=[Path(__file__).resolve(),old.OUT/'final_evidence.json',old.OUT/'validation_evidence.json',old.OUT/'build_evidence.json',RUNNER]
    for c in cs:
        c['finite_jobs']=[];seed=s.LOCAL/c['key']/'generated/state.csv'
        for label,native,state in [('baseline',False,seed),('manual_constants',True,seed),('manual_constants_native_seed',True,v.m.previous.LOCAL/c['key']/'native_referenced_joint/state.csv')]:
            dest=LOCAL/c['key']/label;cfg=old.cfg(c);cfg.update(state_file=str(state),output_state_file=str(dest/'state.csv'))
            a.write(dest/'config.json',cfg);v.post_config(cfg,dest);files += [dest/(name+'.json') for name in ('config','all_row','acceptance_edges')]+v.probes(cfg,dest/'post',dest/'state.csv')+[state]
            c['finite_jobs'].append(dict(label=label,native=native,dest=str(dest)))
    a.write(OUT/'contract.json',dict(cases=cs,DC=12,manual_constants=NATIVE,production_constants=BASE,
        candidate='Simultaneous q/kb/eps0 compatibility convention from manual, not fitted values; coherent constants in all 47 core TUs.',
        qualification='All original 1814 row/KCL/global-source gates; dual phi<=1e-6 V, density relative<=1e-4, Id relative<=1e-6.',
        physics_gate='High-NWell absolute Id improves, low-NWell absolute nonworse, paired log error improves at both Vd; full-field diagnostics retain all 907 free Si nodes.',
        untouched='Generated box, element-box Masetti, signed Si Poisson charge, original SRH source volume, material inputs, biases and solver settings.',
        calibrated_scope='Vg=.8, n19/n23 x Vd=.05/1; no extension to Vg=1 or full curves.',native_new_runs=False))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json']);print('Frozen 12 finite-constant DC',flush=True)

def run():
    a.verify(OUT/'freeze.json')
    def one(c):
        rows=[]
        for job in c['finite_jobs']:
            dest=Path(job['dest']);status=execute(dest/'config.json',job['native'])
            r=dict(key=c['key'],label=job['label'],qualified=False,failure=status.get('failure_reason',''))
            if (dest/'state.csv').exists():
                for name in ('all_row','acceptance_edges'):execute(dest/(name+'.json'),job['native'])
                r.update(s.prior.old.qualify(c,dest))
                for name in ('functional','edges','terms','mobility'):execute(dest/'post'/(name+'.json'),job['native'])
            rows.append(r);print(c['key'],job['label'],r['qualified'],r.get('max_row_ratio'),flush=True)
        return rows
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(one,a.read(OUT/'contract.json')['cases']) for r in group]
    keys=list(dict.fromkeys(k for r in rows for k in r));a.write_csv(OUT/'dc.csv',[{k:r.get(k,'') for k in keys} for r in rows])

def native_fields(c):
    result={};root=s.prior.old.l.raw_root(c)/'fields'
    for field,name in [('psi','ElectrostaticPotential'),('phin','eQuasiFermiPotential'),('phip','hQuasiFermiPotential'),('electrons_m3','eDensity'),('holes_m3','hDensity'),('SRH','srhRecombination')]:
        result[field]={int(x['node_id']):Decimal(x['component0']) for x in a.rows(root/(name+'_region0.csv'))}
    return result

def field_metrics(c,label,dest,volume_ratios=None,nodes=(792,1057,1000,1009)):
    geo,mask=v.m.previous.prior.support(c);ids=np.flatnonzero(mask);native=native_fields(c)
    state=d.ordered(dest/'state.csv',geo.count);terms=d.ordered(dest/'all_row.csv',geo.count);edges=a.rows(dest/'acceptance_edges.csv')
    ratio=next(float(e['electron_particle_line_flux_per_m_s'])/float(e['electron_flux']) for e in edges if float(e['electron_flux'])!=0)
    metrics=[];local=[];values={}
    for field in ('psi','phin','phip'):
        diff=np.array([float(v.m.previous.prior.physical(state[i],field)-native[field][i]) for i in ids]);values[field]=diff
        metrics.append(dict(key=c['key'],label=label,field=field,units='V',max_absolute=float(max(abs(diff))),rms=float(np.linalg.norm(diff)/np.sqrt(len(ids))),worst_node=int(ids[np.argmax(abs(diff))])))
    for field in ('electrons_m3','holes_m3'):
        ref=np.array([float(native[field][i])*1e6 for i in ids]);actual=np.array([float(state[i][field]) for i in ids]);diff=actual/ref-1;values[field]=diff
        metrics.append(dict(key=c['key'],label=label,field=field,units='relative',max_absolute=float(max(abs(diff))),rms=float(np.linalg.norm(diff)/np.sqrt(len(ids))),worst_node=int(ids[np.argmax(abs(diff))])))
    rates=np.array([float(terms[i]['electron_recombination'])*ratio/(geo.volumes['all_cell'][i]*(volume_ratios or {}).get(int(i),1.))/1e6 for i in ids]);ref=np.array([float(native['SRH'][i]) for i in ids]);diff=rates-ref
    metrics.append(dict(key=c['key'],label=label,field='SRH',units='cm^-3 s^-1',max_absolute=float(max(abs(diff))),rms=float(np.linalg.norm(diff)/np.sqrt(len(ids))),worst_node=int(ids[np.argmax(abs(diff))])))
    for k,i in enumerate(ids):
        if i not in nodes:continue
        local.append(dict(key=c['key'],label=label,node_id=int(i),x_um=geo.coords[i][0],y_um=geo.coords[i][1],
                          **{f+'_delta_V':float(values[f][k]) for f in ('psi','phin','phip')},electron_density_relative=float(values['electrons_m3'][k]),hole_density_relative=float(values['holes_m3'][k]),SRH_cm3_s=float(rates[k]),native_SRH_cm3_s=float(ref[k]),SRH_relative=float(diff[k]/ref[k]) if ref[k] else 'zero_native',SRH_volume_ratio=(volume_ratios or {}).get(int(i),1.)))
    return metrics,local

def analyze():
    a.verify(OUT/'freeze.json');dc=a.rows(OUT/'dc.csv');comparisons=[];dual=[];fields=[];locals=[];ports=[]
    for c in cases():
        root=LOCAL/c['key'];geo,mask=v.m.previous.prior.support(c);rows={r['label']:r for r in dc if r['key']==c['key']}
        for label,r in rows.items():
            Id=float(r['current_A_per_um']);comparisons.append(dict(key=c['key'],device=c['device'],vd=c['vd'],vg=c['vg'],label=label,current_A_per_um=Id,native_current_A_per_um=c['native_Id_A_per_um'],relative_error=Id/c['native_Id_A_per_um']-1,qualified=r['qualified']=='True'))
            met,loc=field_metrics(c,label,root/label);fields+=met;locals+=loc
            status=a.read(root/label/'post/functional.status.json');error=max(abs(status['current_A_per_um']/status['contact_current_extractor_A_per_um']-1),abs(Id/status['current_A_per_um']-1));ports.append(dict(key=c['key'],label=label,relative=error,qualified=error<=1e-8))
        states=[d.ordered(root/lab/'state.csv',geo.count) for lab in ('manual_constants','manual_constants_native_seed')]
        delta=v.m.previous.prior.delta_states(*states,mask);ide=abs(float(rows['manual_constants']['current_A_per_um'])/float(rows['manual_constants_native_seed']['current_A_per_um'])-1)
        dual.append(dict(key=c['key'],**delta,Id_relative=ide,qualified=max(delta[k] for k in ('psi_max_V','phin_max_V','phip_max_V'))<=1e-6 and delta['density_max_relative']<=1e-4 and ide<=1e-6))
    pairing=[]
    for vd in (.05,1.):
        rows={(r['device'],r['label']):r for r in comparisons if r['vd']==vd};rel=lambda dev,lab:rows[(dev,lab)]['relative_error']
        high=abs(rel('n23','manual_constants'))<abs(rel('n23','baseline'));low=abs(rel('n19','manual_constants'))<=abs(rel('n19','baseline'))
        pair=abs(math.log((1+rel('n23','manual_constants'))/(1+rel('n19','manual_constants'))))<=abs(math.log((1+rel('n23','baseline'))/(1+rel('n19','baseline'))))
        pairing.append(dict(vd=vd,high_improved=high,low_nonworse=low,pair_improved=pair,qualified=high and low and pair))
    for name,data in [('comparison',comparisons),('dual',dual),('fields',fields),('local_fields',locals),('ports',ports),('physical_gates',pairing)]:a.write_csv(OUT/(name+'.csv'),data)
    summary=dict(DC=len(dc),qualified_DC=sum(x['qualified']=='True' for x in dc),qualified_dual=sum(x['qualified'] for x in dual),qualified_ports=sum(x['qualified'] for x in ports),qualified_physical_gates=sum(x['qualified'] for x in pairing),
        max_row_ratio=max(float(x['max_row_ratio']) for x in dc),max_manual_absolute_Id_error=max(abs(x['relative_error']) for x in comparisons if x['label']=='manual_constants'),
        qualified=all(x['qualified']=='True' for x in dc) and all(x['qualified'] for x in dual+ports+pairing))
    a.write(OUT/'summary.json',summary)
    d.matrix.freeze(OUT/'validation_evidence.json',[OUT/'freeze.json',OUT/'summary.json']+list(OUT.glob('*.csv'))+[x for c in cases() for x in (LOCAL/c['key']).rglob('*') if x.is_file()]);print(summary,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run','analyze'));globals()[parser.parse_args().action]()
