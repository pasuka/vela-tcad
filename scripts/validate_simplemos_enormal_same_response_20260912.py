"""Same Enormal material perturbation in Vela and Sentaurus, with frozen gates."""
import argparse,copy,math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import validate_simplemos_enormal_candidate_20260912 as q
import simplemos_enormal_material_scale_20260912 as n

a,d=q.a,q.d
L=q.R/'build-release/enormal_same_response_20260912'
O=q.R/'reference_tcad/simplemos_sentaurus2022/enormal_same_response_20260912'

def prepare():
    for path in (q.O/'comparison_evidence.json',q.O/'real_column_evidence.json',n.O/'pilot_evidence.json',n.O/'rest_evidence.json'):a.verify(path)
    assert a.read(q.O/'summary.json')['all_qualified'] and a.read(q.O/'real_column_summary.json')['qualified']
    assert a.read(n.O/'pilot_summary.json')['all_qualified'] and a.read(n.O/'rest_summary.json')['all_qualified']
    assert not (O/'contract.json').exists()
    selected=a.rows(q.O/'selected_states.csv');jobs=[];files=[]
    for cc in a.read(q.O/'inputs.json')['cases']:
        if cc['index']!=40:continue
        base=next(r for r in selected if r['key']==cc['key'] and r['arm']=='phumob_seed');assert base['qualified']=='True'
        for label,delta in [('zero',0.),('minus_large',-.001),('minus_small',-.0005),('plus_small',.0005),('plus_large',.001)]:
            cfg=q.config_case(cc);surface=cfg['solver']['mobility']['surface']
            surface['acoustic_factor']=1+delta;surface['roughness_factor']=1+delta
            dest=L/'inputs'/cc['key']/label;template=dest/'template.json';q.v.write_same(template,cfg)
            seed=Path(base['dest'])/'state.csv'
            jobs.append(dict(**{k:v for k,v in cc.items() if k!='template'},template=str(template),seed=str(seed),label=label,delta=delta,base_Id_A_per_um=float(base['current_A_per_um'])))
            files += [template,seed]
    contract=dict(jobs=jobs,dc_gates=q.v.GATES,response_gates=dict(raw_cross_solver_relative=1e-3,two_amplitude_relative=1e-3,even_over_odd=.01,signal_over_drift=100),gate_provenance='Existing same-source port response gates: analyze_simplemos_native_pair_source_20260909.py. No change to DC, row, source or dual-initialization gates.',direction='Direct Vela inverse acoustic/roughness factors f=1+delta match native B,C,delta,eta divided by f.',scope='Four Vg=.8 V cases, two amplitudes and both signs. Zero before signed; no finite fitted mobility replacement.')
    a.write(O/'contract.json',contract)
    files += [Path(__file__).resolve(),O/'contract.json',q.O/'candidate_freeze.json',q.O/'comparison_evidence.json',q.O/'real_column_evidence.json',n.O/'pilot_evidence.json',n.O/'rest_evidence.json',q.R/'scripts/analyze_simplemos_native_pair_source_20260909.py']
    d.matrix.freeze(O/'freeze.json',files)

def run(stage):
    a.verify(O/'freeze.json');q.configure();q.v.LOCAL=L;q.v.OUT=O
    if stage=='signed':
        a.verify(O/'zero_evidence.json');assert a.read(O/'zero_summary.json')['all_qualified']
    jobs=[j for j in a.read(O/'contract.json')['jobs'] if (j['label']=='zero')==(stage=='zero')]
    def one(j):return q.v.solve_target(j,j['index'],j['label'],Path(j['seed']))
    with ThreadPoolExecutor(max_workers=2) as pool:attempts=[r for group in pool.map(one,jobs) for r in group]
    q.v.csv_union(O/f'{stage}_attempts.csv',attempts)
    if stage=='zero':
        checks=[]
        for j in jobs:
            group=[r for r in attempts if r['case']==j['case']];r=next((r for r in group if r['qualified']),group[-1])
            drift=abs(float(r.get('current_A_per_um',math.nan))/j['base_Id_A_per_um']-1)
            checks.append(dict(case=j['case'],qualified=bool(r['qualified']) and drift<=q.v.GATES['initialization_Id_relative'],zero_Id_relative=drift))
        a.write_csv(O/'zero_checks.csv',checks);a.write(O/'zero_summary.json',dict(all_qualified=all(r['qualified'] for r in checks)))
    files=[O/'freeze.json',O/f'{stage}_attempts.csv']+[f for f in (L/'vela').rglob('*') if f.is_file()]
    if stage=='zero':files += [O/'zero_checks.csv',O/'zero_summary.json']
    d.matrix.freeze(O/f'{stage}_evidence.json',files)

def quantum(x):
    return 10.**(math.floor(math.log10(abs(x)))-14) if x else 0.

def analyze():
    for stage in ('zero','signed'):a.verify(O/f'{stage}_evidence.json')
    contract=a.read(O/'contract.json');attempts=a.rows(O/'zero_attempts.csv')+a.rows(O/'signed_attempts.csv')
    native=a.rows(n.O/'pilot_points.csv')+a.rows(n.O/'rest_points.csv');results=[]
    for case in sorted({j['case'] for j in contract['jobs']}):
        chosen={}
        for label in ('zero','minus_large','minus_small','plus_small','plus_large'):
            group=[r for r in attempts if r['case']==case and r['arm']==label]
            chosen[label]=next((r for r in group if r['qualified']=='True'),group[-1])
        original=next(j for j in contract['jobs'] if j['case']==case)
        nr={r['label']:r for r in native if r['case']==case}
        I0=float(chosen['zero']['current_A_per_um']);N0=float(nr['zero']['Id_A_per_um'])
        drift=abs(I0-original['base_Id_A_per_um']);slopes=[];local=[]
        for label,amp in [('large',.001),('small',.0005)]:
            ip=float(chosen['plus_'+label]['current_A_per_um']);im=float(chosen['minus_'+label]['current_A_per_um'])
            np=float(nr['plus_'+label]['Id_A_per_um']);nm=float(nr['minus_'+label]['Id_A_per_um'])
            slope=(ip-im)/(2*amp);native_slope=(np-nm)/(2*amp);slopes.append(slope)
            odd=abs(ip-im)/2;even=abs(math.fsum((ip-I0,im-I0)))/2
            raw=abs(slope/native_slope-1) if native_slope else math.inf
            normalized=(slope/I0)/(native_slope/N0)-1 if native_slope else math.inf
            noise=max(drift,math.ulp(I0));native_noise=max(quantum(N0),(quantum(np)+quantum(nm))/4)
            native_snr=abs(np-nm)/2/native_noise if native_noise else 0.
            good=all(r['qualified']=='True' for r in chosen.values()) and odd>0 and odd>=100*noise and even<=.01*odd and raw<=1e-3 and native_snr>=100
            local.append(dict(case=case,amplitude=amp,Vela_derivative_A_per_um=slope,native_derivative_A_per_um=native_slope,raw_relative_error=raw,normalized_relative_error=normalized,even_over_odd=even/odd if odd else None,signal_over_zero_noise=odd/noise,native_signal_over_export_noise=native_snr,qualified=good))
        linearity=abs(slopes[1]/slopes[0]-1) if slopes[0] else math.inf
        for r in local:r.update(two_amplitude_relative=linearity);r['qualified'] &= linearity<=1e-3
        results += local
    a.write_csv(O/'comparison.csv',results)
    summary=dict(attempts=len(attempts),failed_attempts=sum(r['qualified']!='True' for r in attempts),response_checks=len(results),qualified=sum(r['qualified'] for r in results),max_raw_relative=max(r['raw_relative_error'] for r in results),max_normalized_relative=max(abs(r['normalized_relative_error']) for r in results),all_qualified=all(r['qualified'] for r in results),acceptance_changed=False)
    a.write(O/'summary.json',summary);d.matrix.freeze(O/'comparison_evidence.json',[O/'zero_evidence.json',O/'signed_evidence.json',O/'comparison.csv',O/'summary.json'])
    print(summary,results,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=('prepare','zero','signed','analyze'));args=ap.parse_args()
    if args.action in ('zero','signed'):run(args.action)
    else:globals()[args.action]()
