"""Recompute solver qualification, centered responses and two-initialization A/B."""
import ast
import copy
from decimal import Decimal
import math
from pathlib import Path
import re
import subprocess
import numpy as np
import validate_simplemos_joint_stable_20260907 as launch

v=launch.v;b=launch.b;p=v.p;a=p.a;d=p.d;OUT=b.OUT;LOCAL=b.LOCAL
REPORT=p.REPO/'docs/validation/simplemos_joint_geometry_self_consistent_2026-09-07.md'


def stripped(cfg):
    cfg=copy.deepcopy(cfg)
    for k in ('state_file','output_state_file'):cfg.pop(k)
    cfg['solver']['local_update_diagnostics'].pop('csv_file');return cfg


def analyze():
    for file in (OUT/'freeze.json',OUT/'preflight_evidence.json',OUT/'numerical_addendum.json'):a.verify(file)
    contract=a.read(OUT/'contract.json');dc=a.rows(OUT/'dc.csv');assert len(dc)==76
    checks=[];fd=[];comparison=[];invariance=[];jvp=[];paired=[];fields=[]
    for c in contract['cases']:
        geo=d.matrix.spatial.m73.Geometry(c['device']);mask=d.matrix.spatial.old.m78.supports(c['device'],geo,.05)[0]['all_si'].copy();mask[geo.contact_nodes]=False
        base=a.read(Path(c['base'])/'config.json');values={};passed={}
        mapped=d.ordered(Path(c['native_initial']),geo.count)
        weights=geo.volumes['signed_si'][mask];weights=weights/np.sum(weights)
        for j in c['jobs']:
            dest=Path(j['config']).parent;cfg=a.read(Path(j['config']));assert stripped(cfg)==stripped(base)
            s=a.read(dest/'config.status.json');q=a.read(dest/'all_row.status.json');terms=d.ordered(dest/'all_row.csv',geo.count)
            post=a.read(dest/'all_row.json');assert post['solver']['global_continuity_closure']==dict(mode='enforce',tolerance=1e-6,source_floor=1e-10)
            ratios=[];zero=0
            for r,keep in zip(terms,mask):
                if not keep:continue
                for car in ('electron','hole'):
                    scale=max(float(r[car+'_flux_abs_sum']),abs(float(r[car+'_recombination'])),abs(float(r[car+'_impact'])));zero+=scale==0
                    ratios.append(abs(float(r[car+'_residual']))/scale if scale else math.inf)
            assert len(ratios)==q['carrier_row_convergence']['qualified_row_count']==1814
            bad=sum(x>1e-6 for x in ratios);assert bad==q['carrier_row_convergence']['violation_count']
            cc=s['contact_currents_A_per_um'];current=cc['drain'];kcl=abs(math.fsum(cc.values()))/abs(current)
            ok=s['exit_code']==q['exit_code']==0 and s['converged'] and bad==zero==0 and q['global_continuity_closure']['satisfied'] and kcl<=1e-8
            reported=next(x for x in dc if x['key']==c['key'] and x['arm']==j['arm'] and x['label']==j['label'])
            assert ok==(reported['qualified']=='True') and current==float(reported['current_A_per_um'])
            values[(j['arm'],j['label'])]=current;passed[(j['arm'],j['label'])]=ok
            checks.append(dict(key=c['key'],arm=j['arm'],label=j['label'],active_rows=len(ratios),violations=bad,zero_scale_rows=zero,max_row_ratio=max(ratios),kcl_over_Id=kcl,qualified=ok))
            linear=a.rows(dest/'linear_summary.csv')
            for iteration in {x['iteration'] for x in linear}:assert [int(x['correction']) for x in linear if x['iteration']==iteration]==list(range(5))
            if j['label'] in ('zero','replacement','replacement_native'):
                status=a.read(dest/'functional.status.json');assert status['exit_code']==0
                porterr=abs(status['current_A_per_um']/status['contact_current_extractor_A_per_um']-1);assert porterr<=1e-8
                jp=launch.check.jvp_metrics(dest/'jvp.csv',c['key'],j['arm']+'/'+j['label'])[0];jvp+=jp
                comparison.append(dict(key=c['key'],device=c['device'],vd=c['vd'],arm=j['arm'],initialization=j['initialization'],
                    current_A_per_um=current,native_Id_A_per_um=c['native_Id_A_per_um'],signed_Id_relative_error=current/c['native_Id_A_per_um']-1,
                    qualified=ok,port_functional_relative=porterr,jvp_changed_blocks_qualified=all(x['qualified'] for x in jp)))
                if j['initialization']=='vela':
                    actual=d.ordered(dest/'state.csv',geo.count)
                    for field in ('psi','phin','phip'):
                        car={'phin':'electron','phip':'hole'}.get(field)
                        def physical(row):
                            if car and car+'_qf_reference_V' in row:return Decimal(row[car+'_qf_reference_V'])+Decimal(row[car+'_qf_increment_V'])
                            return Decimal(row[field])
                        delta=np.array([float(physical(x)-physical(y)) for x,y in zip(actual,mapped)])
                        fields.append(dict(key=c['key'],arm=j['arm'],field=field,qualified_state=ok,free_Si_max_absolute_V=float(max(abs(delta[mask]))),
                            Si_volume_weighted_RMS_V=float(np.sqrt(np.sum(weights*delta[mask]**2))),node_320_delta_V=delta[320],node_324_delta_V=delta[324],node_338_delta_V=delta[338]))
        center=values[('baseline','zero')];drift=abs(center-c['base_Id_A_per_um']);assert abs(center/c['base_Id_A_per_um']-1)<=1e-8
        for arm in ('transport','poisson_volume','joint'):
            pred=float(next(x['unit_prediction_A_per_um'] for x in a.rows(OUT/'predictions.csv') if x['key']==c['key'] and x['arm']==arm))
            odd={n:(values[(arm,'plus_'+n)]-values[(arm,'minus_'+n)])/2 for n in ('full','half')};linearity=abs(odd['full']/(2*odd['half'])-1)
            for name,alpha in [('full',.001),('half',.0005)]:
                plus,minus=values[(arm,'plus_'+name)],values[(arm,'minus_'+name)]
                err=abs(odd[name]/(alpha*pred)-1);even=abs(((plus+minus)/2-center)/odd[name]);snr=abs(odd[name])/max(drift,1e-300);dd=abs(math.log10(center/c['base_Id_A_per_um']))
                signs=(plus-center)*pred>0 and (minus-center)*pred<0
                ok=passed[('baseline','zero')] and all(passed[(arm,s+'_'+name)] for s in ('plus','minus')) and err<=.001 and linearity<=.001 and even<=.01 and snr>=100 and dd<=1e-5 and signs
                fd.append(dict(key=c['key'],arm=arm,amplitude=name,unit_prediction_A_per_um=pred,actual_odd_A_per_um=odd[name],prediction_relative_error=err,two_amplitude_relative=linearity,
                    even_over_odd=even,signal_over_zero_drift=snr,zero_drift_dex=dd,signs_correct=signs,qualified=ok))
            first=d.ordered(LOCAL/c['key']/arm/'replacement/state.csv',geo.count);second=d.ordered(LOCAL/c['key']/arm/'replacement_native/state.csv',geo.count)
            delta={}
            for field in ('psi','phin','phip'):
                car={'phin':'electron','phip':'hole'}.get(field)
                def value(row):
                    if car and car+'_qf_reference_V' in row:return Decimal(row[car+'_qf_reference_V'])+Decimal(row[car+'_qf_increment_V'])
                    return Decimal(row[field])
                delta[field]=max(abs(float(value(x)-value(y))) for x,y,keep in zip(first,second,mask) if keep)
            density=max(abs(float(x[k])/float(y[k])-1) for x,y,keep in zip(first,second,mask) if keep for k in ('electrons_m3','holes_m3'))
            iddiff=abs(values[(arm,'replacement')]/values[(arm,'replacement_native')]-1)
            ok=passed[(arm,'replacement')] and passed[(arm,'replacement_native')] and max(delta.values())<=1e-6 and density<=1e-4 and iddiff<=1e-6
            invariance.append(dict(key=c['key'],arm=arm,psi_max_V=delta['psi'],phin_max_V=delta['phin'],phip_max_V=delta['phip'],density_max_relative=density,Id_relative=iddiff,qualified=ok))
        print('Recomputed all 19 states and responses',c['key'],flush=True)
    promotion=[]
    for arm in ('transport','poisson_volume','joint'):
        allgood=True
        for vd in (.05,1.):
            vals={dev:next(x for x in comparison if x['device']==dev and x['vd']==vd and x['arm']==arm and x['initialization']=='vela') for dev in ('n19','n23')}
            base={dev:next(x for x in comparison if x['device']==dev and x['vd']==vd and x['arm']=='baseline') for dev in ('n19','n23')}
            pair=math.log10(vals['n23']['current_A_per_um']/vals['n19']['current_A_per_um'])-math.log10(vals['n23']['native_Id_A_per_um']/vals['n19']['native_Id_A_per_um'])
            oldpair=math.log10(base['n23']['current_A_per_um']/base['n19']['current_A_per_um'])-math.log10(base['n23']['native_Id_A_per_um']/base['n19']['native_Id_A_per_um'])
            hi=abs(vals['n23']['signed_Id_relative_error'])<abs(base['n23']['signed_Id_relative_error'])
            lo=abs(vals['n19']['signed_Id_relative_error'])<=abs(base['n19']['signed_Id_relative_error'])
            better=abs(pair)<abs(oldpair)
            ok=hi and lo and better and all(x['qualified'] and x['jvp_changed_blocks_qualified'] for x in vals.values());allgood &= ok
            paired.append(dict(arm=arm,vd=vd,pair_error_dex=pair,baseline_pair_error_dex=oldpair,high_NWell_improved=hi,low_NWell_not_worse=lo,pair_improved=better,qualified=ok))
        allgood &= all(x['qualified'] for x in fd if x['arm']==arm) and all(x['qualified'] for x in invariance if x['arm']==arm)
        promotion.append(dict(arm=arm,four_point_gate_passed=bool(allgood),M82_released=False,M83_released=False))
    for name,rows in [('independent_qualification',checks),('response_calibration',fd),('replacement_comparison',comparison),('initialization_invariance',invariance),('replacement_jvp',jvp),('NWell_pairing',paired),('promotion',promotion),('field_comparison',fields)]:a.write_csv(OUT/(name+'.csv'),rows)
    metrics=dict(total_DC=len(checks),qualified_DC=sum(x['qualified'] for x in checks),independently_checked_rows=len(checks)*1814,
        qualified_amplitudes=sum(x['qualified'] for x in fd),total_amplitudes=len(fd),qualified_initialization_pairs=sum(x['qualified'] for x in invariance),
        max_carrier_row_ratio=max(x['max_row_ratio'] for x in checks),max_kcl_over_Id=max(x['kcl_over_Id'] for x in checks),
        max_response_prediction_error=max(x['prediction_relative_error'] for x in fd),max_response_two_amplitude_difference=max(x['two_amplitude_relative'] for x in fd),
        max_initialization_phi_V=max(max(x[k] for k in ('psi_max_V','phin_max_V','phip_max_V')) for x in invariance),
        max_initialization_density_relative=max(x['density_max_relative'] for x in invariance),max_initialization_Id_relative=max(x['Id_relative'] for x in invariance))
    a.write(OUT/'analysis_evidence.json',dict(status='completed_bounded_comparison',metrics=metrics,
        input_hashes={a.rel(f):a.sha(f) for f in [Path(__file__).resolve()]+[x for root in (LOCAL,OUT) for x in root.rglob('*') if x.is_file()]}))
    print(metrics,flush=True)


def seal():
    for f in (OUT/'analysis_evidence.json',OUT/'numerical_addendum.json',OUT/'freeze.json',OUT/'preflight_evidence.json',p.OUT/'validation_evidence.json'):a.verify(f)
    paths=[Path(__file__).resolve(),Path(launch.__file__).resolve(),Path(v.__file__).resolve(),Path(launch.check.__file__).resolve(),Path(b.__file__).resolve(),
        p.REPO/'scripts/test_simplemos_joint_geometry_20260907.py',p.REPO/'scripts/test_simplemos_stable_sg_20260907.py',
        p.REPO/'tests/diagnostics/test_simplemos_native_geometry.cpp',p.REPO/'tests/diagnostics/test_simplemos_stable_sg_derivative.cpp',b.HEADER,b.old.HEADER,REPORT]
    for f in paths:
        if f.suffix=='.py':ast.parse(f.read_text(encoding='utf-8-sig'))
    for root in (LOCAL,OUT,b.old.LOCAL,b.old.OUT):paths += [x for x in root.rglob('*') if x.is_file() and x.name!='validation_evidence.json']
    for target in re.findall(r'\]\(([^)]+)\)',REPORT.read_text(encoding='utf-8')):
        if '://' in target:continue
        f=(REPORT.parent/target.split('#')[0]).resolve()
        if f.name!='validation_evidence.json':assert f.exists(),f;paths.append(f)
    diff=subprocess.run(['git','-c','core.fsmonitor=false','diff','--numstat','--ignore-space-at-eol'],cwd=p.REPO,capture_output=True,text=True)
    assert diff.returncode==0 and diff.stdout.strip()=='127\t1\tsrc/tools/vela_example_runner.cpp'
    assert 'All tests passed' in (b.old.LOCAL/'test_result.log').read_text() and 'All tests passed' in (LOCAL/'test_result.log').read_text()
    a.write(OUT/'validation_evidence.json',dict(status='completed_bounded_joint_geometry_validation',date='2026-09-07',metrics=a.read(OUT/'analysis_evidence.json')['metrics'],
        input_hashes={a.rel(f):a.sha(f) for f in sorted(set(paths))},production_changes=False,acceptance_changes=False,
        isolated_numerical_change='Stable analytical equal-ni Boltzmann SG psi derivatives; original failed preflight retained.',
        limitations=['Four Vg=1 V endpoints only.','No full native mobility, dielectric K, or SRH volume replacement.','No exhaustive arbitrary-column Jacobian qualification.','M82/M83 not released.']))
    a.verify(OUT/'validation_evidence.json');print('Final evidence sealed',flush=True)


if __name__=='__main__':
    import argparse
    q=argparse.ArgumentParser();q.add_argument('action',choices=('analyze','seal'));globals()[q.parse_args().action]()
