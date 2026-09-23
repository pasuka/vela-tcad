"""Audit copied contact-repair sweeps, repeated trajectories and full costs."""
import argparse
import hashlib
import json
from electrothermal_state import read_bound_record
from pathlib import Path
import statistics
import re

from evidence_paths import candidate_path
from analyze_templates_ldmos_predictor_study import state_delta
from run_templates_ldmos_electrothermal_curve import state_gate


def read(p):
    return read_bound_record(p)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def numerical(value):
    """Keep all non-timing output, including rejected contact events."""
    if isinstance(value,dict):
        return {k:numerical(v) for k,v in value.items() if k!='performance' and
                k!='seconds' and not k.endswith('_seconds')}
    if isinstance(value,list):return [numerical(v) for v in value]
    return value


def audit_run(directory, mapping):
    ledger=read(directory/'results/ledger.json')
    if ledger['status']!='complete':raise ValueError('Incomplete run')
    resolve=lambda p:candidate_path(p,mapping)
    totals=dict(updates=0,attempts=len(ledger['runs']),rejected=0,trials=0,
                assemblies=0,factorizations=0,floor_updates=0,contact_attempts=0,
                contact_accepted=0,contact_rejected=0,contact_seconds=0.)
    signatures=[];paths={}
    parent=0.;last=None;recovery_checks=0
    for row in ledger['runs']:
        if row['parent_bias_V']!=parent:raise ValueError('Rejected state leaked into next parent bias')
        path=resolve(row['directory'])/'output.json';data=read(path);paths[str(path)]=sha(path)
        check=state_gate(data,row['bias_V'])
        if check['pass_gate']!=row['gate']['pass_gate']:raise ValueError('Gate ledger mismatch')
        if last is False:recovery_checks+=1
        if check['pass_gate']:parent=row['bias_V']
        last=check['pass_gate']
        totals['updates']+=data['newton_updates'];totals['rejected']+=not last
        totals['trials']+=sum(h['line_search_trials'] for h in data['history'])
        totals['floor_updates']+=sum(h['scaled_l2_before']<1e-9 for h in data['history'])
        for key,outkey in (('assembly_calls','assemblies'),('factorizations','factorizations')):
            totals[outkey]+=data['performance'][key]
        for e in data.get('near_steady_contact_consistency',{}).get('events',[]):
            totals['contact_attempts']+=1;totals['contact_accepted']+=e['accepted']
            totals['contact_rejected']+=not e['accepted'];totals['contact_seconds']+=e['seconds']
        signatures.append(dict(bias=row['bias_V'],parent=row['parent_bias_V'],prediction=numerical(row['prediction']),result=numerical(data)))
    if recovery_checks!=totals['rejected']:raise ValueError('Failed attempt has no following recovery')
    points=[read(resolve(p['result'])) for p in ledger['exact_points']]
    init=[numerical(read(resolve(p['result']))) for p in ledger['initialization_runs']]
    totals['initialization_updates']=sum(p['newton_updates'] for p in ledger['initialization_runs'])
    return dict(totals=totals,signatures=signatures,points=points,init=init,ledger=ledger,hashes=paths)


def audit_native(directory, reference, normalized):
    import h5py
    import numpy as np
    from extract_templates_ldmos_d0_temperature import geometry_hash, values
    from run_templates_ldmos_sentaurus_vm import normalize_plt_files
    from analyze_templates_ldmos_stage4_d5 import read_curve
    gate=int(directory.name.split('_vg')[1].split('_')[0])
    if not normalized.exists():normalize_plt_files(directory,normalized)
    name='IdVd_Vg%d_n4_des_drain_curve.csv'%(1 if gate==4 else 2)
    old=read_curve(reference/'native_fields_r1/normalized'/name);new=read_curve(normalized/name)
    if old!=new or len(new)!=31:raise ValueError('Native curve changed')
    prior=sorted((reference/'native_full_r1/raw').glob('field_vg%d_*_des.tdr'%gate))
    fields=sorted(directory.glob('field_vg%d_*_des.tdr'%gate))
    if len(prior)!=31 or len(fields)!=31:raise ValueError('Native fields missing')
    hashes={}
    for a,b in zip(prior,fields):
        with h5py.File(a,'r') as fa,h5py.File(b,'r') as fb:
            if geometry_hash(fa)!=geometry_hash(fb):raise ValueError('Native geometry changed')
            for name,regions in [('ElectrostaticPotential',(0,1,2)),('eQuasiFermiPotential',(0,)),
                                 ('hQuasiFermiPotential',(0,)),('LatticeTemperature',(0,1,2))]:
                for region in regions:
                    if not np.array_equal(values(fa,name,region),values(fb,name,region)):
                        raise ValueError('Native field changed')
        hashes[str(b)]=sha(b)
    log=(directory/'n4_des.log').read_text(errors='replace');start=log.rfind('Contact drain : 40V')
    if start<0:raise ValueError('Missing native drain stage')
    tables=re.findall(r'Iteration\s+\|Rhs\|[^\n]*\n-+\n(.*?)Finished, because',log[start:],re.S)
    if len(tables)<=30:raise ValueError('Missing native iteration tables')
    counts=[sum(len(p)==8 and p[0].isdigit() and float(p[2])>0
                for p in (line.split() for line in table.splitlines())) for table in tables]
    return dict(curve_and_fields_exact=True,solve_tables=len(tables),
                positive_factor_updates_including_zero=sum(counts),zero_updates=counts[0],hashes=hashes)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix',type=Path,required=True)
    p.add_argument('--path-map',nargs=2,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--first8',type=Path,help='Completed first8 matrix for actual resume/prefix comparison')
    p.add_argument('--native-reference',type=Path,default=Path('reference_staging/templates_ldmos_d0_electrothermal_20260912'))
    args=p.parse_args();summary=read(args.matrix/'summary.json')
    if summary['status']!='complete':raise ValueError('Matrix incomplete')
    if sha(args.matrix/'vela_example_runner')!=summary['runner_sha256']:raise ValueError('Frozen runner changed')
    rows=summary['runs'];data={};report=dict(stage=summary['stage'],status='pass',runs=[],gates={},hashes={})
    variants=('baseline','contact') if summary['stage']=='first8' else ('baseline','contact','native')
    repeats={r['repeat'] for r in rows}
    if repeats!=set(range(len(repeats))) or not repeats:raise ValueError('Missing repeat index')
    if summary['stage']=='full' and len(repeats)<2:raise ValueError('Repeated full timing needs at least two rounds')
    expected={(r,g,v) for r in repeats for g in (4,8) for v in variants}
    observed={(r['repeat'],r['gate_V'],r['variant']) for r in rows}
    if observed!=expected or len(rows)!=len(expected):raise ValueError('Missing or duplicate runs')
    for row in rows:
        directory=candidate_path(row['directory'],args.path_map)
        if row['variant']=='native':
            native=audit_native(directory,args.native_reference,args.output.parent/(directory.name+'_normalized'))
            report['hashes'].update(native.pop('hashes'))
            report['runs'].append(dict(repeat=row['repeat'],gate=row['gate_V'],variant='native',
                external_wall_seconds=row['external_wall_seconds'],child_cpu_seconds=row['child_cpu_seconds'],**native))
            continue
        if sha(directory/'input.json')!=row['input_sha256'] or sha(directory/'results/ledger.json')!=row['ledger_sha256']:
            raise ValueError('Recorded evidence changed')
        value=audit_run(directory,args.path_map);data[row['repeat'],row['gate_V'],row['variant']]=value
        report['hashes'].update(value.pop('hashes'))
        report['runs'].append(dict(repeat=row['repeat'],gate=row['gate_V'],variant=row['variant'],
            costs=value['totals'],external_wall_seconds=row['external_wall_seconds'],child_cpu_seconds=row['child_cpu_seconds']))
    for gate in (4,8):
        repeats=sorted({r['repeat'] for r in rows if r['gate_V']==gate});pairs=[]
        for repeat in repeats:
            base=data[repeat,gate,'baseline'];cand=data[repeat,gate,'contact']
            if base['init']!=cand['init']:raise ValueError('Paired initialization differs')
            if [p['bias_V'] for p in base['ledger']['exact_points']]!=[p['bias_V'] for p in cand['ledger']['exact_points']]:
                raise ValueError('Output voltage changed')
            deltas=[state_delta(a,b) for a,b in zip(base['points'],cand['points'])]
            maxdelta=[max(d[k] for d in deltas) for k in range(4)]
            if any(v>t for v,t in zip(maxdelta,(1e-8,1e-8,1e-8,1e-7))):raise ValueError('Paired state disagreement')
            for variant in ('baseline','contact'):
                v=data[repeat,gate,variant];first=data[0,gate,variant]
                if v['signatures']!=first['signatures'] or v['init']!=first['init']:raise ValueError('Nonrepeatable numerical trajectory')
            by={r['variant']:r for r in rows if r['repeat']==repeat and r['gate_V']==gate}
            ratios={}
            if 'native' in by:
                ratios={v:by[v]['external_wall_seconds']/by['native']['external_wall_seconds'] for v in ('baseline','contact')}
                report['status']='pass' if report['status']=='pass' and all(r<=1.5 for r in ratios.values()) else 'performance_gate_failed'
            pairs.append(dict(repeat=repeat,max_state_delta=maxdelta,wall_ratios_to_native=ratios,
                wall_seconds={v:r['external_wall_seconds'] for v,r in by.items()},
                cpu_seconds={v:r['child_cpu_seconds'] for v,r in by.items()}))
        prefix=None
        if args.first8:
            first=read(args.first8/'summary.json');prefix={}
            for variant in ('baseline','contact'):
                row=next(r for r in first['runs'] if r['gate_V']==gate and r['variant']==variant)
                old=audit_run(candidate_path(row['directory'],args.path_map),args.path_map);new=data[0,gate,variant]
                if old['signatures']!=new['signatures'][:len(old['signatures'])] or old['init']!=new['init']:
                    raise ValueError('Paused first8 differs from fresh full prefix')
                prefix[variant]=dict(exact=True,pause_resume=row.get('pause_resume',False))
        report['gates'][str(gate)]=dict(pairs=pairs,repeated_numerics_exact=len(repeats)>1,repeats_checked=len(repeats),first8_prefix=prefix,
            median_wall_seconds={v:statistics.median(p['wall_seconds'][v] for p in pairs) for v in pairs[0]['wall_seconds']})
    report['hashes'][str(args.matrix/'summary.json')]=sha(args.matrix/'summary.json')
    args.output.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(dict(status=report['status'],gates=report['gates']),indent=2))


if __name__=='__main__':main()
