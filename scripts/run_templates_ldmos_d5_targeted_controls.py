"""Bounded T470p D5 controls; never changes the frozen package or production defaults."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import sys
import time

RUNNER_HASH = '75d7a34556f6448c4e869707c8661c9717550c3058f9580387ebcc303a7274bc'


def targets(start, stop, step):
    if not (0 < step and start < stop): raise ValueError('Invalid interval')
    out = []
    while start < stop-1e-10:
        start = min(start+step, stop)
        out.append(start)
    return out


def secant_rows(current, previous, ratio):
    if not 0 < ratio <= 2+1e-10: raise ValueError('Unprotected extrapolation ratio')
    if len(current) != len(previous): raise ValueError('State lengths differ')
    result = []
    for row, old in zip(current, previous):
        if row['node_id'] != old['node_id']: raise ValueError('State order differs')
        new = dict(row)
        for field in ('psi', 'phin', 'phip'):
            value = float(row[field])+ratio*(float(row[field])-float(old[field]))
            if not math.isfinite(value): raise ValueError('Nonfinite prediction')
            new[field] = format(value, '.17g')
        for field, ref, inc in (('phin','electron_qf_reference_V','electron_qf_increment_V'),
                                ('phip','hole_qf_reference_V','hole_qf_increment_V')):
            new[ref], new[inc] = new[field], '0'
        result.append(new)
    return result


def relocate(value, source, dest):
    if isinstance(value, dict): return {k:relocate(v,source,dest) for k,v in value.items()}
    if isinstance(value, list): return [relocate(v,source,dest) for v in value]
    if isinstance(value, str): return value.replace(str(source),str(dest)).replace(source.as_posix(),dest.as_posix())
    return value


def set_target(cfg, target, step, qf_limit, frame, seed):
    # The outer voltage increment and Newton QF limit are independent factors.
    for c in cfg['contacts']:
        if c['name'] == 'drain': c['bias'] = target-frame
    cfg['sweep'].update(start=target-frame, stop=target-frame, bias_points=[target-frame],
        step=step, initial_step=step, min_step=step, max_step=step, initial_state_file=str(seed))
    cfg['solver']['quasi_fermi_update_limit_V'] = qf_limit


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--mode',choices=['steps','diagnostics'],required=True)
    a=parser.parse_args()
    sys.path.insert(0,str(a.package/'scripts'))
    import run_templates_ldmos_linked_d5 as linked
    from dc_worker_client import DCWorker
    from analyze_templates_ldmos_d5_newton_cost import trace_cost
    runner=a.package/'build-release/vela_example_runner.exe'
    if linked.digest(runner)!=RUNNER_HASH: raise ValueError('Frozen executable mismatch')
    a.output.mkdir(parents=True,exist_ok=False)
    env=os.environ.copy()
    env.update(VELA_LINEAR_SOLVER='umfpack',VELA_LINEAR_THREADS='1',VELA_BLAS_THREADS='1',
        OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',OMP_DYNAMIC='FALSE',OMP_MAX_ACTIVE_LEVELS='1',
        VELA_LINEAR_FACTOR_STATS='0')
    summary=dict(status='running',mode=a.mode,runner_sha256=RUNNER_HASH,
                 script_sha256=linked.digest(Path(__file__)),cases=[])
    def save(): linked.write(a.output/'summary.json',summary)
    def execute(worker,src,dest,cfg,target):
        linked.write(dest/'control.json',cfg)
        start=time.perf_counter()
        response,cpu=worker.run(dest/'control.json',timeout=240)
        elapsed=time.perf_counter()-start
        status=dict(returncode=response['returncode'],curve=linked.rows(dest/'curve.csv'))
        trace=linked.rows(dest/'newton_iterations.csv')
        counts,events=trace_cost(trace,cfg['solver']['block_absolute_convergence'],
                                 cfg['solver']['carrier_row_convergence']['eps_row'])
        attempts=linked.rows(dest/'newton_attempts.csv')
        good=linked.good(status,dest,target)
        record=dict(directory=str(dest),target_V=target,passed=good,wall_seconds=elapsed,cpu=cpu,
            failed_attempts=sum(x['status']!='accepted' for x in attempts),events=events,
            config_sha256=linked.digest(dest/'control.json'),**counts)
        print(a.mode,dest.name,target,counts['updates'],good,flush=True)
        return record
    for gate in (4,8):
        curve=a.package/f'results/t470p/d5_p31_r0_vg{gate}_U1'
        ledger=linked.read(curve/'fixed/ledger.json')
        if a.mode=='steps':
            for start in (20.,32.):
                pick=lambda v:min(ledger['exact_points'],key=lambda p:abs(p['bias_V']-v))
                previous,current,final=map(pick,(start-4/3,start,start+4/3))
                for maximum,qf in ((.2,.2),(.4,.2),(.8,.2),(4/3,.2),(.8,.4)):
                    name=f'vg{gate}_from{start:g}_step{maximum:.6f}_qf{qf:g}'
                    dest=a.output/name; dest.mkdir()
                    case=dict(name=name,gate=gate,start_V=current['bias_V'],stop_V=final['bias_V'],
                        max_step_V=maximum,qf_limit_V=qf,initial_history='previous_exact_output',runs=[])
                    summary['cases'].append(case);save()
                    src=curve/final['case'];base=linked.read(src/'control.json')
                    frame=linked.frame_of(src)
                    if abs(linked.frame_of(Path(current['state']).parent)-frame)>1e-10:
                        raise ValueError('Window crosses a reference-frame boundary')
                    old_path,parent=Path(previous['state']),Path(current['state'])
                    for point in (previous,current,final):
                        if linked.digest(Path(point['state']))!=point['sha256']: raise ValueError('Frozen state changed')
                    old_v,parent_v=previous['bias_V'],current['bias_V']
                    worker=DCWorker(runner,dest,env,linked.cpu_times,reuse_linear_analysis=True)
                    try:
                        for index,target in enumerate(targets(parent_v,final['bias_V'],maximum)):
                            child=dest/f'point_{index:03d}'; child.mkdir()
                            ratio=(target-parent_v)/(parent_v-old_v)
                            seed=child/'prediction.csv'
                            linked.csv_out(seed,secant_rows(linked.rows(parent),linked.rows(old_path),ratio))
                            cfg=relocate(base,src,child)
                            set_target(cfg,target,target-parent_v,qf,frame,seed)
                            record=execute(worker,src,child,cfg,target)
                            record.update(ratio=ratio,parent_sha256=linked.digest(parent),
                                          history_sha256=linked.digest(old_path),seed_sha256=linked.digest(seed))
                            case['runs'].append(record);save()
                            if not record['passed']: break
                            old_path,parent=parent,child/'state.csv';old_v,parent_v=parent_v,target
                    finally: worker.close()
                    case['completed']=abs(parent_v-final['bias_V'])<1e-9
                    if case['completed']:
                        delta=linked.state_difference(parent,Path(final['state']))
                        case['state_difference']=delta
                        case['equivalent']=max(delta[k]['max_absolute'] for k in ('psi','phin','phip'))<=1e-8
                        now,kcl=linked.ports(parent.parent,final['bias_V']-frame)
                        ref,_=linked.ports(Path(final['state']).parent,final['bias_V']-frame)
                        case['current_relative_difference']=max(abs(now[k]-ref[k]) for k in ref)/max(map(abs,ref.values()))
                    save()
        else:
            predicted=[r for r in ledger['runs'] if r.get('outer_predictor')]
            for voltage in (.8368462059131809,1.0368462059131809,20.235,33.5683333333333):
                r=min(predicted,key=lambda p:abs(p['target_V']-voltage))
                src=curve/r['case']; base=linked.read(src/'control.json')
                violations=linked.rows(src/'carrier_violations.csv')
                nodes=list(dict.fromkeys(int(x['node_id']) for x in sorted(violations,
                           key=lambda x:float(x['ratio']),reverse=True)))[:16]
                # Full node traces permit common fixed weights and actual worst-row localization.
                all_nodes=[int(x['node_id']) for x in linked.rows(src/'state.csv')]
                for variant in ('guarded','none') if voltage<2 else ('guarded',):
                    dest=a.output/f'vg{gate}_{src.name}_{variant}';dest.mkdir()
                    cfg=relocate(base,src,dest)
                    seed=Path(r['parent_state'] if variant=='guarded' else r['outer_predictor']['current_state'])
                    cfg['sweep']['initial_state_file']=str(seed)
                    cfg['solver']['carrier_row_convergence'].update(trace_csv=str(dest/'row_trace.csv'),
                        trace_nodes=all_nodes,trace_first_iterations=160,trace_every_iterations=1)
                    cfg['solver']['local_update_diagnostics']=dict(enabled=True,csv_file=str(dest/'local_updates.csv'),
                        nodes=nodes,first_iterations=160,every_iterations=1)
                    worker=DCWorker(runner,dest,env,linked.cpu_times,reuse_linear_analysis=True)
                    try: record=execute(worker,src,dest,cfg,r['target_V'])
                    finally: worker.close()
                    record.update(gate=gate,variant=variant,source=str(src),nodes=nodes,
                        initial_seed_sha256=linked.digest(seed),trace_sha256=linked.digest(dest/'row_trace.csv'))
                    if record['passed']:
                        delta=linked.state_difference(dest/'state.csv',src/'state.csv')
                        record['state_difference']=delta
                        record['equivalent']=max(delta[k]['max_absolute'] for k in ('psi','phin','phip'))<=1e-8
                    summary['cases'].append(record);save()
    summary['status']='completed'
    save()


if __name__=='__main__':main()
