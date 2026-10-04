"""Serial qualification of electrothermal sparse-structure reuse.

Requires the explicitly supplied prior HDF5 qualification directory, including
its audited D0 inputs, native reference fields and frozen-300 K trajectories.
Never modifies that evidence. An unsuccessful gate stops subsequent stages.
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

from electrothermal_state import read_bound_record, write as write_state
from state_archive import mesh_identity
from run_templates_ldmos_electrothermal_curve import state_gate
from analyze_templates_ldmos_d0 import analyze
if os.name=='nt':
    from run_templates_ldmos_cost_probe import windows_peak_rss
    from run_templates_ldmos_linked_d5 import cpu_times

FIELDS=('state_interleaved','referenced_state_interleaved','electron_qf_reference_V','hole_qf_reference_V','temperature_K')

def read(path): return json.loads(Path(path).read_text(encoding='utf-8'))
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,value):
    path=Path(path);temp=path.with_suffix('.tmp');temp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8');temp.replace(path)

def compare(a,b):
    errors={key:max(abs(x-y) for x,y in zip(a[key],b[key])) for key in FIELDS}
    if any(len(a[k])!=len(b[k]) for k in FIELDS) or not all(v<=1e-8 for v in errors.values()):
        raise ValueError('Saved physical/reference state disagreement: '+str(errors))
    return dict(maximum_absolute_errors=errors,bitwise_values_equal=all(a[k]==b[k] for k in FIELDS),
                updates=[a['newton_updates'],b['newton_updates']])

def electrical_gate(result,bias,cold=False):
    if not cold:
        gate=state_gate(result,bias)
        if not gate['pass_gate']:raise ValueError(str(gate))
    else:
        if result['diagnostic_stop']!='diagnostic_scaled_residual' or not result['carrier_row_gate']['satisfied'] or not all(v['satisfied'] for v in result['electrical_block_gates']):
            raise ValueError('Frozen-temperature electrical gate failed')
        currents=[c['total_outflow_A_per_m'] for c in result['contacts']]
        if abs(sum(currents))/max(max(map(abs,currents)),1e-24)>1e-3:raise ValueError('Frozen-temperature KCL failed')
        if not all(t==300. for t in result['temperature_K']):raise ValueError('Temperature not frozen')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('runner','prior','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--stage',choices=('points','short','full','repeats'),required=True)
    args=p.parse_args();runner=args.runner.resolve();probe=runner.parent/('electrothermal_probe.exe' if os.name=='nt' else 'electrothermal_probe')
    prior=args.prior.resolve();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).resolve().parents[1]
    status=dict(status='running',pid=os.getpid(),stage=args.stage,started=time.time(),runs=[],comparisons=[],joint=[])
    frozen={str(v):sha(v) for v in (runner,probe,Path(__file__),prior/'final_curves/d0/vg4.json',prior/'final_curves/d0/vg8.json')}
    status['frozen_sha256']=frozen
    os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_DYNAMIC='FALSE',VELA_LINEAR_THREADS='1',VELA_BLAS_THREADS='1')
    def save():write(out/'summary.json',status)
    def run(name,command,expected=0):
        status.update(current=name,command=list(map(str,command)));save();start=time.perf_counter()
        with (out/(name+'.log')).open('x') as log:
            child=subprocess.Popen(list(map(str,command)),stdout=log,stderr=subprocess.STDOUT,cwd=source)
            status['child_pid']=child.pid;save();code=child.wait()
        row=dict(name=name,wall_seconds=time.perf_counter()-start,returncode=code)
        if os.name=='nt':row.update(peak_working_set_bytes=windows_peak_rss(child),cpu_seconds=cpu_times(child))
        status['runs'].append(row);status.pop('child_pid',None);save()
        if code!=expected:raise RuntimeError(f'{name}: exit {code}, expected {expected}')
        return row
    def point_cfg(path,enabled):
        cfg=read_bound_record(path);cfg['reuse_jacobian_structure']=enabled;cfg['performance_profiling']=True
        return cfg
    def state_input(path,cfg):
        metadata=dict(mode='electrothermal',potential_origin_V=cfg.get('potential_origin_V',0.),
            mesh_sha256=mesh_identity(read(cfg['mesh_file']),cfg.get('coordinate_to_metres',1.)))
        write_state(path,cfg,metadata)
    save()
    try:
        if args.stage=='points':
            for cold in (True,False):
                for gate in (4,8):
                    directory=prior/(f'cold300_electrical/vg{gate}' if cold else f'final_curves/d0/results_vg{gate}')
                    ledger=read(directory/'ledger.json')
                    for index in (1,15,30):
                        point=ledger['exact_points'][index];results=[]
                        for enabled in (False,True):
                            name=f'{"cold" if cold else "hot"}_vg{gate}_p{index}_{int(enabled)}'
                            folder=out/name;folder.mkdir();cfg=point_cfg(Path(point['result']).parent/'input.json',enabled)
                            state_input(folder/'input.json',cfg);run(name,[probe,folder/'input.json',folder/'output.json'])
                            result=read_bound_record(folder/'output.json');electrical_gate(result,point['bias_V'],cold);results.append(result)
                            for key in FIELDS[:-1]:cfg[key]=result[key]
                            for key in ('diagnostic_predictor_candidates','diagnostic_tangent_predictor','temperature_K'):cfg.pop(key,None)
                            state_input(folder/'restart.json',cfg);run(name+'_restart',[probe,folder/'restart.json',folder/'restarted.json'])
                            repeated=read_bound_record(folder/'restarted.json');electrical_gate(repeated,point['bias_V'],cold)
                            status['comparisons'].append(dict(case=name,scope='restart',**compare(result,repeated)))
                        status['comparisons'].append(dict(gate=gate,index=index,cold=cold,scope='toggle',**compare(*results)));save()
        else:
            rounds=3 if args.stage=='repeats' else 1;points=8 if args.stage=='short' else 31
            previous={}
            for repetition in range(rounds):
                paired={}
                for gate in (4,8):
                    for enabled in ((False,True) if (repetition+gate//4)%2 else (True,False)):
                        name=f'r{repetition}_vg{gate}_{int(enabled)}';folder=out/name;folder.mkdir()
                        deck=read(prior/f'final_curves/d0/vg{gate}.json')
                        cfg=point_cfg(prior/f'final_curves/d0/input_vg{gate}.json',enabled)
                        state_input(folder/'input.json',cfg)
                        deck.update(input_file=str(folder/'input.json'),output_directory=str(folder/'results'),resume=False,pause_after_attempts=0,
                                    reuse_jacobian_structure=enabled)
                        deck['sweep']['bias_points_V']=deck['sweep']['bias_points_V'][:points]
                        write(folder/'deck.json',deck)
                        row=run(name,[runner,'--config',folder/'deck.json'])
                        ledger=read(folder/'results/ledger.json')
                        if ledger['status']!='complete' or len(ledger['exact_points'])!=points:raise ValueError('Incomplete curve')
                        paired[gate,enabled]=folder/'results'
                        counters={};total_updates=0;failures=0
                        for attempt in ledger['runs']:
                            result=read_bound_record(Path(attempt['directory'])/'output.json');total_updates+=result['newton_updates']
                            failures+=not attempt['gate']['pass_gate']
                            if attempt['gate']['pass_gate']:electrical_gate(result,attempt['bias_V'])
                            for key in ('assembly_seconds','full_assembly_seconds','residual_assembly_seconds','assembly_calls','factorization_seconds','factorizations','symbolic_analyses','linear_solve_seconds','jacobian_structure_builds','jacobian_structure_hits','heat_structure_builds','heat_structure_hits'):
                                counters[key]=counters.get(key,0)+result.get('performance',{}).get(key,0)
                        row.update(gate=gate,enabled=enabled,repetition=repetition,points=points,updates=total_updates,failed_attempts=failures,counters=counters)
                        for i,entry in enumerate(ledger['exact_points']):
                            result=read_bound_record(entry['result']);electrical_gate(result,entry['bias_V'])
                            if (gate,enabled,i) in previous:
                                status['comparisons'].append(dict(scope='repeat',gate=gate,index=i,enabled=enabled,**compare(previous[gate,enabled,i],result)))
                            else:previous[gate,enabled,i]=result
                        save()
                    left=read(paired[gate,False]/'ledger.json');right=read(paired[gate,True]/'ledger.json')
                    for a,b in zip(left['exact_points'],right['exact_points']):
                        if a['bias_V']!=b['bias_V']:raise ValueError('Exact bias mismatch')
                        status['comparisons'].append(dict(scope='toggle',gate=gate,bias_V=a['bias_V'],repetition=repetition,
                            **compare(read_bound_record(a['result']),read_bound_record(b['result']))))
                    row['trajectory_equal']=[(v['bias_V'],v['gate']['pass_gate']) for v in left['runs']]==[(v['bias_V'],v['gate']['pass_gate']) for v in right['runs']];save()
                if points==31:
                    for enabled in (False,True):
                        report=analyze(prior/'d0_native',{g:paired[g,enabled] for g in (4,8)},
                            read(source/'reference_tcad/templates_ldmos_sentaurus2022/thermal/d0_thermal_acceptance.json'),15)
                        reportfile=out/f'joint_r{repetition}_{int(enabled)}.json';write(reportfile,report)
                        status['joint'].append(dict(repetition=repetition,enabled=enabled,status=report['status'],file=str(reportfile)))
                        if report['status']!='pass':raise ValueError('Native joint gates failed')
            totals={str(enabled):sum(v['wall_seconds'] for v in status['runs'] if v.get('enabled')==enabled) for enabled in (False,True)}
            status['paired_wall_seconds']=totals;status['wall_change_percent']=100*(totals['True']/totals['False']-1)
        if not all(sha(path)==value for path,value in frozen.items()):raise ValueError('Frozen executable/input changed')
        status.update(status='completed',finished=time.time());save()
    except BaseException as e:
        status.update(status='failed',error=repr(e),traceback=traceback.format_exc(),finished=time.time());save();raise

if __name__=='__main__':main()
