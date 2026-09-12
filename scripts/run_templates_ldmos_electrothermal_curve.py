"""Continue an experimental four-equation LDMOS curve through exact bias points.
Accepts only original electrical block/row/KCL gates and approved heat balance.
Native electrical curve and temperature-field comparison remains a separate gate.
"""
import argparse,csv,hashlib,json,math,os,subprocess,time
from pathlib import Path


def load(path):return json.loads(path.read_text(encoding='utf-8'))
def save(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp');temporary.write_text(json.dumps(value,indent=2),encoding='utf-8');temporary.replace(path)
def hash_file(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def state_gate(result,bias):
    reasons=[]
    if result.get('diagnostic_stop')!='diagnostic_scaled_residual':reasons.append('four_equation_convergence')
    blocks=result.get('electrical_block_gates',[])
    if len(blocks)!=3 or not all(b['satisfied'] for b in blocks):reasons.append('electrical_blocks')
    if not result.get('carrier_row_gate',{}).get('satisfied',False):reasons.append('carrier_rows')
    contacts=result.get('contacts',[])
    if {c['contact'] for c in contacts}!={'gate','drain','source','substrate'} or len(contacts)!=4:
        return dict(pass_gate=False,reasons=reasons+['terminal_set'])
    currents=[c['total_outflow_A_per_m'] for c in contacts]
    kcl=abs(sum(currents))/max(max(map(abs,currents)),1e-24)
    if not math.isfinite(kcl) or kcl>1e-3:reasons.append('kcl')
    source=result['lattice_source_W_per_m'];boundary=result['boundary_heat_W_per_m']
    temperature=result['temperature_K']
    if not temperature or not all(math.isfinite(t) and t>=50. for t in temperature):reasons.append('temperature')
    if bias==0.:
        balance=None
        if source!=0. or boundary!=0. or any(t!=300. for t in temperature):reasons.append('zero_power_equilibrium')
    else:
        scale=max(abs(source),abs(boundary));balance=abs(source-boundary)/scale if scale>0. else None
        if balance is None or not math.isfinite(balance) or balance>1e-3:reasons.append('heat_balance')
    return dict(pass_gate=not reasons,reasons=reasons,kcl_relative=kcl,heat_balance_relative=balance)


def write_curves(root,points):
    with (root/'curve.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['bias_V','current_total_A_per_um','peak_temperature_K','mean_temperature_K']);w.writeheader()
        for point in points:
            r=load(Path(point['result']));d=next(c for c in r['contacts'] if c['contact']=='drain');areas=r['nodal_area_m2']
            w.writerow(dict(bias_V=point['bias_V'],current_total_A_per_um=d['total_outflow_A_per_m']*1e-6,peak_temperature_K=max(r['temperature_K']),mean_temperature_K=sum(a*t for a,t in zip(areas,r['temperature_K']))/sum(areas)))
    with (root/'terminal_balance.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['point_index','bias_V','contact','current_total_A_per_um']);w.writeheader()
        for index,point in enumerate(points):
            for contact in load(Path(point['result']))['contacts']:
                w.writerow(dict(point_index=index,bias_V=point['bias_V'],contact=contact['contact'],current_total_A_per_um=contact['total_outflow_A_per_m']*1e-6))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('probe','input','bias-points','output'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--initial-step',type=float,default=.1)
    parser.add_argument('--maximum-step',type=float,default=4./3.)
    parser.add_argument('--minimum-step',type=float,default=1e-4)
    parser.add_argument('--max-newton',type=int,default=60)
    args=parser.parse_args();base=load(args.input);biases=load(args.bias_points)
    if not biases or biases[0]!=0. or any(not math.isfinite(v) for v in biases) or any(a>=b for a,b in zip(biases,biases[1:])):raise ValueError('Exact biases must increase strictly from zero')
    if not 0.<args.minimum_step<=args.initial_step<=args.maximum_step:raise ValueError('Invalid step limits')
    if args.max_newton<1:raise ValueError('Positive Newton budget required')
    args.output.mkdir(parents=True,exist_ok=False)
    mesh=load(Path(base['mesh_file']));drain=set(next(c for c in mesh['contacts'] if c['name']=='drain')['node_ids']);origin=base.get('potential_origin_V',0.)
    for b in base['boundaries']:
        if b['node'] in drain and b['kind']=='neutral_contact' and b['value']+origin!=0.:raise ValueError('Initial drain bias must be zero')
    env=dict(os.environ)
    if os.name=='nt':env['PATH']='D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+env['PATH']
    ledger=dict(scope='Experimental D0 continuation; full native curve/field gates scored separately',status='running',input_sha256=hash_file(args.input),probe_sha256=hash_file(args.probe),biases=biases,runs=[],exact_points=[],accepted_bias_V=0.,accepted_result=None)
    save(args.output/'ledger.json',ledger)
    reference_keys=('referenced_state_interleaved','electron_qf_reference_V','hole_qf_reference_V')
    reference_pack={k:base[k] for k in reference_keys if k in base}
    state=base['state_interleaved'];current=0.;step=args.initial_step;attempt=0;index=0
    while index<len(biases):
        target=0. if index==0 else min(biases[index],current+step)
        directory=args.output/f'step_{attempt:04d}';directory.mkdir();cfg=dict(base)
        cfg['boundaries']=[dict(b) for b in base['boundaries']]
        for b in cfg['boundaries']:
            if b['node'] in drain and b['kind']=='neutral_contact':b['value']=target-origin
        cfg['state_interleaved']=state;cfg['diagnostic_newton_max_iterations']=args.max_newton
        for key in reference_keys:cfg.pop(key,None)
        cfg.update(reference_pack)
        save(directory/'input.json',cfg);start=time.perf_counter()
        with (directory/'run.log').open('w',encoding='utf-8') as log:
            process=subprocess.Popen([str(args.probe.resolve()),str((directory/'input.json').resolve()),str((directory/'output.json').resolve())],stdout=log,stderr=subprocess.STDOUT,env=env)
            ledger['active_child']=dict(pid=process.pid,directory=str(directory.resolve()),bias_V=target);save(args.output/'ledger.json',ledger)
            returncode=process.wait()
        record=dict(parent_bias_V=current,bias_V=target,directory=str(directory.resolve()),returncode=returncode,wall_seconds=time.perf_counter()-start)
        result_path=directory/'output.json'
        result=load(result_path) if returncode==0 and result_path.exists() else None
        gate=state_gate(result,target) if result else dict(pass_gate=False,reasons=['probe_failed'])
        record.update(gate=gate,newton_updates=result['newton_updates'] if result else None)
        ledger['runs'].append(record);ledger['active_child']=None;attempt+=1
        if gate['pass_gate']:
            current=target;state=result['state_interleaved'];reference_pack={k:result[k] for k in reference_keys if k in result};ledger['accepted_bias_V']=current;ledger['accepted_result']=str(result_path.resolve())
            if abs(current-biases[index])<1e-12:
                ledger['exact_points'].append(dict(bias_V=biases[index],result=str(result_path.resolve())));index+=1;write_curves(args.output,ledger['exact_points'])
            if result['newton_updates']<=8:step=min(args.maximum_step,step*1.5)
            elif result['newton_updates']>20:step=max(args.minimum_step,step*.5)
        else:
            step=(target-current)*.5
            if target==0. or step<args.minimum_step:
                ledger['status']='failed';save(args.output/'ledger.json',ledger);raise RuntimeError(f'Gate failure at {target}: {gate}')
        save(args.output/'ledger.json',ledger)
        print(json.dumps(dict(accepted_bias_V=current,exact_points=len(ledger['exact_points']),attempted_bias_V=target,gate=gate,newton_updates=record['newton_updates'])),flush=True)
    ledger['status']='complete';save(args.output/'ledger.json',ledger)

if __name__=='__main__':main()
