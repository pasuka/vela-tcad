"""One frozen restart of each rejected finite state, with unchanged criteria."""
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import validate_simplemos_srh_finite_20260908 as t
a=t.a;d=t.d;v=t.v;LOCAL=t.LOCAL/'restart';OUT=t.OUT/'restart'

def prepare():
    a.verify(t.OUT/'validation_evidence.json');failed=[r for r in a.rows(t.OUT/'finite_dc.csv') if r['qualified']!='True'];jobs=[]
    files=[Path(__file__).resolve(),t.OUT/'validation_evidence.json',t.RUNNER]
    for r in failed:
        old=t.LOCAL/r['key']/r['axis']/r['label'];dest=LOCAL/r['key']/r['axis']/r['label'];cfg=a.read(old/'config.json');cfg.update(state_file=str(old/'state.csv'),output_state_file=str(dest/'state.csv'))
        status=a.read(old/'config.status.json')
        nodes=sorted({int(x['node_id']) for x in status['carrier_row_convergence'].get('violations',[])})
        cfg['solver']['local_update_diagnostics']=dict(enabled=True,csv_file=str(dest/'updates.csv'),nodes=nodes,first_iterations=200,every_iterations=1)
        a.write(dest/'config.json',cfg);v.post_config(cfg,dest)
        files += [old/'state.csv',dest/'config.json',dest/'all_row.json',dest/'acceptance_edges.json']+v.probes(cfg,dest/'post',dest/'state.csv')
        jobs.append(dict(key=r['key'],axis=r['axis'],label=r['label'],dest=str(dest),previous_max_ratio=float(r['max_row_ratio'])))
    a.write(OUT/'contract.json',dict(jobs=jobs,maximum_restarts_per_failed_path=1,policy='Same finite parameter and numerical solver configuration, exported failed state as new initial condition; read-only raw/capped/applied update diagnostics enabled at violated rows. No failed artifact overwritten, no tolerance/line-search change.'))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json']);print('Frozen one restart per failed path:',len(jobs),flush=True)

def run():
    a.verify(OUT/'freeze.json');cases={c['key']:c for c in a.read(t.OUT/'contract.json')['cases']}
    def one(job):
        case=cases[job['key']];dest=Path(job['dest']);status=t.execute(dest/'config.json',case,job['axis'],1.)
        row=dict(key=job['key'],axis=job['axis'],label=job['label'],qualified=False,failure=status.get('failure_reason',''))
        if (dest/'state.csv').exists():
            for name in ('all_row','acceptance_edges'):t.execute(dest/(name+'.json'),case,job['axis'],1.)
            row.update(t.c.s.prior.old.qualify(case,dest))
            for name in ('functional','edges','terms'):t.execute(dest/'post'/(name+'.json'),case,job['axis'],1.)
        print(job['key'],job['axis'],job['label'],row['qualified'],row.get('max_row_ratio'),flush=True);return row
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(one,a.read(OUT/'contract.json')['jobs']))
    keys=list(dict.fromkeys(k for r in rows for k in r));a.write_csv(OUT/'dc.csv',[{k:r.get(k,'') for k in keys} for r in rows])

def analyze():
    initial=a.rows(t.OUT/'finite_dc.csv');rerun=a.rows(OUT/'dc.csv');selected=[];pairs=[];fields=[];local=[];failures=[]
    for case in a.read(t.OUT/'contract.json')['cases']:
        geo,mask=v.m.previous.prior.support(case);nodes=[case['mapped_nodes'][tag]['node'] for tag in ('792','1057')]+[1000,1009]
        for axis in t.AXES:
            states=[];records=[]
            tags=['792','1057'] if axis=='joint' else [axis[5:]];ratios={case['mapped_nodes'][tag]['node']:case['mapped_nodes'][tag]['ratio'] for tag in tags}
            for label in ('finite','finite_native_seed'):
                base=next(r for r in initial if r['key']==case['key'] and r['axis']==axis and r['label']==label)
                recovered=next((r for r in rerun if r['key']==case['key'] and r['axis']==axis and r['label']==label),None)
                row=recovered or base;dest=(LOCAL if recovered else t.LOCAL)/case['key']/axis/label
                r=dict(row,path=a.rel(dest),restarted=recovered is not None,initial_qualified=base['qualified']=='True');selected.append(r);records.append(r)
                states.append(d.ordered(dest/'state.csv',geo.count));met,loc=t.c.field_metrics(case,axis+'/'+label,dest,ratios,nodes)
                for rr in met+loc:rr['qualified_state']=row['qualified']=='True'
                fields+=met;local+=loc
                if row['qualified']!='True':
                    status=a.read(dest/'config.status.json');rv=status['carrier_row_convergence'];violations=rv.get('violations',[])
                    for violation in violations:
                        failures.append(dict(key=case['key'],axis=axis,label=label,iteration=status['iterations'],failure=status['failure_reason'],psi_block=status['final_block_residuals']['psi'],**violation))
            delta=v.m.previous.prior.delta_states(*states,mask);ie=abs(float(records[0]['current_A_per_um'])/float(records[1]['current_A_per_um'])-1)
            pairs.append(dict(key=case['key'],axis=axis,**delta,Id_relative=ie,qualified=all(r['qualified']=='True' for r in records) and max(delta[k] for k in ('psi_max_V','phin_max_V','phip_max_V'))<=1e-6 and delta['density_max_relative']<=1e-4 and ie<=1e-6))
    for name,rows in [('selected_states',selected),('dual',pairs),('fields',fields),('local_fields',local)]:a.write_csv(OUT/(name+'.csv'),rows)
    if failures:
        keys=list(dict.fromkeys(k for r in failures for k in r));a.write_csv(OUT/'remaining_violations.csv',[{k:r.get(k,'') for k in keys} for r in failures])
    summary=dict(restarts=len(rerun),qualified_restarts=sum(r['qualified']=='True' for r in rerun),qualified_endpoints=sum(r['qualified']=='True' for r in selected),endpoints=len(selected),qualified_dual=sum(r['qualified'] for r in pairs),initial_failures_retained=sum(r['qualified']!='True' for r in initial),qualified=all(r['qualified']=='True' for r in selected) and all(r['qualified'] for r in pairs))
    a.write(OUT/'summary.json',summary);d.matrix.freeze(OUT/'validation_evidence.json',[OUT/'freeze.json',OUT/'summary.json']+list(OUT.glob('*.csv'))+[x for x in LOCAL.rglob('*') if x.is_file()]);print(summary,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run','analyze'));globals()[parser.parse_args().action]()
