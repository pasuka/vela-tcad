"""Qualify the explicit PhuMob box candidate without fitting the native G floor."""
import argparse
import copy
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal as D, localcontext
from pathlib import Path
import numpy as np
import validate_simplemos_phumob_chain_fix_20260909 as previous
import audit_simplemos_phumob_cross_derivatives_20260908 as hp

a,d,q=previous.a,previous.d,previous.q
REPO=previous.p.REPO
LOCAL=REPO/'build-release/phumob_box_20260909'
OUT=REPO/'reference_tcad/simplemos_sentaurus2022/phumob_box_20260909'
NODES=hp.NODES
VT=previous.cross.VT


def prepare():
    a.verify(OUT/'baseline.json')
    jobs=[];fixed=[];files=[OUT/'baseline.json',OUT/'scope.json',Path(__file__).resolve(),Path(hp.__file__)]
    for old in a.read(q.OUT/'validation_contract.json')['jobs']:
        if old['model']!='old_slotboom' or old['index'] not in (40,50):continue
        job=dict(old);job['model']='phumob';key=f"{old['case']}_vg_{old['index']:03d}"
        cfg=a.read(Path(old['original_config']))
        cfg['solver']['mobility'].update(model='phumob',edge_averaging='element_box_phumob')
        path=LOCAL/'inputs'/key/old['arm']/'config.json'
        a.write(path,cfg);job['original_config']=str(path);jobs.append(job)
        files += [path,Path(job['seed'])]+[Path(cfg[k]) for k in ('mesh_file','materials_file','node_doping_file')]
        if old['arm']!='vela':continue
        mesh=a.read(Path(cfg['mesh_file']))
        contacts={i for c in mesh['contacts'] for i in c['node_ids']}
        assert not contacts.intersection(NODES)
        dest=LOCAL/'fixed'/key
        probe=copy.deepcopy(cfg);probe['state_file']=job['seed']
        probe.pop('output_state_file',None)
        probe['solver']['recombination']=['none'];probe['solver']['srh_doping_dependence']['enabled']=False
        for name,simulation in [('edges','sg_edge_flux_probe'),('cells','element_box_probe'),('mobility','edge_mobility_probe')]:
            pc=copy.deepcopy(probe);pc.update(simulation_type=simulation,output_csv=str(dest/(name+'.csv')))
            a.write(dest/(name+'.json'),pc);files.append(dest/(name+'.json'))
        neighbors=set(NODES)
        for cell in mesh['triangles']:
            nodes=cell['nodes'] if 'nodes' in cell else cell['node_ids']
            if set(nodes).intersection(NODES):neighbors.update(nodes)
        _,mask=q.run.V.m.previous.prior.support(job)
        pc=copy.deepcopy(probe);pc.update(simulation_type='newton_jvp_probe',output_csv=str(dest/'jvp.csv'),row_output_csv=str(dest/'rows.csv'))
        pc['directions']=[dict(name=f'{mode}_{node}_{step:.0e}',mode=mode,node_ids=[node],amplitude_V=step,exclude_contacts=True)
            for node in NODES for mode in ('phin','phip') for step in (1e-5,3e-6)]
        pc['directions'] += [dict(name=f'block_{mode}_{step:.0e}',mode=mode,node_ids=np.where(mask)[0][::3].tolist(),amplitude_V=step,exclude_contacts=True)
            for mode in ('psi','phin','phip') for step in (1e-4,3e-5,1e-5,3e-6)]
        pc['sample_rows']=[dict(block=block,node_id=node) for block in ('phin','phip') for node in sorted(neighbors-contacts)]
        a.write(dest/'jvp.json',pc);files.append(dest/'jvp.json')
        fixed.append(dict(**job,key=key,dir=str(dest),contacts=sorted(contacts)))
    assert len(jobs)==16 and len(fixed)==8
    a.write(OUT/'contract.json',dict(jobs=jobs,fixed=fixed,gates=a.read(OUT/'scope.json')['gates'],
        derivative_gates=dict(cell_and_edge_mu=1e-12,cross_relative=1e-5,base_closure=1e-10,block_relative=1e-4),
        native_qualification='Independent exact Vela formula qualification permits diagnostic self-consistency. Native G-floor mismatch remains separately failed; not full native model acceptance.'))
    d.matrix.freeze(OUT/'input_freeze.json',files+[OUT/'contract.json'])
    print('Prepared eight points, sixteen initialization jobs',flush=True)


def fixed():
    a.verify(OUT/'input_freeze.json')
    sources=[REPO/n for n in ('src/physics/MobilityModel.cpp','src/equation/CoupledDDAssembler.cpp',
        'include/vela/equation/AssemblerUtils.h','src/post/ContactCurrent.cpp','src/tools/vela_example_runner.cpp')]
    d.matrix.freeze(OUT/'fixed_freeze.json',[OUT/'input_freeze.json',q.run.RUNNER]+sources)
    def one(job):
        rows=[]
        for name in ('cells','edges','mobility','jvp'):
            path=Path(job['dir'])/(name+'.json');status=q.run.V.execute(path,q.run.RUNNER,q.run.V.environment())
            rows.append(dict(key=job['key'],probe=name,**status))
            if status['exit_code']!=0:break
        return rows
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(one,a.read(OUT/'contract.json')['fixed']) for r in group]
    a.write(OUT/'fixed_execution.json',dict(results=rows))
    d.matrix.freeze(OUT/'fixed_raw_evidence.json',[OUT/'fixed_freeze.json',OUT/'fixed_execution.json']+[p for p in (LOCAL/'fixed').rglob('*') if p.is_file()])
    assert len(rows)==32 and all(r['exit_code']==0 for r in rows),rows
    print('All 32 fixed-state probes completed',flush=True)


def analyze_fixed():
    a.verify(OUT/'fixed_raw_evidence.json')
    cellchecks=[];edgechecks=[];cross=[];closure=[];blocks=[];third=[]
    with localcontext() as ctx:
        ctx.prec=100
        for job in a.read(OUT/'contract.json')['fixed']:
            path=Path(job['dir']);cfg=a.read(path/'edges.json')
            doping={int(r['node_id']):r for r in a.rows(Path(cfg['node_doping_file']))}
            edges=a.rows(path/'edges.csv');populations={};cells=defaultdict(list)
            for e in edges:
                for k in range(2):populations[int(e[f'node{k}'])]=[D(e[f'{car}_density{k}_m3'])*D('1e-6') for car in ('electron','hole')]
            for row in a.rows(path/'cells.csv'):cells[int(row['cell_id'])].append(row)
            active={cid:rs for cid,rs in cells.items() if float(rs[0]['electron_mobility_m2_V_s'])>0 or float(rs[0]['hole_mobility_m2_V_s'])>0}
            nodal={}
            for node in {int(r['node_id']) for rs in active.values() for r in rs}:
                state=[D(doping[node]['donors_cm3']),D(doping[node]['acceptors_cm3']),*populations[node]]
                nodal[node]=[hp.hp(state,car) for car in (0,1)]
            cellmu={};parts=defaultdict(list)
            for cid,rs in active.items():
                rs=sorted(rs,key=lambda r:int(r['local_vertex']));volume=sum(D(r['measure_m2']) for r in rs)
                assert volume>0
                weights={int(r['node_id']):D(r['measure_m2'])/volume for r in rs}
                cellmu[cid]=[sum(w*nodal[n][car][0] for n,w in weights.items()) for car in (0,1)]
                for car,name in enumerate(('electron','hole')):
                    error=float(abs(cellmu[cid][car]/(D(rs[0][name+'_mobility_m2_V_s'])*D('1e4'))-1))
                    cellchecks.append(dict(key=job['key'],cell=cid,carrier=name,relative=error))
                for k,r in enumerate(rs):
                    edge=tuple(sorted((int(r['node_id']),int(rs[(k+1)%3]['node_id']))))
                    parts[edge].append((cid,D(r['coefficient_next']),weights))
            ref=defaultdict(lambda:D(0));fluxes=defaultdict(list)
            mobility={int(r['edge_id']):r for r in a.rows(path/'mobility.csv')}
            for e in edges:
                i,j=int(e['node0']),int(e['node1']);part=parts[tuple(sorted((i,j)))];geometry=sum(g for _,g,_ in part)
                for car,(name,block,mode,pop,sign) in enumerate([('electron','phin','phip',1,1),('hole','phip','phin',0,-1)]):
                    flux=D(e[name+'_flux']);fluxes[block,i].append(flux);fluxes[block,j].append(-flux)
                    mu=D(e[name+'_mobility_m2_V_s'])*D('1e4')
                    if not mu:
                        assert flux==0;continue
                    assert geometry>0
                    expected=sum(g*cellmu[cid][car] for cid,g,_ in part)/geometry
                    edgechecks.append(dict(key=job['key'],edge=e['edge_id'],carrier=name,relative=float(abs(mu/expected-1))))
                    for cid,g,weights in part:
                        if not g:continue
                        for node,w in weights.items():
                            if node not in NODES:continue
                            derivative=g/geometry*w*nodal[node][car][1][pop]
                            value=flux/expected*derivative*D(sign)/D(str(VT))
                            ref[mode,node,block,i]+=value;ref[mode,node,block,j]-=value
                            if node not in (i,j) and value:
                                third.append(dict(key=job['key'],edge=e['edge_id'],input_node=node,carrier=name,reference=str(value)))
            for row in a.rows(path/'rows.csv'):
                if row['direction'].startswith('block_'):continue
                mode,node_s,_=row['direction'].split('_');node=int(node_s);block=row['row_block'];outnode=int(row['row_node'])
                if row['direction_mode']==block:continue
                vals=fluxes[block,outnode];scale=sum(abs(f) for f in vals)
                closure.append(float(abs(sum(vals)-D(row['base_residual']))/scale) if scale else 0.)
                target=ref[mode,node,block,outnode];actual=D(row['analytic_derivative'])
                if not target:
                    assert actual==0,(job['key'],row);continue
                error=float(abs(actual/target-1))
                cross.append(dict(key=job['key'],direction=row['direction'],output_block=block,output_node=outnode,
                    reference=str(target),production=str(actual),relative=error,qualified=error<=1e-5,production_zero=actual==0))
            for row in a.rows(path/'jvp.csv'):
                if not row['direction'].startswith('block_'):continue
                for block in ('psi','phin','phip'):
                    ana=float(row[f'analytic_{block}_norm']);fd=float(row[f'finite_difference_{block}_norm'])
                    error=float(row[f'{block}_relative_error'])*max(1.,fd)/max(ana,fd,1e-300)
                    weak=(row['mode'],block) in (('phin','phip'),('phip','phin'))
                    gated=not weak and float(row['amplitude_V'])<1e-4
                    blocks.append(dict(key=job['key'],input=row['mode'],output=block,step_V=row['amplitude_V'],relative=error,gated=gated,qualified=error<=1e-4))
            print('Analyzed',job['key'],flush=True)
    summary=dict(points=8,cell_checks=len(cellchecks),cell_mu_max_relative=max(r['relative'] for r in cellchecks),
        edge_checks=len(edgechecks),edge_mu_max_relative=max(r['relative'] for r in edgechecks),
        cross_checks=len(cross),cross_failed=sum(not r['qualified'] for r in cross),cross_max_relative=max(r['relative'] for r in cross),
        third_vertex_nonzero_contributions=len(third),base_closure_max_relative=max(closure),
        gated_blocks=sum(r['gated'] for r in blocks),gated_block_failures=sum(r['gated'] and not r['qualified'] for r in blocks),
        gated_block_max_relative=max(r['relative'] for r in blocks if r['gated']))
    for name,rows in [('cell_checks.csv',cellchecks),('edge_checks.csv',edgechecks),('cross_checks.csv',cross),('blocks.csv',blocks),('third_vertex.csv',third)]:a.write_csv(OUT/name,rows)
    a.write(OUT/'fixed_summary.json',summary)
    d.matrix.freeze(OUT/'fixed_analysis_evidence.json',[OUT/'fixed_raw_evidence.json',Path(__file__).resolve(),Path(hp.__file__)]+[OUT/n for n in ('cell_checks.csv','edge_checks.csv','cross_checks.csv','blocks.csv','third_vertex.csv','fixed_summary.json')])
    print(summary,flush=True)
    assert summary['cell_mu_max_relative']<=1e-12 and summary['edge_mu_max_relative']<=1e-12
    assert summary['cross_failed']==0 and summary['gated_block_failures']==0 and summary['base_closure_max_relative']<=1e-10


def solve(arm):
    a.verify(OUT/'fixed_analysis_evidence.json');s=a.read(OUT/'fixed_summary.json')
    assert s['cross_failed']==s['gated_block_failures']==0
    assert s['cell_mu_max_relative']<=1e-12 and s['edge_mu_max_relative']<=1e-12
    q.run.LOCAL=LOCAL/arm;q.run.OUT=OUT/arm
    jobs=[]
    for old in a.read(OUT/'contract.json')['jobs']:
        job=dict(old);cfg=a.read(Path(old['original_config']))
        if arm=='legacy':cfg['solver']['mobility']['edge_averaging']='legacy'
        path=LOCAL/arm/'inputs'/f"{old['case']}_{old['index']}"/old['arm']/'config.json';a.write(path,cfg)
        job['original_config']=str(path);jobs.append(job)
    a.write(q.run.OUT/'validation_contract.json',dict(jobs=jobs))
    d.matrix.freeze(q.run.OUT/'validation_freeze.json',[OUT/'fixed_analysis_evidence.json',q.run.RUNNER]+[Path(j['original_config']) for j in jobs])
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(q.run.solve_job,jobs) for r in group]
    q.run.v.csv_union(q.run.OUT/'attempts.csv',rows)
    d.matrix.freeze(q.run.OUT/'dc_evidence.json',[q.run.OUT/'validation_freeze.json',q.run.OUT/'attempts.csv']+[p for p in (q.run.LOCAL/'dc').rglob('*') if p.is_file()])
    print(arm,'qualified attempts',sum(r['qualified'] for r in rows),'/',len(rows),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','fixed','analyze_fixed','legacy','candidate'))
    action=parser.parse_args().action
    solve(action) if action in ('legacy','candidate') else globals()[action]()
