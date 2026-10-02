"""Read-only fixed-state mobility audit. No Newton solve or threshold changes.

Cell probes consume rounded physical fields. Their edge reconstruction is
compared separately with the frozen complete-split production edge ledger.
Native node-to-edge averaging is deliberately labelled incommensurate.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def stats(values):
    values = sorted(abs(float(x)) for x in values)
    require(bool(values) and all(math.isfinite(x) for x in values), 'Empty/nonfinite statistic')
    return dict(count=len(values), median=values[len(values)//2],
                p95=values[min(len(values)-1, int(.95*len(values)))], max=values[-1])


def scalar(path, key='node_id'):
    data = rows(path)
    result = {int(r[key]): float(r['component0']) for r in data}
    require(len(result)==len(data), 'Duplicate scalar indices')
    require(all(math.isfinite(x) for x in result.values()), 'Nonfinite scalar field')
    return result


def cell_data(path, mesh):
    data = rows(path); groups = defaultdict(list)
    for r in data:
        groups[int(r['cell_id'])].append(r)
    expected = {int(c['id']): c for c in mesh['triangles']}
    require(set(groups)==set(expected), 'Incomplete cell probe')
    silicon = {int(r['id']) for r in mesh['regions'] if r['material']=='Si'}
    cells = {}; parts = defaultdict(list)
    for cid, group in groups.items():
        group.sort(key=lambda r:int(r['local_vertex']))
        require([int(r['local_vertex']) for r in group]==[0,1,2], 'Cell local index mismatch')
        ns = [int(r['node_id']) for r in group]
        require(ns==expected[cid]['node_ids'], 'Cell topology mismatch')
        if expected[cid]['region_id'] not in silicon:
            continue
        mu = [float(group[0][f'{c}_mobility_m2_V_s']) for c in ('electron','hole')]
        require(all(math.isfinite(x) and x>0 for x in mu), 'Invalid Si mobility')
        for r in group:
            require([float(r[f'{c}_mobility_m2_V_s']) for c in ('electron','hole')]==mu, 'Inconsistent repeated cell mobility')
        measure = [float(r['measure_m2']) for r in group]
        require(sum(measure)>0, 'Nonpositive cell measure sum')
        cells[cid] = dict(nodes=ns, mu=mu, measure=measure)
        for k,r in enumerate(group):
            g = float(r['coefficient_next'])
            require(math.isfinite(g) and g>=0, 'Invalid transferred box coefficient')
            parts[tuple(sorted((ns[k],ns[(k+1)%3])))].append((cid,g))
    return cells, dict(parts)


def edge_reconstruct(parts, cell_values):
    result = {}
    for edge, terms in parts.items():
        denom=math.fsum(g for _,g in terms)
        if denom>0:
            result[edge]=math.fsum(g*cell_values[cid] for cid,g in terms)/denom
    return result


def probe(runner, cfg, target, stem):
    config=target/(stem+'.json'); output=target/(stem+'.csv')
    cfg=dict(cfg); cfg.pop('output_state_file',None); cfg.pop('sweep',None)
    cfg.update(simulation_type='element_box_probe',output_csv=str(output))
    write(config,cfg)
    status_path=target/(stem+'.status.json')
    identity=dict(config=sha(config),state=sha(Path(cfg['state_file'])),binary=sha(runner))
    if status_path.exists():
        old=read(status_path);require(old['identity']==identity and old['exit']==0,'Cached probe identity changed')
        require(old['output']==sha(output),'Cached probe output changed'); return
    env=os.environ.copy();env.update(VELA_LINEAR_SOLVER='sparselu',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VELA_LINEAR_THREADS='1')
    with (target/(stem+'.stdout')).open('w') as out,(target/(stem+'.stderr')).open('w') as err:
        proc=subprocess.run([str(runner),'--config',str(config),'--log','off'],env=env,stdout=out,stderr=err)
    status=dict(identity=identity,exit=proc.returncode)
    if proc.returncode==0:
        record=json.loads((target/(stem+'.stdout')).read_text().strip().splitlines()[-1])
        require(record['state_evaluation']=='physical_double_fields' and record['split_low_components_used'] is False,'Probe evaluation semantics missing')
        status.update(output=sha(output),record=record)
    write(status_path,status);require(proc.returncode==0,'Cell probe failed: '+str(target))
    require(identity['state']==sha(Path(cfg['state_file'])),'Input state mutated')


def audit_point(point, target, inputs, runner, native=None):
    cfg=read(point/'config.json'); case=point.parent.name; index=int(point.name[3:]); device=case.split('_')[0]
    require(read(point/'comparison.json')['passed'] is True,'Unqualified source point')
    cfg['state_file']=str((point/'state.h5').resolve());cfg['state_format']='hdf5'
    for name in ('mesh_file','node_doping_file','materials_file'):
        cfg[name]=str((inputs/(device+'/mesh.json' if name=='mesh_file' else device+'/doping.csv' if name=='node_doping_file' else 'materials.json')).resolve())
    mesh=read(inputs/device/'mesh.json')
    probe(runner,cfg,target,'vela_cells')
    cells,parts=cell_data(target/'vela_cells.csv',mesh)
    ledger={tuple(sorted((int(r['node0']),int(r['node1'])))):r for r in rows(point/'acceptance_edges.csv')}
    original=point/'native_export'; result=dict(case=case,index=index,vg=index*.05,cells=len(cells),
        state_sha256=sha(point/'state.h5'),edge_sha256=sha(point/'acceptance_edges.csv'),carriers={})
    for b,car in enumerate(('electron','hole')):
        projected=edge_reconstruct(parts,{cid:c['mu'][b] for cid,c in cells.items()})
        node=scalar(original/'fields'/('eMobility_region0.csv' if b==0 else 'hMobility_region0.csv'))
        projection=[];wrong=[];weighted=0.;weight=0.;worst=None
        for edge,value in projected.items():
            production=float(ledger[edge][car+'_mobility_m2_V_s'])
            require(production>0,'Invalid frozen production mobility')
            error=value/production-1;projection.append(error)
            if worst is None or abs(error)>abs(worst['relative']):worst=dict(edge=edge,relative=error)
            naive=(node[edge[0]]+node[edge[1]])*.5*1e-4
            wrong.append(naive/production-1)
            w=abs(float(ledger[edge][car+'_flux']));weighted+=w*abs(naive/production-1);weight+=w
        result['carriers'][car]=dict(projected_cell_to_split_edge=stats(projection),projection_worst=worst,
            incommensurate_node_average_to_edge=stats(wrong),
            incommensurate_flux_weighted_relative=weighted/weight if weight else None)
    if native is not None:
        seed_cfg=dict(cfg);seed_cfg['state_file']=str((point/'native_seed.h5').resolve())
        probe(runner,seed_cfg,target,'native_seed_cells')
        seed_cells,seed_parts=cell_data(target/'native_seed_cells.csv',mesh)
        require(seed_parts==parts,'State changed mesh geometry')
        # The importer must map original node coordinates and topology exactly.
        for name in ('nodes.csv','elements.csv'):
            require(rows(native/name)==rows(original/name),'Native Load/Plot topology changed: '+name)
        identity={}
        for name in ('ElectrostaticPotential','eQuasiFermiPotential','hQuasiFermiPotential','eDensity','hDensity'):
            a=scalar(original/'fields'/(name+'_region0.csv'));z=scalar(native/'fields'/(name+'_region0.csv'))
            require(set(a)==set(z),'Native loaded field coverage')
            density=name.endswith('Density')
            delta=max(abs(z[k]/a[k]-1) if density else abs(z[k]-a[k]) for k in a)
            identity[name]=delta;require(delta<=(1e-10 if density else 1e-12),'Native Load changed state: '+name)
        result['native_load_identity']=identity
        # Enormal also consumes the dielectric-side electrostatic potential.
        for region in (1,2,3):
            name=f'ElectrostaticPotential_region{region}.csv'
            a=scalar(original/'fields'/name); z=scalar(native/'fields'/name)
            require(set(a)==set(z),'Dielectric potential coverage')
            delta=max(abs(z[k]-a[k]) for k in a)
            identity[name]=delta
            require(delta<=1e-12,'Native Load changed dielectric potential')
        element_lookup={tuple(sorted(int(r[k]) for k in ('node0','node1','node2'))):int(r['id']) for r in rows(native/'elements.csv') if r['material']=='Si'}
        require(len(element_lookup)==len(cells),'Native cell coverage')
        for b,car in enumerate(('electron','hole')):
            native_mu=scalar(native/'fields'/('eMobility_region0_cells.csv' if b==0 else 'hMobility_region0_cells.csv'),'cell_id')
            require(set(native_mu)==set(element_lookup.values()),'Native mobility cell coverage')
            mapped={cid:native_mu[element_lookup[tuple(sorted(c['nodes']))]]*1e-4 for cid,c in cells.items()}
            require(all(x>0 for x in mapped.values()),'Nonpositive native mobility')
            edge=edge_reconstruct(parts,mapped)
            result['carriers'][car].update(
                vela_state_cell_vs_native_state_cell=stats(c['mu'][b]/mapped[cid]-1 for cid,c in cells.items()),
                native_seed_cell_vs_native_cell=stats(seed_cells[cid]['mu'][b]/mapped[cid]-1 for cid in cells),
                vela_cell_vs_native_seed_cell=stats(c['mu'][b]/seed_cells[cid]['mu'][b]-1 for cid,c in cells.items()),
                same_export_cell_gate_1e_7_passed=all(abs(seed_cells[cid]['mu'][b]/mapped[cid]-1)<=1e-7 for cid in cells),
                split_edge_vs_native_cell_reconstruction=stats(float(ledger[e][car+'_mobility_m2_V_s'])/mu-1 for e,mu in edge.items()))
            cid=max(cells,key=lambda i:abs(seed_cells[i]['mu'][b]/mapped[i]-1))
            result['carriers'][car]['same_export_worst_cell']=dict(cell=cid,nodes=cells[cid]['nodes'],
                coordinates_um=[mesh['nodes'][i] for i in cells[cid]['nodes']],
                relative=seed_cells[cid]['mu'][b]/mapped[cid]-1,
                native_m2_V_s=mapped[cid],vela_m2_V_s=seed_cells[cid]['mu'][b])
        # Node-display reconstruction hypotheses, all retained, no best-fit promotion.
        for b,car in enumerate(('electron','hole')):
            nodes=scalar(original/'fields'/('eMobility_region0.csv' if b==0 else 'hMobility_region0.csv'))
            native_mu=scalar(native/'fields'/('eMobility_region0_cells.csv' if b==0 else 'hMobility_region0_cells.csv'),'cell_id')
            modes={}
            for mode in ('equal','cell_area','box_vertex'):
                numer=defaultdict(float);denom=defaultdict(float)
                for cid,c in cells.items():
                    mu=native_mu[element_lookup[tuple(sorted(c['nodes']))]]
                    xy=[mesh['nodes'][i] for i in c['nodes']]
                    area=abs((xy[1]['x']-xy[0]['x'])*(xy[2]['y']-xy[0]['y'])-(xy[1]['y']-xy[0]['y'])*(xy[2]['x']-xy[0]['x']))*.5e-12
                    for k,n in enumerate(c['nodes']):
                        w=1. if mode=='equal' else area if mode=='cell_area' else c['measure'][k]
                        numer[n]+=w*mu;denom[n]+=w
                modes[mode]=stats(numer[n]/denom[n]/nodes[n]-1 for n in nodes)
            result['carriers'][car]['native_node_display_hypotheses']=modes
    write(target/'result.json',result);return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--runner',type=Path,required=True);parser.add_argument('--cases',nargs='+',required=True)
    parser.add_argument('--indices',nargs='+',type=int,default=list(range(51)))
    parser.add_argument('--native-exports',type=Path)
    args=parser.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    results=[];start=time.time()
    for case in args.cases:
        for index in args.indices:
            point=args.base/'matrix'/case/f'vg_{index:03d}'
            native=args.native_exports/case/f'vg_{index:03d}' if args.native_exports else None
            result=audit_point(point,args.output/case/f'vg_{index:03d}',args.base/'data/inputs',args.runner,native)
            results.append(result);write(args.output/'progress.json',dict(points=len(results),last=dict(case=case,index=index)))
            print(case,index,flush=True)
    write(args.output/'summary.json',dict(points=len(results),results=results,elapsed_seconds=time.time()-start,
        new_dc_solves=0,gate_changes=False,native_cells_available=args.native_exports is not None,
        scope='Fixed-state field comparison; projected cell probe and full-split edge values separated'))


if __name__=='__main__':
    main()
