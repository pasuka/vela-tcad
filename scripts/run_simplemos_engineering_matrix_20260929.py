"""Gated cold/native initialization and physical-field matrix on T470p.

Every native seed is independently exported, then self-consistently solved.
Rounded native fields are never used for a fixed-state current attribution.
"""
import argparse
import copy
from decimal import Decimal, localcontext
import math
from pathlib import Path
import subprocess
import numpy as np

import simplemos_hfs_cloud_20260926 as h
import state_archive
from simplemos_original_matrix_cloud_20260926 import Matrix
from check_simplemos_cold_overlap_20260927 import load, potential
from simplemos_srh_cloud_20260926 import signed_si
from sentaurus_import import parse_quoted_list, parse_values_block


def scalar(folder,name,all_regions=False):
    files=sorted((folder/'fields').glob(name+'_region*.csv')) if all_regions else [folder/'fields'/f'{name}_region0.csv']
    values={}
    for path in files:
        for row in h.rows(path):
            i,v=int(row['node_id']),float(row['component0'])
            if i in values and values[i]!=v:raise ValueError(f'Conflicting regional field {path}:{i}')
            if not math.isfinite(v):raise ValueError('Nonfinite native field')
            values[i]=v
    return values


def native_seed(export,mesh,output):
    nodes=h.rows(export/'nodes.csv');n=len(mesh['nodes'])
    if len(nodes)!=n:raise ValueError('Node count mismatch')
    for row,node in zip(nodes,mesh['nodes']):
        if int(row['id'])!=node['id'] or abs(float(row['x_um'])-node['x'])>1e-12 or abs(float(row['y_um'])-node['y'])>1e-12:
            raise ValueError('Native mesh requires a verified coordinate map')
    fields={name:np.zeros(n) for name in ('psi','phin','phip','electrons_m3','holes_m3')}
    for name,native in [('psi','ElectrostaticPotential'),('phin','eQuasiFermiPotential'),('phip','hQuasiFermiPotential'),('electrons_m3','eDensity'),('holes_m3','hDensity')]:
        values=scalar(export,native,name=='psi')
        if name=='psi' and len(values)!=n:raise ValueError('Missing oxide potential')
        for i,v in values.items():fields[name][i]=v*(1e6 if name.endswith('_m3') else 1.)
    meta=dict(mode='dd',mesh_sha256=state_archive.mesh_identity(mesh,1e-6),potential_origin_V=0.,
              provenance='M60 native exported binary64 fields; independent Newton seed, not a full-precision native state')
    if output.exists():
        old,m=state_archive.read(output,n,meta['mesh_sha256'])
        if m!=meta or any(not np.array_equal(old[k],v) for k,v in fields.items()):raise ValueError('Seed identity changed')
    else:state_archive.write(output,fields,meta)


def qualify_native(matrix,c,vg,seed,dest):
    cfg=matrix.config(c);cfg.pop('sweep',None)
    for contact in cfg['contacts']:contact['bias']=c['vd'] if contact['name']=='drain' else vg if contact['name']=='gate' else 0.
    cfg.update(simulation_type='newton_solve_from_state',state_format='hdf5',state_file=str(seed),output_state_file=str(dest/'state.h5'))
    h.write(dest/'config.json',cfg);status=matrix.executor.execute(dest/'config.json')
    if status['exit_code']!=0 or not status.get('converged'):raise ValueError(f'Native reclosure failed: {c["case"]}/{vg}')
    probe=copy.deepcopy(cfg);probe.pop('output_state_file');probe.update(simulation_type='newton_carrier_term_probe',state_file=str(dest/'state.h5'),output_csv=str(dest/'all_row.csv'),carrier_term_probe=dict(solved_equation_terms=True))
    probe['solver']['stable_merit_comparison']=False;probe['solver']['carrier_row_convergence']['mode']='report'
    probe['solver']['global_continuity_closure']=dict(mode='enforce',tolerance=1e-6,source_floor=1e-10)
    edge=copy.deepcopy(probe);edge.pop('carrier_term_probe');edge.update(simulation_type='sg_edge_flux_probe',output_csv=str(dest/'acceptance_edges.csv'))
    functional=copy.deepcopy(edge);functional.pop('output_csv');functional.update(simulation_type='terminal_current_functional_probe',contact='drain',residual_output_csv=str(dest/'port_residual.csv'),contact_edge_output_csv=str(dest/'port_edges.csv'))
    for name,deck in [('all_row',probe),('acceptance_edges',edge),('functional',functional)]:
        h.write(dest/(name+'.json'),deck);matrix.executor.execute(dest/(name+'.json'))
    result=h.qualify(dest,h.read(matrix.inputs/c['device']/'geometry.json'))
    result.update(case=c['case'],device=c['device'],vd=c['vd'],vg=vg,index=round(vg/.05),seed_sha256=h.sha(seed),state_sha256=h.sha(dest/'state.h5'))
    h.write(dest/'result.json',result)
    if not result['qualified']:raise ValueError(f'Native all-row/port qualification failed: {result}')
    return result


def compare_pair(cold_state,warm_state,cold,warm,geo):
    a,b=load(cold_state),load(warm_state)
    if a[1]['mesh_sha256']!=b[1]['mesh_sha256']:raise ValueError('Dual mesh mismatch')
    with localcontext() as ctx:
        ctx.prec=100
        delta={f+'_max_V':max(float(abs(potential(a,f,i)-potential(b,f,i))) for i in geo['free_si']) for f in ('psi','phin','phip')}
    density=max(abs(float(a[0][f][i])/float(b[0][f][i])-1) for f in ('electrons_m3','holes_m3') for i in geo['free_si'])
    current=abs(cold['current_A_per_um']/warm['current_A_per_um']-1)
    ok=cold['qualified'] and warm['qualified'] and all(math.isfinite(v) and v<=1e-6 for v in delta.values()) and math.isfinite(density) and density<=1e-4 and math.isfinite(current) and current<=1e-6
    return dict(**delta,density_relative=density,Id_relative=current,qualified=ok)


def fields(dest,export,mesh,geo,cfg):
    state=load(dest/'state.h5');terms=h.ordered(dest/'all_row.csv',geo['count']);ids=geo['all_si']
    area=signed_si(mesh);weights=np.zeros(len(mesh['nodes']))
    silicon={r['id'] for r in mesh['regions'] if r['material']=='Si'}
    for cell in mesh['triangles']:
        if cell['region_id'] not in silicon:continue
        a,b,c=(mesh['nodes'][i] for i in cell['node_ids'])
        measure=abs((b['x']-a['x'])*(c['y']-a['y'])-(b['y']-a['y'])*(c['x']-a['x']))/6
        for i in cell['node_ids']:weights[i]+=measure
    rows=[];rates={};density_errors=[]
    native_names={'psi':'ElectrostaticPotential','phin':'eQuasiFermiPotential','phip':'hQuasiFermiPotential','electrons_m3':'eDensity','holes_m3':'hDensity'}
    with localcontext() as ctx:
        ctx.prec=100;D=lambda x:Decimal.from_float(float(x));vt=D(state[1]['packed_potential_scale_V'])
        for i in ids:
            psi,fn,fp=(potential(state,k,i) for k in ('psi','phin','phip'))
            ni=D(terms[i]['ni_eff_m3']);n=ni*((psi-fn)/vt).exp();p=ni*((fp-psi)/vt).exp()
            density_errors.extend([abs(float(n)*1e6/float(state[0]['electrons_m3'][i])-1),abs(float(p)*1e6/float(state[0]['holes_m3'][i])-1)])
            doping=float(terms[i]['donors_m3'])+float(terms[i]['acceptors_m3']);tau={}
            for car in ('electron','hole'):
                par=cfg['solver']['srh_doping_dependence'][car]
                tau[car]=D(par['tau_min_s']+(par['tau_max_s']-par['tau_min_s'])/(1+(doping/par['reference_doping_m3'])**par['gamma']))
            rates[i]=float(ni*ni*(((fp-fn)/vt).exp()-1)/(tau['hole']*(n+ni)+tau['electron']*(p+ni)))
        for field,native_name in native_names.items():
            native=scalar(export,native_name)
            if set(native)!=set(ids):raise ValueError('Native Si mask mismatch')
            delta=[math.log10(float(state[0][field][i])*1e-6/native[i]) if field.endswith('_m3') else float(potential(state,field,i)-D(native[i])) for i in ids]
            rows.append(dict(field=field,unit='dex' if field.endswith('_m3') else 'V',max_abs=max(map(abs,delta)),weighted_rms=math.sqrt(math.fsum(weights[i]*v*v for i,v in zip(ids,delta))/math.fsum(weights[i] for i in ids))))
    native=scalar(export,'srhRecombination');denom=math.fsum(weights[i]*abs(native[i]) for i in ids)
    rate_error=math.fsum(weights[i]*abs(rates[i]-native[i]) for i in ids)/denom if denom else None
    edges=h.rows(dest/'acceptance_edges.csv');edge=max(edges,key=lambda r:abs(float(r['electron_flux'])))
    factor=float(edge['electron_particle_line_flux_per_m_s'])*1.602176634e-19*1e-6/float(edge['electron_flux'])
    actual=[float(terms[i]['electron_recombination'])*factor for i in ids]
    expected=[rates[i]*area[i]*1.602176634e-19 for i in ids]
    source_error=math.fsum(abs(a-b) for a,b in zip(actual,expected))/math.fsum(map(abs,actual))
    return dict(fields=rows,SRH_weighted_L1=rate_error,density_reconstruction=max(density_errors),source_reconstruction=source_error,
                reconstruction_passed=max(density_errors)<=1e-12 and source_error<=1e-7,
                native_field_acceptance='descriptive_only_no_new_threshold',mobility_comparison='not inferred from node-versus-edge averages')


def run(base,case_name,runner=None,importer=None):
    base=base.resolve();root=base/'matrix';root.mkdir(exist_ok=True)
    if not h.read(base/'preflight/summary.json')['passed']:raise ValueError('Platform preflight not qualified')
    runner=Path(runner).resolve() if runner else base/'source/build-release/vela_example_runner.exe';inputs=base/'data/inputs'
    importer=Path(importer).resolve() if importer else base/'source/build-release/sentaurus_import.exe'
    for rel,expected in h.read(inputs/'hashes.json').items():
        if h.sha(inputs/rel)!=expected:raise ValueError(f'Input changed: {rel}')
    c=next(c for c in h.read(inputs/'contract.json')['cases'] if c['case']==case_name)
    native_dir=base/'native/bundle'/c['device'];prefix='m60fields_'+case_name
    if (base/'native'/(prefix+'.exitcode')).read_text().strip()!='0':raise ValueError('Native run incomplete')
    tdrs=sorted(native_dir.glob(prefix+'_state_*_des.tdr'))
    if len(tdrs)!=51:raise ValueError(f'Expected 51 native states, found {len(tdrs)}')
    plt=native_dir/('IdVg_'+prefix+'_des.plt');text=plt.read_text();names=parse_quoted_list(text,'datasets')
    native_points=[dict(zip(names,r)) for r in parse_values_block(text,len(names))]
    if len(native_points)!=51 or any(abs(r['gate OuterVoltage']-i*.05)>1e-10 or abs(r['drain OuterVoltage']-c['vd'])>1e-10 for i,r in enumerate(native_points)):
        raise ValueError('Native curve bias grid mismatch')
    seal=dict(binary=h.sha(runner),importer=h.sha(importer),driver=h.sha(__file__),profile=h.sha(base/'data/engineering_contract.json'),inputs=h.sha(inputs/'hashes.json'),native_PLT=h.sha(plt),native={p.name:h.sha(p) for p in tdrs})
    dest=root/case_name;dest.mkdir(exist_ok=True)
    if (dest/'seal.json').exists():
        if h.read(dest/'seal.json')!=seal:raise ValueError('Matrix identity changed')
    else:h.write(dest/'seal.json',seal)
    # Reuse qualified sweep and all-row instrumentation without its obsolete
    # original-M8 current acceptance as the engineering M60 acceptance.
    matrix=object.__new__(Matrix);matrix.root=root/'cold';matrix.inputs=inputs;matrix.runner=runner;matrix.contract=h.read(inputs/'contract.json')
    matrix.executor=object.__new__(h.Run);matrix.executor.runner=runner
    gate=matrix.root/case_name/'gate'
    if not (dest/'cold_results.json').exists():
        cold=matrix.case(c)
        if len(cold)!=51 or not all(r['qualified'] for r in cold):raise ValueError('Cold numerical qualification failed')
        h.write(dest/'cold_results.json',cold)
    cold=h.read(dest/'cold_results.json');mesh=h.read(inputs/c['device']/'mesh.json');geo=h.read(inputs/c['device']/'geometry.json')
    references={int(r['index']):float(r['tight_default_Id_A_per_um']) for r in h.rows(base/'data/m60_joined.csv') if r['case']==case_name}
    if set(references)!=set(range(51)):raise ValueError('Missing M60 reference')
    results=[]
    for index,tdr in enumerate(tdrs):
        target=dest/f'vg_{index:03d}';target.mkdir(exist_ok=True);export=target/'native_export'
        if not (export/'field_manifest.json').exists():
            export.mkdir(exist_ok=True)
            with (target/'import.stdout.json').open('w') as out,(target/'import.stderr.txt').open('w') as err:
                subprocess.run([str(importer),'--tdr',str(tdr),'--export-dir',str(export)],stdout=out,stderr=err,check=True)
        seed=target/'native_seed.h5';native_seed(export,mesh,seed)
        warm=qualify_native(matrix,c,index*.05,seed,target)
        old=next(r for r in cold if r['index']==index)
        pair=compare_pair(gate/('state_bias_'+format(index*.05,'.6f').replace('.','p')+'.h5'),target/'state.h5',old,warm,geo)
        field=fields(target,export,mesh,geo,matrix.config(c));h.write(target/'fields.json',field)
        errors={arm:100*(row['current_A_per_um']/references[index]-1) for arm,row in [('cold',old),('native',warm)]}
        record=dict(case=case_name,index=index,vg=index*.05,dual=pair,Id_error_percent=errors,reconstruction_passed=field['reconstruction_passed'],
                    supplemental_native_Id_change_percent=100*(native_points[index]['drain TotalCurrent']/references[index]-1))
        record['passed']=pair['qualified'] and field['reconstruction_passed'] and all(math.isfinite(v) and abs(v)<=2 for v in errors.values())
        h.write(target/'comparison.json',record);results.append(record);h.write(dest/'progress.json',dict(points=len(results),last=record))
        print(case_name,index,'passed',record['passed'],flush=True)
        if not record['passed']:raise ValueError(f'Matrix gate failed: {record}')
    h.write(dest/'summary.json',dict(points=len(results),passed=len(results)==51 and all(r['passed'] for r in results),results=results))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--base',type=Path,required=True);parser.add_argument('--case',required=True);parser.add_argument('--runner',type=Path);parser.add_argument('--importer',type=Path);args=parser.parse_args()
    try:run(args.base,args.case,args.runner,args.importer)
    except Exception as exc:
        h.write(args.base/'matrix'/args.case/'failure.json',dict(error=repr(exc)));raise
