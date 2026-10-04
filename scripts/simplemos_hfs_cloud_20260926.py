"""Portable merged-source HFS baseline; big states and analysis stay on compute.

No numerical settings are tuned here. A sealed input package and binary define
the experiment. Controls gate curves; every failed/reloaded attempt is retained.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import csv
from decimal import Decimal, localcontext
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time
import state_archive
import migrate_state_seed_to_hdf5


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Encode failed nonfinite metrics as null, never as a JSON NaN success.
    def clean(x):
        if isinstance(x, float) and not math.isfinite(x): return None
        if isinstance(x, dict): return {k:clean(v) for k,v in x.items()}
        if isinstance(x, (tuple,list)): return [clean(v) for v in x]
        return x
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(clean(obj), indent=2, allow_nan=False)+'\n', encoding='utf-8')
    tmp.replace(path)


def rows(path):
    with Path(path).open(encoding='utf-8', newline='') as f: return list(csv.DictReader(f))


def csvout(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    keys=list(dict.fromkeys(k for row in data for k in row))
    with path.open('w', encoding='utf-8', newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(data)


def ordered(path, count):
    result=rows(path)
    assert len(result)==count, (path,len(result),count)
    result.sort(key=lambda r:int(r['node_id']))
    assert [int(r['node_id']) for r in result]==list(range(count)), path
    return result


def qualify(dest, geo):
    """Same frozen all-row, global-source, KCL and port equations as 2026-09-14."""
    statuses={name:read(dest/(name+'.status.json')) for name in
              ('config','all_row','acceptance_edges','functional')}
    terms=ordered(dest/'all_row.csv',geo['count']);edges=rows(dest/'acceptance_edges.csv')
    ratios=[];bad_rows=[];closure={};contacts=set(geo['contacts'])
    for i in geo['free_si']:
        r=terms[i]
        for car in ('electron','hole'):
            numbers=[float(r[car+'_'+k]) for k in ('flux_abs_sum','recombination','impact','residual')]
            scale=max(numbers[0],abs(numbers[1]),abs(numbers[2]))
            ratio=abs(numbers[3])/scale if scale>0 and all(map(math.isfinite,numbers)) else math.inf
            ratios.append(ratio)
            if not math.isfinite(ratio) or ratio>1e-6:
                bad_rows.append(dict(node_id=i,carrier=car,ratio=ratio))
    assert len(ratios)==2*len(geo['free_si']) and len(ratios)>0
    for car in ('electron','hole'):
        flux=math.fsum((int(int(e['node0']) in contacts)-int(int(e['node1']) in contacts))*float(e[car+'_flux']) for e in edges)
        source=math.fsum(float(r[car+'_recombination'])+float(r[car+'_impact']) for i,r in enumerate(terms) if i not in contacts)
        ratio=abs(flux-source)/max(abs(flux),abs(source),1e-10)
        closure[car]=dict(flux=flux,source=source,ratio=ratio,
             active=abs(source)>=1e-10,satisfied=all(map(math.isfinite,(flux,source,ratio))) and (abs(source)<1e-10 or ratio<=1e-6))
    status=statuses['config'];cc=status.get('contact_currents_A_per_um',{})
    current=cc.get('drain',math.nan)
    kcl=abs(math.fsum(cc.values()))/abs(current) if current else math.inf
    port=statuses['functional'];other=port.get('contact_current_extractor_A_per_um',0.)
    error=abs(port.get('current_A_per_um',math.nan)/other-1) if other else math.inf
    ok=(all(s.get('exit_code')==0 for s in statuses.values()) and status.get('converged',False)
        and not bad_rows and all(x['satisfied'] for x in closure.values())
        and math.isfinite(kcl) and kcl<=1e-8 and math.isfinite(error) and error<=1e-8)
    return dict(qualified=bool(ok),current_A_per_um=current,max_row_ratio=max(ratios),
                carrier_rows=len(ratios),row_violations=len(bad_rows),kcl_over_Id=kcl,
                port_relative=error,iterations=status.get('iterations'),
                failure=status.get('failure_reason',''),closure=closure,bad_rows=bad_rows)


def physical(row,field):
    # Keep the exact historical dual-initialization comparison contract.
    car={'phin':'electron','phip':'hole'}.get(field)
    if car and car+'_qf_reference_V' in row:
        return Decimal(row[car+'_qf_reference_V'])+Decimal(row[car+'_qf_increment_V'])
    return Decimal(row[field])


def dual(first,second,geo):
    x=ordered(Path(first['dest'])/'state.csv',geo['count'])
    y=ordered(Path(second['dest'])/'state.csv',geo['count'])
    with localcontext() as ctx:
        ctx.prec=100
        diffs={f+'_max_V':max(abs(float(physical(x[i],f)-physical(y[i],f))) for i in geo['free_si'])
               for f in ('psi','phin','phip')}
    density=max(abs(float(x[i][k])/float(y[i][k])-1) for i in geo['free_si'] for k in ('electrons_m3','holes_m3'))
    current=abs(first['current_A_per_um']/second['current_A_per_um']-1)
    ok=first['qualified'] and second['qualified'] and all(math.isfinite(v) and v<=1e-6 for v in diffs.values())
    ok=ok and math.isfinite(density) and density<=1e-4 and math.isfinite(current) and current<=1e-6
    return dict(**diffs,density_max_relative=density,Id_relative=current,qualified=bool(ok))


class Run:
    def __init__(self,root,runner,workers):
        self.root=root.resolve();self.inputs=self.root/'inputs';self.runner=runner.resolve();self.workers=workers
        self.contract=read(self.inputs/'contract.json');self.cases=self.contract['cases']
        self.summary=self.root/'summary';self.summary.mkdir(exist_ok=True)
        self.geo={c['case']:read(self.inputs/c['case']/'geometry.json') for c in self.cases}

    def migrate(self):
        reports=[]
        for c in self.cases:
            parent=self.inputs/c['case'];mesh=parent/'mesh_file.json'
            for source in [*sorted(parent.glob('enormal_*.csv')),*sorted(parent.glob('native/*/state.csv'))]:
                target=self.root/'seeds'/source.relative_to(self.inputs).with_suffix('.h5')
                record=target.with_suffix('.migration.json')
                if record.exists():
                    result=read(record)
                    assert sha(target)==result['output_sha256'] and sha(source)==result['source_sha256']
                else:
                    result=migrate_state_seed_to_hdf5.migrate(source,mesh,target,'unit_scaling',0.)
                    write(record,result)
                reports.append(result)
        assert len(reports)==220
        write(self.summary/'seed_migration.json',dict(states=len(reports),all_decoded_fields_equal=all(r['all_decoded_fields_equal'] for r in reports),reports=reports))
        write(self.root/'seed_manifest.json',{p.relative_to(self.root).as_posix():sha(p) for p in sorted((self.root/'seeds').rglob('*')) if p.is_file()})

    def export_state(self,dest,case):
        geo=self.geo[case['case']]
        expected=state_archive.mesh_identity(read(self.inputs/case['case']/'mesh_file.json'),1e-6)
        fields,meta=state_archive.read(dest/'state.h5',geo['count'],expected)
        table=state_archive.fields_to_rows(fields)
        for row in table:
            for key in ('packed_potential_scale_V','split_state_schema','split_mesh_fingerprint'):
                if key in meta:row[key]=meta[key]
        csvout(dest/'state.csv',table)

    def seal(self):
        for rel,h in read(self.inputs/'manifest.json').items(): assert sha(self.inputs/rel)==h, rel
        for rel,h in read(self.inputs/'source_hashes.json').items(): assert sha(self.root/'source'/rel)==h, rel
        for rel,h in read(self.root/'seed_manifest.json').items(): assert sha(self.root/rel)==h,rel
        expected=dict(binary=str(self.runner),binary_sha256=sha(self.runner),
                      input_manifest_sha256=sha(self.inputs/'manifest.json'),driver_sha256=sha(__file__),
                      seed_manifest_sha256=sha(self.root/'seed_manifest.json'),
                      state_helpers={Path(m.__file__).name:sha(m.__file__) for m in (state_archive,migrate_state_seed_to_hdf5)},
                      solver_environment=dict(VELA_LINEAR_SOLVER='sparselu',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VELA_LINEAR_THREADS='1'))
        target=self.root/'run_seal.json'
        if target.exists(): assert read(target)==expected,'Frozen run identity changed'
        else: write(target,expected)
        fixture=self.inputs/'acceptance_fixture'
        actual=qualify(fixture,self.geo[read(fixture/'case.json')['case']]);historical=read(fixture/'result.json')
        for key in ('qualified','current_A_per_um','max_row_ratio','row_violations','kcl_over_Id','port_relative'):
            assert actual[key]==historical[key],(key,actual[key],historical[key])
        write(self.summary/'harness_validation.json',dict(historical_fixture_equivalent=True,compared_fields=6,source_files=len(read(self.inputs/'source_hashes.json'))))

    def execute(self,path):
        output=path.with_suffix('.status.json')
        if output.exists(): return read(output)
        # A process killed before writing status leaves evidence; do not erase it.
        for suffix in ('.stdout.txt','.stderr.txt'):
            log=path.with_suffix(suffix)
            if log.exists(): log.rename(log.with_name(log.name+'.interrupted-'+str(time.time_ns())))
        env=os.environ.copy()
        for key in list(env):
            if key.startswith(('VELA_DIAGNOSTIC_','VELA_VALIDATE_','VELA_MINORITY_','VELA_SIMPLEMOS_','VELA_TEST_','VELA_NATIVE_GEOMETRY_','VELA_CANDIDATE_')):env.pop(key)
        env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VELA_LINEAR_THREADS='1',VELA_LINEAR_SOLVER='sparselu')
        start=time.monotonic()
        with path.with_suffix('.stdout.txt').open('w') as out,path.with_suffix('.stderr.txt').open('w') as err:
            proc=subprocess.run([str(self.runner),'--config',str(path),'--log','off'],stdout=out,stderr=err,env=env)
        try: status=json.loads(path.with_suffix('.stdout.txt').read_text().strip().splitlines()[-1])
        except (ValueError,IndexError): status=dict(failure_reason='No JSON status; inspect preserved logs')
        status.update(exit_code=proc.returncode,elapsed_seconds=time.monotonic()-start)
        write(output,status);return status

    def attempt(self,case,index,arm,seed,attempt,phase):
        dest=self.root/phase/case['case']/arm/f'vg_{index:03d}'/f'attempt_{attempt}'
        dest.mkdir(parents=True,exist_ok=True)
        cfg=read(self.inputs/case['case']/'template.json')
        for key in ('mesh_file','materials_file','node_doping_file'):cfg[key]=str(self.inputs/cfg[key])
        for contact in cfg['contacts']:
            if contact['name']=='gate':contact['bias']=index/50
        cfg.update(state_file=str(seed),output_state_file=str(dest/'state.h5'))
        probe=copy.deepcopy(cfg);probe.pop('output_state_file');probe['solver'].pop('local_update_diagnostics',None)
        probe.update(simulation_type='newton_carrier_term_probe',state_file=str(dest/'state.h5'),output_csv=str(dest/'all_row.csv'),carrier_term_probe=dict(solved_equation_terms=True))
        probe['solver']['carrier_row_convergence']['mode']='report'
        probe['solver']['global_continuity_closure']=dict(mode='enforce',tolerance=1e-6,source_floor=1e-10)
        probe['solver']['stable_merit_comparison']=False
        edge=copy.deepcopy(probe);edge.pop('carrier_term_probe');edge.update(simulation_type='sg_edge_flux_probe',output_csv=str(dest/'acceptance_edges.csv'))
        fun=copy.deepcopy(edge);fun.pop('output_csv');fun.update(simulation_type='terminal_current_functional_probe',contact='drain',residual_output_csv=str(dest/'port_residual.csv'),contact_edge_output_csv=str(dest/'port_edges.csv'))
        configs=dict(config=cfg,all_row=probe,acceptance_edges=edge,functional=fun)
        for name,deck in configs.items():
            path=dest/(name+'.json')
            if path.exists():assert read(path)==deck,path
            else:write(path,deck)
        identity=dict(seed_sha256=sha(seed),run_seal_sha256=sha(self.root/'run_seal.json'),
                      configs={k:sha(dest/(k+'.json')) for k in configs})
        seal=dest/'input_seal.json'
        if seal.exists():assert read(seal)==identity,seal
        else:write(seal,identity)
        if (dest/'result.json').exists():
            result=read(dest/'result.json')
            for rel,h in result['output_hashes'].items():assert sha(dest/rel)==h,(dest,rel)
            return result
        # Keep any interrupted state before re-executing a DC attempt.
        if not (dest/'config.status.json').exists():
            for name in ('state.h5','state.csv'):
                path=dest/name
                if path.exists():path.rename(dest/(path.stem+'.interrupted-'+str(time.time_ns())+path.suffix))
        status=self.execute(dest/'config.json')
        result=dict(**case,index=index,vg=index/50,arm=arm,attempt=attempt,dest=str(dest),qualified=False,
                    exit_code=status['exit_code'],failure=status.get('failure_reason',''),elapsed_seconds=status['elapsed_seconds'])
        if (dest/'state.h5').exists():
            self.export_state(dest,case)
            for name in ('all_row','acceptance_edges','functional'):self.execute(dest/(name+'.json'))
            try: result.update(qualify(dest,self.geo[case['case']]))
            except (KeyError,ValueError,AssertionError,FileNotFoundError,ZeroDivisionError) as error:
                result.update(qualified=False,qualification_error=repr(error))
        result['output_hashes']={p.name:sha(p) for p in dest.iterdir() if p.is_file() and p.suffix in ('.csv','.h5')}
        write(dest/'result.json',result)
        print(phase,case['case'],index,arm,attempt,'qualified',result['qualified'],'row',result.get('max_row_ratio'),'seconds',round(result['elapsed_seconds'],2),flush=True)
        return read(dest/'result.json')

    def target(self,case,index,arm,seed,phase):
        result=self.attempt(case,index,arm,seed,0,phase)
        if not result['qualified'] and (Path(result['dest'])/'state.h5').exists():
            result=self.attempt(case,index,arm,Path(result['dest'])/'state.h5',1,phase)
        return result

    def controls(self):
        jobs=[(c,i,arm,self.root/'seeds'/c['case']/(f'enormal_{i:03d}.h5' if arm=='enormal' else f'native/vg_{i:03d}/state.h5'))
              for c in self.cases for i in self.contract['controls'] for arm in ('enormal','native')]
        with ThreadPoolExecutor(max_workers=self.workers) as pool:
            results=list(pool.map(lambda j:self.target(*j,'controls'),jobs))
        write(self.summary/'control_selected.json',results)
        compare=[]
        for c in self.cases:
            for i in self.contract['controls']:
                pair=[r for r in results if r['case']==c['case'] and r['index']==i]
                row=dict(**c,index=i,vg=i/50,qualified=False)
                if all(r['qualified'] for r in pair):row.update(dual(*pair,self.geo[c['case']]))
                compare.append(row)
        csvout(self.summary/'control_dual.csv',compare)
        ok=len(results)==32 and all(r['qualified'] for r in results) and all(r['qualified'] for r in compare)
        write(self.summary/'control_summary.json',dict(states=len(results),qualified_states=sum(r['qualified'] for r in results),paired_points=len(compare),qualified_pairs=sum(r['qualified'] for r in compare),all_qualified=ok))
        assert ok,'Control gate failed; curves not launched; all attempts retained'

    def curves(self):
        assert read(self.summary/'control_summary.json')['all_qualified']
        selected=read(self.summary/'control_selected.json')
        def path(job):
            c,arm=job
            if arm=='native':
                for i in range(51):self.target(c,i,arm,self.root/'seeds'/c['case']/f'native/vg_{i:03d}/state.h5','curves')
            else:
                anchor=next(r for r in selected if r['case']==c['case'] and r['index']==40 and r['arm']=='enormal')
                result=self.target(c,40,arm,Path(anchor['dest'])/'state.h5','curves')
                assert result['qualified'],'Merged continuation anchor failed'
                seed=Path(result['dest'])/'state.h5'
                for indices in (range(39,-1,-1),range(41,51)):
                    current=seed
                    for i in indices:
                        result=self.target(c,i,arm,current,'curves')
                        if result['qualified']:current=Path(result['dest'])/'state.h5'
        with ThreadPoolExecutor(max_workers=self.workers) as pool:
            list(pool.map(path,[(c,arm) for c in self.cases for arm in ('continuation','native')]))

    def summarize(self):
        attempts=[read(p) for phase in ('controls','curves') for p in (self.root/phase).glob('*/*/vg_*/attempt_*/result.json')]
        flat=[{k:v for k,v in r.items() if not isinstance(v,(list,dict))} for r in attempts]
        csvout(self.summary/'attempts.csv',flat)
        failures=[r for r in attempts if not r['qualified']]
        write(self.summary/'failures.json',failures)
        native={(r['case'],int(r['index'])):r for r in rows(self.inputs/'native_points.csv')}
        comparisons=[];selected=[]
        for c in self.cases:
            for i in range(51):
                pair=[]
                for arm in ('continuation','native'):
                    group=[r for r in attempts if '/curves/' in r['dest'].replace('\\','/') and r['case']==c['case'] and r['index']==i and r['arm']==arm]
                    if group:
                        chosen=next((r for r in group if r['qualified']),max(group,key=lambda r:r['attempt']))
                        pair.append(chosen);selected.append(chosen)
                if len(pair)!=2:continue
                row=dict(**c,index=i,vg=i/50,qualified=False,native_Id_A_per_um=float(native[c['case'],i]['Id_A_per_um']))
                if all(r['qualified'] for r in pair):
                    row.update(dual(*pair,self.geo[c['case']]))
                    row.update(vela_Id_A_per_um=pair[0]['current_A_per_um'],Id_error_percent=100*(pair[0]['current_A_per_um']/row['native_Id_A_per_um']-1))
                comparisons.append(row)
        csvout(self.summary/'curve_comparison.csv',comparisons)
        complete=len(selected)==408 and all(r['qualified'] for r in selected) and len(comparisons)==204 and all(r['qualified'] for r in comparisons)
        result=dict(attempts=len(attempts),failed_attempts=len(failures),curve_selected_states=len(selected),curve_qualified_states=sum(r['qualified'] for r in selected),paired_points=len(comparisons),qualified_pairs=sum(r['qualified'] for r in comparisons),numerical_baseline_complete=complete)
        write(self.summary/'progress.json',result)
        if complete:
            field_analysis(self,selected)
            result.update(field_analysis_complete=True,max_abs_Id_error_percent=max(abs(r['Id_error_percent']) for r in comparisons),points_over_1percent=sum(abs(r['Id_error_percent'])>1 for r in comparisons))
            write(self.summary/'progress.json',result)
        write(self.summary/'artifact_hashes.json',{p.relative_to(self.root).as_posix():sha(p) for p in sorted(self.root.rglob('*')) if p.is_file() and p!=self.summary/'artifact_hashes.json' and ('source' not in p.relative_to(self.root).parts) and p.name not in ('run.log','launcher.log')})
        print(json.dumps(result),flush=True)
        return complete


def field_analysis(run,selected):
    """Same Decimal100 density/SRH reconstruction and common Si-volume metrics."""
    fields=[];sources=[];checks=[]
    names=('ElectrostaticPotential','eQuasiFermiPotential','hQuasiFermiPotential','eDensity','hDensity','srhRecombination')
    D=lambda v:Decimal.from_float(float(v))
    for r in selected:
        if r['arm']!='continuation':continue
        dest=Path(r['dest']);geo=run.geo[r['case']];ids=geo['all_si']
        state=ordered(dest/'state.csv',geo['count']);terms=ordered(dest/'all_row.csv',geo['count']);cfg=read(dest/'config.json')
        native=read(run.inputs/r['case']/f"native/vg_{r['index']:03d}/fields.json")
        values={n:[] for n in names};density_errors=[];node_rows=[]
        with localcontext() as ctx:
            ctx.prec=100
            for i in ids:
                s,t=state[i],terms[i];scale=D(s['packed_potential_scale_V']);vt=D(run.contract['VT']);assert scale==vt
                psi=scale*(D(s['packed_psi'])+D(s['packed_psi_low']))
                fn=D(s['electron_qf_reference_V'])+scale*(D(s['packed_electron_qf_increment'])+D(s['packed_electron_qf_increment_low']))
                fp=D(s['hole_qf_reference_V'])+scale*(D(s['packed_hole_qf_increment'])+D(s['packed_hole_qf_increment_low']))
                ni=D(t['ni_eff_m3']);en=ni*((psi-fn)/vt).exp();hp=ni*((fp-psi)/vt).exp()
                density_errors.extend((abs(float(en)*1e6/float(s['electrons_m3'])-1),abs(float(hp)*1e6/float(s['holes_m3'])-1)))
                total=float(t['donors_m3'])+float(t['acceptors_m3']);tau={}
                for car in ('electron','hole'):
                    p=cfg['solver']['srh_doping_dependence'][car]
                    tau[car]=D(p['tau_min_s']+(p['tau_max_s']-p['tau_min_s'])/(1+(total/p['reference_doping_m3'])**p['gamma']))
                rate=ni*ni*(((fp-fn)/vt).exp()-1)/(tau['hole']*(en+ni)+tau['electron']*(hp+ni))
                node=dict(node_id=i)
                for name,val in zip(names,(psi,fn,fp,en,hp,rate)):
                    values[name].append(float(val));node[name+'_vela']=float(val);node[name+'_native']=native[name][str(i)]
                node_rows.append(node)
        # Large node tables remain remote, outside the downloaded summary directory.
        csvout(run.root/'field_nodes'/r['case']/f"vg_{r['index']:03d}.csv",node_rows)
        weights=[geo['barycentric_si'][i] for i in ids];assert min(weights)>0
        edges=rows(dest/'acceptance_edges.csv');edge=max(edges,key=lambda e:abs(float(e['electron_flux'])))
        factor=float(edge['electron_particle_line_flux_per_m_s'])*1.602176634e-19*1e-6/float(edge['electron_flux'])
        actual=[float(terms[i]['electron_recombination'])*factor for i in ids]
        expected=[v*geo['all_cell'][i]*1.602176634e-19 for i,v in zip(ids,values['srhRecombination'])]
        denom=math.fsum(map(abs,actual))
        srherr=math.fsum(abs(x-y) for x,y in zip(actual,expected))/denom if denom else max(map(abs,expected))
        label={k:r[k] for k in ('case','device','vd','vg','index')}
        checks.append(dict(**label,density_reconstruction_max=max(density_errors),SRH_integral_reconstruction_L1=srherr,qualified=max(density_errors)<=1e-12 and srherr<=1e-7))
        for name in names[:-1]:
            nr=[native[name][str(i)] for i in ids];vr=values[name]
            diff=[math.log10(v/n) if 'Density' in name else v-n for v,n in zip(vr,nr)]
            peak=max(range(len(ids)),key=lambda k:abs(diff[k]));volume=math.fsum(weights)
            fields.append(dict(**label,field=name,units='dex' if 'Density' in name else 'V',mean_signed=math.fsum(w*x for w,x in zip(weights,diff))/volume,weighted_rms=math.sqrt(math.fsum(w*x*x for w,x in zip(weights,diff))/volume),max_abs=abs(diff[peak]),max_node=ids[peak]))
        nr=[native['srhRecombination'][str(i)] for i in ids];vr=values['srhRecombination']
        base=math.fsum(w*abs(v) for w,v in zip(weights,nr))
        sources.append(dict(**label,native_charge_equivalent_A_per_um=1.602176634e-19*math.fsum(w*v for w,v in zip(weights,nr)),vela_charge_equivalent_A_per_um=1.602176634e-19*math.fsum(w*v for w,v in zip(weights,vr)),normalized_L1=math.fsum(w*abs(v-n) for w,v,n in zip(weights,vr,nr))/base if base else None))
    for name,data in (('fields',fields),('srh',sources),('reconstruction',checks)):csvout(run.summary/(name+'.csv'),data)
    write(run.summary/'field_summary.json',dict(points=len(checks),qualified=sum(r['qualified'] for r in checks),all_reconstruction_qualified=all(r['qualified'] for r in checks),density_units='cm^-3; state CSV m^-3',SRH_units='cm^-3 s^-1',SRH_meaning='Common positive barycentric Si integration; not terminal current or own source volume. Production source reconstruction uses unchanged all_cell volume.'))
    assert len(checks)==204 and all(r['qualified'] for r in checks),'Field reconstruction failed'


def main():
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('migrate','seal','controls','curves','summarize','all'))
    parser.add_argument('--root',type=Path,required=True);parser.add_argument('--runner',type=Path,required=True);parser.add_argument('--workers',type=int,default=4)
    args=parser.parse_args();run=Run(args.root,args.runner,args.workers)
    if args.action=='migrate':run.migrate();return
    run.seal()
    if args.action=='seal':return
    if args.action=='summarize':run.summarize();return
    try:
        if args.action in ('controls','all'):run.controls()
        if args.action in ('curves','all'):run.curves()
    finally:run.summarize()


if __name__=='__main__':main()
