"""Independently check the bounded runtime/SG/source calibration and seal evidence."""
import argparse
import ast
import copy
import math
from pathlib import Path
import re
import subprocess

import numpy as np
import validate_simplemos_distributed_transport_20260907 as v

p=v.p; a=p.a; d=p.d; OUT=p.OUT; LOCAL=p.LOCAL
REPORT=p.REPO/'docs/validation/simplemos_native_mobility_flux_calibration_2026-09-07.md'
ARMS={'combined':(v.OUT,v.LOCAL),
      'geometry':(OUT/'factor_resolved_response/geometry',LOCAL/'factor_resolved_response/geometry'),
      'mobility':(OUT/'factor_resolved_response/mobility',LOCAL/'factor_resolved_response/mobility')}


def same(x,y):
    assert math.isclose(x,y,rel_tol=1e-11,abs_tol=1e-26),(x,y)


def record_csv(path,rows):
    if path.exists():assert a.rows(path)==[{k:str(x) for k,x in row.items()} for row in rows]
    else:a.write_csv(path,rows)


def manifests():
    return [p.prior.OUT/'validation_evidence.json']+[OUT/n for n in (
        'native_freeze.json','box_mobility_scope.json','response_analysis_scope.json',
        'supported_runtime/freeze.json','geometry_runtime/freeze.json','geometry_runtime/evidence.json',
        'distributed_response/freeze.json','factor_response/geometry/freeze.json','factor_response/mobility/freeze.json',
        'factor_resolved_response/geometry/freeze.json','factor_resolved_response/mobility/freeze.json')]


def physics(config):
    config=copy.deepcopy(config)
    config.pop('state_file');config.pop('output_state_file')
    config['solver']['local_update_diagnostics'].pop('csv_file')
    return config


def row_source(rows,count):
    rhs=np.zeros(count);absolute=np.zeros(count)
    for row in rows:
        i,j=int(row['node0']),int(row['node1']);f=float(row['source_normalized'])
        rhs[i]+=f;rhs[j]-=f;absolute[i]+=abs(f);absolute[j]+=abs(f)
    return rhs,absolute


def check():
    for file in manifests():a.verify(file)
    checks=[];fd=[];preflight=[];summary=[];decomposition=[];original_failures=[]
    contracts={arm:a.read(paths[0]/'contract.json') for arm,paths in ARMS.items()}
    sources={arm:a.rows(paths[0]/'source_edges.csv') for arm,paths in ARMS.items()}
    geometry_cache={device:v.mobility.geom.geometry(device) for device in ('n19','n23')}
    native_weights={device:v.mobility.geometry(device)[1] for device in geometry_cache}
    for arm,(out,local) in ARMS.items():
        dc=a.rows(out/'dc.csv');cal=a.rows(out/'calibration.csv');assert len(dc)==20 and len(cal)==8
        assert contracts[arm]['frozen_gates']==a.read(OUT/'native_contract.json')['eventual_response_gates']
        for c in contracts[arm]['cases']:
            geo,xy,cells,vol,K,parts,info,modified=geometry_cache[c['device']]
            mask=d.matrix.spatial.old.m78.supports(c['device'],geo,.05)[0]['all_si'].copy();mask[geo.contact_nodes]=False
            assert 2*int(sum(mask))==1814
            base=Path(c['base']);baseline=a.read(base/'config.json');mesh=a.read(Path(baseline['mesh_file']))
            drain=set(next(x['node_ids'] for x in mesh['contacts'] if x['name']=='drain'))
            rows=[x for x in sources[arm] if x['key']==c['key']]
            edges={int(x['edge_id']):x for x in a.rows(p.prior.LOCAL/'vela'/c['key']/'strict/edges.csv')}
            file_source={int(x.split()[0]):float(x.split()[1]) for x in Path(c['source']).read_text().splitlines()}
            assert len(rows)==len(edges)==len(file_source)
            for row in rows:
                edge=edges[int(row['edge'])];i,j=int(row['node0']),int(row['node1'])
                assert (i,j)==(int(edge['node0']),int(edge['node1']))
                same(float(row['source_injection_native_units']),file_source[int(row['edge'])])
                same(float(row['source_particle_per_m_s']),.01*file_source[int(row['edge'])])
                g0=float(edge['couple_m'])/float(edge['length_m']);mu0=float(edge['electron_mobility_m2_V_s'])
                gn=math.fsum(x['coefficient'] for x in parts[(i,j)] if x['material']=='Si')
                if g0<=0 or mu0<=0:
                    assert gn==0 and native_weights[c['device']][(i,j)]['e']==0;gr=fr=1.
                else:
                    gr=gn/g0;fr=native_weights[c['device']][(i,j)]['e']*1e-4/(g0*mu0)
                    if abs(gr-1)<=1e-10:gr=1.
                delta={'combined':fr-1,'geometry':gr-1,'mobility':fr-gr}[arm]
                same(float(row['source_normalized']),float(edge['electron_flux'])*delta)
                same(float(row['source_particle_per_m_s']),float(edge['electron_particle_line_flux_per_m_s'])*delta)
            rhs,absolute=row_source(rows,geo.count);rhs[geo.contact_nodes]=0.
            adj=d.ordered(p.prior.LOCAL/'vela'/c['key']/'adjoint.csv',geo.count)
            feedback=-math.fsum(float(row['lambda_electron'])*float(rhs[i]) for i,row in enumerate(adj))
            direct=-d.fixed.Q*1e-6*math.fsum((int(int(row['node0']) in drain)-int(int(row['node1']) in drain))*float(row['source_particle_per_m_s']) for row in rows)
            same(feedback,c['unit_feedback_A_per_um']);same(direct,c['unit_direct_A_per_um'])
            prediction=feedback+direct;same(prediction,c['unit_prediction_A_per_um'])
            root=local/c['key'];z=d.ordered(root/'preflight/zero/residual.csv',geo.count)
            unit=d.ordered(root/'preflight/unit/residual.csv',geo.count)
            old=d.ordered(p.prior.LOCAL/'vela'/c['key']/'strict/residual.csv',geo.count)
            pre_ratios=[]
            for i,(x,y,o) in enumerate(zip(z,unit,old)):
                for k in ('psi_residual','phin_residual','phip_residual'):assert x[k]==o[k]
                for k in ('psi_residual','phip_residual'):assert x[k]==y[k]
                error=abs(float(y['phin_residual'])-float(x['phin_residual'])-rhs[i])
                bound=max(1e-8*absolute[i],1e-25);assert error<=bound;pre_ratios.append(error/bound)
            zz=a.read(root/'preflight/zero/config.status.json');uu=a.read(root/'preflight/unit/config.status.json')
            assert zz['exit_code']==uu['exit_code']==0
            assert abs(uu['current_A_per_um']-zz['current_A_per_um']-direct)<=max(1e-8*abs(direct),1e-14*abs(c['base_Id_A_per_um']))
            preflight.append(dict(arm=arm,key=c['key'],zero_residual_bit_identical=True,max_error_over_frozen_bound=max(pre_ratios),qualified=True))
            currents={}
            for job in c['jobs']:
                path=Path(job['config']);dest=path.parent;cfg=a.read(path);assert physics(cfg)==physics(baseline)
                assert cfg['state_file']==str(base/'state.csv')
                post=a.read(dest/'all_row.json');assert post['solver']['carrier_row_convergence']['mode']=='report'
                assert post['solver']['global_continuity_closure']==dict(mode='enforce',tolerance=1e-6,source_floor=1e-10)
                status=a.read(dest/'config.status.json');q=a.read(dest/'all_row.status.json')
                ratios=[]
                for row,keep in zip(d.ordered(dest/'all_row.csv',geo.count),mask):
                    if not keep:continue
                    for carrier in ('electron','hole'):
                        scale=max(float(row[carrier+'_flux_abs_sum']),abs(float(row[carrier+'_recombination'])),abs(float(row[carrier+'_impact'])))
                        assert scale>0;ratios.append(abs(float(row[carrier+'_residual']))/scale)
                g=q['carrier_row_convergence'];assert len(ratios)==g['qualified_row_count']==1814
                bad=sum(x>1e-6 for x in ratios);assert bad==g['violation_count']==0;same(max(ratios),g['max_ratio'])
                cc=status['contact_currents_A_per_um'];current=cc['drain'];kcl=abs(math.fsum(cc.values()))/abs(current)
                passed=status['exit_code']==q['exit_code']==0 and status['converged'] and q['global_continuity_closure']['satisfied'] and kcl<=1e-8
                reported=next(x for x in dc if x['key']==c['key'] and x['label']==job['label'])
                assert passed and reported['qualified']=='True';same(current,float(reported['current_A_per_um']))
                linear=a.rows(dest/'linear_summary.csv')
                for iteration in {x['iteration'] for x in linear}:assert [int(x['correction']) for x in linear if x['iteration']==iteration]==list(range(5))
                currents[job['label']]=current
                checks.append(dict(arm=arm,key=c['key'],label=job['label'],active_rows=len(ratios),violations=bad,zero_scale_rows=0,
                    max_row_ratio=max(ratios),kcl_over_Id=kcl,frozen_physics_and_solver_unchanged=True,qualified=True))
            center=currents['zero'];odd={n:(currents['plus_'+n]-currents['minus_'+n])/2 for n in ('full','half')}
            linearity=abs(odd['full']/(2*odd['half'])-1);drift=abs(center-c['base_Id_A_per_um'])
            for name,amplitude in [('full',.001),('half',.0005)]:
                error=abs(odd[name]/(amplitude*prediction)-1);even=abs(((currents['plus_'+name]+currents['minus_'+name])/2-center)/odd[name])
                snr=abs(odd[name])/max(drift,1e-300);dd=abs(math.log10(center/c['base_Id_A_per_um']))
                signs=(currents['plus_'+name]-center)*prediction>0 and (currents['minus_'+name]-center)*prediction<0
                reported=next(x for x in cal if x['key']==c['key'] and x['amplitude']==name)
                assert abs(error-float(reported['prediction_relative_error']))<=8*np.finfo(float).eps
                assert abs(linearity-float(reported['two_amplitude_relative']))<=8*np.finfo(float).eps
                assert error<=.001 and linearity<=.001 and even<=.01 and snr>=100 and dd<=1e-5 and signs and reported['qualified']=='True'
                fd.append(dict(arm=arm,key=c['key'],amplitude=name,prediction_relative_error=error,two_amplitude_relative=linearity,
                    even_over_odd=even,signal_over_zero_drift=snr,zero_drift_dex=dd,signs_correct=signs,qualified=True))
            summary.append(dict(arm=arm,key=c['key'],base_Id_A_per_um=c['base_Id_A_per_um'],unit_dId_dalpha_A_per_um=prediction,
                unit_relative_Id_derivative=prediction/c['base_Id_A_per_um'],direct_fraction=direct/prediction))
            print('Independently checked',arm,c['key'],'5 states / 2 amplitudes',flush=True)

    for c in contracts['combined']['cases']:
        geo=geometry_cache[c['device']][0]
        keyed={arm:{int(x['edge']):x for x in sources[arm] if x['key']==c['key']} for arm in ARMS}
        closure=[]
        for edge,row in keyed['combined'].items():
            g=keyed['geometry'][edge];m=keyed['mobility'][edge]
            x=float(row['source_normalized']);y=float(g['source_normalized'])+float(m['source_normalized'])
            err=abs(x-y)/max(abs(x),abs(float(g['source_normalized'])),abs(float(m['source_normalized'])),1e-300)
            assert err<=1e-12;closure.append(err)
        original=a.read(OUT/'factor_response/geometry/contract.json')
        oc=next(x for x in original['cases'] if x['key']==c['key']);oldrows=[x for x in a.rows(OUT/'factor_response/geometry/source_edges.csv') if x['key']==c['key']]
        root=LOCAL/'factor_response/geometry'/c['key'];z=d.ordered(root/'preflight/zero/residual.csv',geo.count);u=d.ordered(root/'preflight/unit/residual.csv',geo.count)
        rhs,absolute=row_source(oldrows,geo.count);rhs[geo.contact_nodes]=0.
        bad=sum(abs(float(y['phin_residual'])-float(x['phin_residual'])-rhs[i])>max(1e-8*absolute[i],1e-25) for i,(x,y) in enumerate(zip(z,u)))
        assert bad>0 and not list(root.glob('*/state.csv'))
        original_failures.append(dict(key=c['key'],failed_preflight_rows=bad,self_consistent_DC_run=False,original_frozen_gates_retained=True))
        pred={arm:next(x['unit_dId_dalpha_A_per_um'] for x in summary if x['arm']==arm and x['key']==c['key']) for arm in ARMS}
        same(pred['geometry']+pred['mobility'],pred['combined'])
        resolved=keyed['geometry'];oldkeyed={int(x['edge']):x for x in oldrows}
        adj=d.ordered(p.prior.LOCAL/'vela'/c['key']/'adjoint.csv',geo.count);lam=np.array([float(x['lambda_electron']) for x in adj]);lam[geo.contact_nodes]=0
        mesh=a.read(Path(a.read(Path(c['base'])/'config.json')['mesh_file']));drain=set(next(x['node_ids'] for x in mesh['contacts'] if x['name']=='drain'))
        bound=[];count=0
        for edge,g in resolved.items():
            o=oldkeyed[edge];df=float(g['source_normalized'])-float(o['source_normalized']);dp=float(g['source_particle_per_m_s'])-float(o['source_particle_per_m_s'])
            i,j=int(g['node0']),int(g['node1']);count+=df!=0 or dp!=0
            bound.append(abs((lam[i]-lam[j])*df)+abs(d.fixed.Q*1e-6*(int(i in drain)-int(j in drain))*dp))
        decomposition.append(dict(key=c['key'],source_sum_max_relative_roundoff=max(closure),
            unit_prediction_sum_relative=abs((pred['geometry']+pred['mobility'])/pred['combined']-1),
            resolved_geometry_changed_edges=count,resolution_unit_prediction_change_over_Id=(pred['geometry']-oc['unit_prediction_A_per_um'])/c['base_Id_A_per_um'],
            resolution_sum_absolute_edge_response_bound_over_Id=math.fsum(bound)/abs(c['base_Id_A_per_um'])))

    # Recheck native return codes, returned input bytes, and raw port closure.
    native=[]
    for version,local in [('unsupported_edge',LOCAL),('unsupported_potential',LOCAL/'supported'),('geometry',LOCAL/'geometry_runtime')]:
        for f in (local/'bundle').rglob('*'):
            if f.is_file():assert a.sha(f)==a.sha(local/'native_raw/bundle'/f.relative_to(local/'bundle'))
        for c in a.read(OUT/'native_contract.json')['jobs']:
            root=local/'native_raw/bundle'/c['key'];code=int((root/'exit_code.txt').read_text());log=(root/'console.log').read_text(errors='replace')
            if version=='geometry':
                assert code==0 and 'Good Bye' in log and 'T-2022.03-SP2' in log
                endpoint=p.prior.prev.p.exporter.pltrows(root/'native_des.plt')[-1]
                cc={name:endpoint[name+' TotalCurrent'] for name in ('drain','source','gate','substrate')}
                drift=abs(cc['drain']/c['native_Id_A_per_um']-1);kcl=abs(math.fsum(cc.values()))/abs(cc['drain']);assert drift<=1e-8 and kcl<=1e-8
            else:
                assert code==5 and ('eMobility' if version=='unsupported_edge' else 'Potential') in log
                drift=kcl='unqualified_runtime_export'
            native.append(dict(version=version,key=c['key'],exit_code=code,Id_relative_drift=drift,kcl_over_Id=kcl))
    runtime=a.rows(OUT/'geometry_runtime/checks.csv');assert len(runtime)==4 and all(x['qualified']=='True' for x in runtime)
    identity=a.rows(OUT/'box_mobility_identity.csv');assert sum(int(x['cells']) for x in identity)==7000 and all(x['qualified']=='True' for x in identity)
    ports=a.rows(OUT/'native_sg_port_checks.csv');assert len(ports)==8
    sg=a.rows(OUT/'native_sg_edges.csv')
    for row in ports:
        c=next(x for x in contracts['combined']['cases'] if x['key']==row['key']);mesh=a.read(Path(a.read(Path(c['base'])/'config.json')['mesh_file']))
        edges=[e for e in sg if e['key']==c['key'] and e['mode']==row['mode']];current={}
        for contact in mesh['contacts']:
            nodes=set(contact['node_ids']);current[contact['name']]=-v.mobility.Q*1e-4*math.fsum((int(int(e['node0']) in nodes)-int(int(e['node1']) in nodes))*(float(e['electron_particle_flux_per_cm_s'])-float(e['hole_particle_flux_per_cm_s'])) for e in edges)
        same(current['drain'],float(row['Id_A_per_um']))
        assert abs(current['drain']/float(row['target_A_per_um'])-1)<=1e-6 and abs(math.fsum(current.values()))/abs(current['drain'])<=1e-8
    for name,rows in [('independent_qualification',checks),('independent_calibration',fd),('independent_preflight',preflight),
                      ('factor_summary',summary),('factor_additivity_and_resolution',decomposition),('retained_preflight_failures',original_failures),('native_attempts',native)]:
        record_csv(OUT/(name+'.csv'),rows)
    metrics=dict(qualified_states=len(checks),independently_recomputed_carrier_rows=len(checks)*1814,qualified_response_amplitudes=len(fd),
        max_row_ratio=max(x['max_row_ratio'] for x in checks),max_kcl_over_Id=max(x['kcl_over_Id'] for x in checks),
        max_prediction_relative_error=max(x['prediction_relative_error'] for x in fd),max_two_amplitude_relative=max(x['two_amplitude_relative'] for x in fd),
        max_even_over_odd=max(x['even_over_odd'] for x in fd),minimum_signal_over_zero_drift=min(x['signal_over_zero_drift'] for x in fd),
        max_geometry_resolution_response_bound_over_Id=max(x['resolution_sum_absolute_edge_response_bound_over_Id'] for x in decomposition),
        native_successful_DC=4,native_failed_runtime_attempts=8,original_failed_preflight_pairs=4)
    a.write(OUT/'independent_check.json',dict(status='passed',metrics=metrics,
        input_hashes={a.rel(f):a.sha(f) for f in sorted(set([Path(__file__).resolve()]+manifests()+[x for root in (LOCAL,OUT) for x in root.rglob('*') if x.is_file()]))}))
    print(metrics,flush=True)


def seal():
    a.verify(OUT/'independent_check.json')
    for file in manifests():a.verify(file)
    scripts=[p.REPO/'scripts'/name for name in (
        'prepare_simplemos_masetti_runtime_20260907.py','prepare_simplemos_runtime_supported_20260907.py',
        'prepare_simplemos_runtime_geometry_only_20260907.py','audit_simplemos_runtime_geometry_20260907.py',
        'audit_simplemos_masetti_box_mobility_20260907.py','build_simplemos_distributed_flux_20260907.py',
        'validate_simplemos_distributed_transport_20260907.py','analyze_simplemos_transport_calibration_20260907.py',
        'validate_simplemos_transport_factors_20260907.py','validate_simplemos_transport_factors_resolved_20260907.py',
        'seal_simplemos_transport_calibration_20260907.py')]
    scripts.append(p.REPO/'tests/regression/test_simplemos_native_sg_calibration.py')
    for file in scripts:ast.parse(file.read_text(encoding='utf-8-sig'),filename=str(file))
    files=set(scripts+manifests()+[REPORT,OUT/'independent_check.json',p.REPO/'src/tools/vela_example_runner.cpp',
        p.REPO/'src/equation/CoupledDDAssembler.cpp',p.REPO/'src/physics/MobilityModel.cpp'])
    for root in (LOCAL,OUT):files.update(x for x in root.rglob('*') if x.is_file() and x.name!='validation_evidence.json')
    for target in re.findall(r'\]\(([^)]+)\)',REPORT.read_text(encoding='utf-8')):
        if '://' in target:continue
        file=(REPORT.parent/target.split('#')[0]).resolve()
        if file.name!='validation_evidence.json':assert file.exists(),file;files.add(file)
    diff=subprocess.run(['git','-c','core.fsmonitor=false','diff','--numstat','--ignore-space-at-eol'],cwd=p.REPO,capture_output=True,text=True)
    assert diff.returncode==0 and diff.stdout.strip()=='127\t1\tsrc/tools/vela_example_runner.cpp',diff.stdout
    a.write(OUT/'validation_evidence.json',dict(schema_version=1,completed_date='2026-09-07',status='completed_bounded_native_mobility_flux_calibration',
        metrics=a.read(OUT/'independent_check.json')['metrics'],input_hashes={a.rel(f):a.sha(f) for f in sorted(files)},
        verified_manifests=len(manifests()),python_syntax_checked=len(scripts),
        observed_unit_tests=dict(command='D:/msys64/ucrt64/bin/python.exe tests/regression/test_simplemos_native_sg_calibration.py',tests=5,result='passed',source='Executed in this task; output recorded in conversation.'),
        production_changes=False,acceptance_changes=False,m82_released=False,m83_released=False,
        limitations=['Four Vg=1 V endpoints; no new full 0-1 V field sweep.',
            'Electron distributed sources only; no hole source or Poisson source response calibration.',
            'The alpha=1 coefficient replacement was not simulated; directional derivatives are not final current improvements.',
            'Conservative SG flux was calibrated through terminal and continuity replay; plotted current vectors are not claimed identical to edge flux.',
            'Eight failed native runtime exports and four failed original geometry preflight pairs remain retained.',
            'Near-identity geometry ratios were resolved to one within 1e-10 only in a newly frozen diagnostic source; no solver gate was changed.']))
    a.verify(OUT/'validation_evidence.json');print('Sealed',len(files),'files',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('check','seal'));globals()[parser.parse_args().action]()
