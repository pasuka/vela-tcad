"""Opt-in production cell-material Poisson migration, immutable eight-point audit.

Native mobility remains an isolated, separately calibrated diagnostic input.
No historical evidence, acceptance threshold or global default is rewritten.
"""
import argparse, copy, math, subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import validate_simplemos_remaining_dielectric_20260907 as old
import build_simplemos_native_mobility_candidate_20260907 as mob
import build_simplemos_joint_geometry_20260907 as jbuild

a=old.a;d=old.d;v=old.v;p=old.p
LOCAL=p.REPO/'build-release/simplemos_permittivity_production_20260907'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/permittivity_production_20260907'
OVERLAY=LOCAL/'mobility_diagnostic';RUNNER=OVERLAY/'runner.exe'

def cases():return a.read(old.OUT/'contract.json')['cases']

def build():
    OVERLAY.mkdir(exist_ok=False)
    header=(p.REPO/'include/vela/equation/AssemblerUtils.h').read_text()
    start=header.index('inline Real edgeMobility(');end=header.index('/// Return average model mobility',start)
    body=header[start:end];marker='    return sum / static_cast<Real>(contributingCells);';assert body.count(marker)==1
    body=body.replace(marker,'    return (sum / static_cast<Real>(contributingCells)) * vela_candidate::factor(edgeId, mesh.edges().size(), carrier == CarrierType::Electron);')
    target=OVERLAY/'include/vela/equation/AssemblerUtils.h';target.parent.mkdir(parents=True)
    target.write_text('#include <fstream>\n#include <cstdlib>\n#include <array>\n#include <vector>\n#include <cmath>\n#include <stdexcept>\n'+mob.HOOK+header[:start]+body+header[end:])
    source=(p.REPO/'src/equation/CoupledDDAssembler.cpp').read_text()
    marker='    return sum / static_cast<Real>(lowFieldMobilities.size());';assert source.count(marker)==1
    source=source.replace(marker,'    return (sum / static_cast<Real>(lowFieldMobilities.size())) * vela_candidate::factor(edgeId, mesh_.numEdges(), carrier == CarrierType::Electron);')
    (OVERLAY/'CoupledDDAssembler.cpp').write_text(source)
    runner=(p.REPO/'src/tools/vela_example_runner.cpp').read_text();marker='int main(int argc, char** argv)';assert runner.count(marker)==1
    runner=runner.replace(marker,jbuild.FUNCTION+'\n'+marker)
    marker='        } else if (type == "newton_residual_probe") {';assert runner.count(marker)==1
    runner=runner.replace(marker,'        } else if (type == "joint_geometry_jvp") {\n            status.update(runJointGeometryJvp(configFile,cfg));\n'+marker)
    (OVERLAY/'runner.cpp').write_text(runner)
    jobs=[]
    for cmd in a.read(old.q.oldLOCAL/'compile_commands.json'):
        cmd=list(cmd);stem=Path(cmd[cmd.index('-o')+1]).stem
        cmd[1]='-I'+str(OVERLAY/'include');cmd[cmd.index('-o')+1]=str(OVERLAY/(stem+'.o'))
        if stem=='CoupledDDAssembler':cmd[cmd.index('-c')+1]=str(OVERLAY/'CoupledDDAssembler.cpp')
        if stem=='vela_example_runner':cmd[cmd.index('-c')+1]=str(OVERLAY/'runner.cpp')
        jobs.append(cmd)
    a.write(OVERLAY/'compile_commands.json',jobs)
    def compile(cmd):
        r=subprocess.run(cmd,cwd=p.REPO/'build-release',env=v.environment(),capture_output=True,text=True)
        (OVERLAY/(Path(cmd[cmd.index('-o')+1]).stem+'.build.log')).write_text(r.stdout+r.stderr)
        assert r.returncode==0,r.stderr[-2500:]
        print('Diagnostic compiled',Path(cmd[cmd.index('-o')+1]).name,flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(compile,jobs))
    prior=a.read(old.q.LOCAL/'link_command.json');libs=[x for x in prior if x.endswith('.a')]
    cmd=['D:/msys64/ucrt64/bin/g++.exe']+[x[x.index('-o')+1] for x in jobs]+libs+['-o',str(RUNNER)]
    a.write(OVERLAY/'link_command.json',cmd)
    r=subprocess.run(cmd,env=v.environment(),capture_output=True,text=True);(OVERLAY/'link.log').write_text(r.stdout+r.stderr)
    assert r.returncode==0,r.stderr[-2500:]

def prepare():
    # Previous artifacts are immutable; only the source/binary preimages supersede.
    a.verify(old.OUT/'validation_evidence.json');a.verify(old.l.OUT/'final_evidence.json')
    entries={};geometry=[];files=[Path(__file__).resolve(),OUT/'prechange.json',RUNNER,v.RUNNER,old.OUT/'validation_evidence.json',old.l.OUT/'final_evidence.json']
    for device in ('n19','n23'):
        geo,xy,elements,vol,K,parts,*_=old.l.g.geometry(device)
        case=next(c for c in cases() if c['device']==device)
        mesh=a.read(Path(v.config(case)['mesh_file']));data=[];error=0.
        for cell in mesh['triangles']:
            nodes=cell['node_ids'];cid=cell['id'];assert set(nodes)==set(elements[cid]['nodes'])
            values=[]
            for k in range(3):
                i,j=sorted((nodes[k],nodes[(k+1)%3]));part=next(x for x in parts[(i,j)] if x['cell']==cid)
                values.append(part['coefficient'])
            data.append(dict(cell_id=cid,node_ids=nodes,coefficients=values))
        assert len(data)==len(elements)
        entries[device]=data;path=LOCAL/(device+'_cell_coefficients.json');a.write(path,data);files.append(path)
        geometry.append(dict(device=device,cells=len(data),coefficients=3*len(data),minimum=min(x for r in data for x in r['coefficients']),maximum=max(x for r in data for x in r['coefficients']),node_identity_verified=True))
    allcases=cases()
    for c in allcases:
        c['migration_jobs']=[];c['migration_probes']=[];cfg=v.config(c)
        for name in ('mesh_file','materials_file','node_doping_file'):files.append(Path(cfg[name]))
        initial=v.LOCAL/c['key']/'from_vela/state.csv'
        native_initial=v.m.previous.LOCAL/c['key']/'native_referenced_joint/state.csv'
        for arm,mu,seed in [('legacy',False,initial),('cell_local',False,initial),
                             ('native_cells',False,initial),('native_cells_native_seed',False,native_initial),
                             ('native_cells_native_mu',True,old.q.LOCAL/c['key']/'replacement/state.csv'),
                             ('native_cells_native_mu_native_seed',True,native_initial)]:
            deck=copy.deepcopy(cfg);geometry=deck.setdefault('mesh_geometry',{})
            if arm!='legacy':geometry['poisson_permittivity_policy']='cell_material'
            if arm.startswith('native_cells'):geometry['poisson_cell_edge_coefficients']=entries[c['device']]
            dest=LOCAL/c['key']/arm;deck.update(state_file=str(seed),output_state_file=str(dest/'state.csv'))
            a.write(dest/'config.json',deck);v.post_config(deck,dest);probes=v.probes(deck,dest/'post',dest/'state.csv')
            files += [seed,dest/'config.json',dest/'all_row.json',dest/'acceptance_edges.json']+probes
            c['migration_jobs'].append(dict(arm=arm,mu=mu,dest=str(dest)))
        for label in ('zero','all_native_K'):
            deck=copy.deepcopy(cfg)
            if label=='all_native_K':deck.setdefault('mesh_geometry',{}).update(poisson_permittivity_policy='cell_material',poisson_cell_edge_coefficients=entries[c['device']])
            dest=LOCAL/c['key']/'fixed'/label
            paths=v.probes(deck,dest,old.q.LOCAL/c['key']/'replacement/state.csv')
            jvp=copy.deepcopy(a.read(dest/'edges.json'));jvp.update(simulation_type='joint_geometry_jvp',output_csv=str(dest/'jvp.csv'));a.write(dest/'jvp.json',jvp)
            files+=paths+[dest/'jvp.json'];c['migration_probes'] += [str(dest/(n+'.json')) for n in ('functional','edges','jvp')]
    a.write(OUT/'geometry_mapping.json',dict(devices=[dict(device=k,cells=len(x),coefficient_count=3*len(x),node_order_verified=True) for k,x in entries.items()]))
    a.write(OUT/'contract.json',dict(cases=allcases,DC=48,gates=a.read(old.OUT/'contract.json')['gates'],
        purpose='Production cell-material Poisson assembly, raw-cell control and native per-cell geometry; dual initialization on native geometry with production and diagnostic native mobility.',
        fixed='Prior signed Si Poisson charge volumes, native Si transport support, original SRH volume, Vela constants, physics and all acceptance gates.',
        controls='8 legacy replays; 8 raw local-cell runs; 16 native-cell runs with production mobility; 16 native-cell runs with calibrated isolated native mobility.',
        default_changed=False,full_curve=False))
    files += [OUT/'contract.json',OUT/'geometry_mapping.json',p.REPO/'build-release/libvela_core.a']
    files += [x for base in ('src','include') for x in (p.REPO/base).rglob('*') if x.suffix in ('.cpp','.h')]
    files += [x for x in OVERLAY.rglob('*') if x.is_file()]
    d.matrix.freeze(OUT/'freeze.json',files);print('Frozen 48 DC, original eight points and gates',flush=True)

def execute(path,c,mu):return v.execute(path,RUNNER if mu else v.RUNNER,old.l.env(c,1) if mu else v.environment())

def preflight():
    a.verify(OUT/'freeze.json');checks=[];jvps=[]
    for c in a.read(OUT/'contract.json')['cases']:
        for path in c['migration_probes']:
            result=execute(Path(path),c,True);assert result['exit_code']==0,(path,result)
        geo,*_=old.l.g.geometry(c['device'])
        for label in ('zero','all_native_K'):
            new=LOCAL/c['key']/'fixed'/label;prior=old.LOCAL/c['key']/'fixed'/label
            nr=d.array(d.ordered(new/'residual.csv',geo.count),('psi_residual','phin_residual','phip_residual'))
            pr=d.array(d.ordered(prior/'residual.csv',geo.count),('psi_residual','phin_residual','phip_residual'))
            # Use the actual dielectric source as scale, not a converged residual.
            # The legacy replay is required to remain bit-identical.
            base=d.array(d.ordered(old.LOCAL/c['key']/'fixed/zero/residual.csv',geo.count),('psi_residual',))[0]
            err=(0. if np.array_equal(nr[0],pr[0]) else math.inf) if label=='zero' else float(np.linalg.norm((nr-pr)[0,geo.free])/max(np.linalg.norm((pr[0]-base)[geo.free]),1e-300))
            carrier=bool(np.array_equal(nr[1:],pr[1:]));ne=a.rows(new/'edges.csv');pe=a.rows(prior/'edges.csv')
            transport=all(x[k]==y[k] for x,y in zip(ne,pe) for k in ('electron_flux','hole_flux'))
            s=a.read(new/'functional.status.json');ps=a.read(prior/'functional.status.json')
            current=abs(s['current_A_per_um']/ps['current_A_per_um']-1)
            _,base_jvp=old.jc.jvp_metrics(prior/'jvp.csv',c['key'],label)
            jp,_=old.jc.jvp_metrics(new/'jvp.csv',c['key'],label,base_jvp);jvps+=jp
            checks.append(dict(key=c['key'],arm=label,poisson_relative=err,carrier_bitwise=carrier,transport_bitwise=transport,current_relative=current,
                               qualified=err<=1e-8 and carrier and transport and current<=1e-8 and all(r['qualified'] for r in jp)))
        print('Preflight',c['key'],checks[-1]['qualified'],flush=True)
    a.write_csv(OUT/'preflight.csv',checks);a.write_csv(OUT/'jvp.csv',jvps)
    d.matrix.freeze(OUT/'preflight_evidence.json',[OUT/'freeze.json',OUT/'preflight.csv',OUT/'jvp.csv']+[x for c in cases() for x in (LOCAL/c['key']/'fixed').rglob('*') if x.is_file()])
    assert all(x['qualified'] for x in checks)

def run():
    a.verify(OUT/'freeze.json');a.verify(OUT/'preflight_evidence.json');assert all(r['qualified']=='True' for r in a.rows(OUT/'preflight.csv'))
    def one(c):
        rows=[]
        for j in c['migration_jobs']:
            dest=Path(j['dest']);s=execute(dest/'config.json',c,j['mu'])
            row=dict(key=c['key'],device=c['device'],vg=c['vg'],vd=c['vd'],arm=j['arm'],qualified=False)
            if (dest/'state.csv').exists():
                for name in ('all_row','acceptance_edges'):execute(dest/(name+'.json'),c,j['mu'])
                row.update(old.qualify(c,dest))
                for name in ('functional','edges','terms'):execute(dest/'post'/(name+'.json'),c,j['mu'])
            else:row['failure']=s.get('failure_reason','missing state')
            rows.append(row);print(c['key'],j['arm'],row['qualified'],row.get('max_row_ratio'),flush=True)
        return rows
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(one,a.read(OUT/'contract.json')['cases']) for r in group]
    keys=list(dict.fromkeys(k for r in rows for k in r));a.write_csv(OUT/'dc.csv',[{k:r.get(k,'') for k in keys} for r in rows])

def analyze():
    a.verify(OUT/'freeze.json');rows=a.rows(OUT/'dc.csv');comparison=[];duals=[];ports=[];replay=[]
    for c in cases():
        selected={r['arm']:r for r in rows if r['key']==c['key']};geo,mask=v.m.previous.prior.support(c)
        native=float(c['native_Id_A_per_um']);historical=next(r for r in a.rows(old.OUT/'comparison.csv') if r['key']==c['key'] and r['axis']=='all_native_K')
        for arm in selected:
            r=selected[arm];Id=float(r['current_A_per_um']);dest=LOCAL/c['key']/arm
            comparison.append(dict(key=c['key'],device=c['device'],vg=c['vg'],vd=c['vd'],arm=arm,current_A_per_um=Id,native_current_A_per_um=native,relative_error=Id/native-1,qualified=r['qualified']=='True'))
            f=a.read(dest/'post/functional.status.json');port=abs(f['current_A_per_um']/f['contact_current_extractor_A_per_um']-1)
            ports.append(dict(key=c['key'],arm=arm,relative_difference=port,qualified=port<=1e-8))
        for arm in ('native_cells','native_cells_native_mu'):
            s0=d.ordered(LOCAL/c['key']/arm/'state.csv',geo.count);s1=d.ordered(LOCAL/c['key']/(arm+'_native_seed')/'state.csv',geo.count)
            values=v.m.previous.prior.delta_states(s0,s1,mask);difference=abs(float(selected[arm]['current_A_per_um'])/float(selected[arm+'_native_seed']['current_A_per_um'])-1)
            passed=max(values[k] for k in ('psi_max_V','phin_max_V','phip_max_V'))<=1e-6 and values['density_max_relative']<=1e-4 and difference<=1e-6
            duals.append(dict(key=c['key'],arm=arm,**values,Id_relative=difference,qualified=passed))
        legacy=float(selected['legacy']['current_A_per_um']);base=a.read(v.LOCAL/c['key']/'from_vela/config.status.json')['contact_currents_A_per_um']['drain']
        drift0=abs(legacy/base-1);drift1=abs(float(selected['native_cells_native_mu']['current_A_per_um'])/float(historical['candidate_Id_A_per_um'])-1)
        replay.append(dict(key=c['key'],legacy_current_drift=drift0,native_cell_migration_current_drift=drift1,qualified=max(drift0,drift1)<=1e-8))
    a.write_csv(OUT/'comparison.csv',comparison);a.write_csv(OUT/'dual_initialization.csv',duals);a.write_csv(OUT/'ports.csv',ports);a.write_csv(OUT/'migration_equivalence.csv',replay)
    print('DC',sum(r['qualified']=='True' for r in rows),'/',len(rows),'dual sample',duals[0],flush=True)
    d.matrix.freeze(OUT/'validation_evidence.json',[OUT/'freeze.json',OUT/'preflight_evidence.json']+list(OUT.glob('*.csv'))+[x for c in cases() for x in (LOCAL/c['key']).rglob('*') if x.is_file()])

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('build','prepare','preflight','run','analyze'))
    getattr(__import__(__name__),parser.parse_args().action)()
