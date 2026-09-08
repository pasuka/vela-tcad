"""Validate computed cell geometry and mobility against sealed native witnesses."""
import argparse,copy,math,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import validate_simplemos_permittivity_production_20260907 as prior
import audit_simplemos_masetti_box_mobility_20260907 as native

a=prior.a;d=prior.d;v=prior.v;p=prior.p
LOCAL=p.REPO/'build-release/simplemos_generated_box_mobility_20260907'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/generated_box_mobility_20260907'
RUNNER=LOCAL/'jvp_runner.exe'

def build():
    source=(p.REPO/'src/tools/vela_example_runner.cpp').read_text();marker='int main(int argc, char** argv)';assert source.count(marker)==1
    source=source.replace(marker,prior.jbuild.FUNCTION+'\n'+marker)
    marker='        } else if (type == "newton_residual_probe") {';assert source.count(marker)==1
    source=source.replace(marker,'        } else if (type == "joint_geometry_jvp") {\n            status.update(runJointGeometryJvp(configFile,cfg));\n'+marker)
    (LOCAL/'jvp_runner.cpp').write_text(source)
    cmd=list(a.read(prior.OVERLAY/'compile_commands.json')[-1]);assert 'vela_example_runner' in cmd[cmd.index('-o')+1]
    cmd=[x for x in cmd if not x.startswith('-I'+str(prior.OVERLAY/'include'))]
    cmd[cmd.index('-c')+1]=str(LOCAL/'jvp_runner.cpp');cmd[cmd.index('-o')+1]=str(LOCAL/'jvp_runner.o')
    a.write(LOCAL/'compile_command.json',cmd);r=subprocess.run(cmd,cwd=p.REPO/'build-release',env=v.environment(),capture_output=True,text=True)
    (LOCAL/'jvp_build.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-2000:]
    libs=[x for x in a.read(prior.OVERLAY/'link_command.json') if x.endswith('.a')]
    cmd=['D:/msys64/ucrt64/bin/g++.exe',str(LOCAL/'jvp_runner.o')]+libs+['-o',str(RUNNER)]
    a.write(LOCAL/'link_command.json',cmd);r=subprocess.run(cmd,env=v.environment(),capture_output=True,text=True)
    (LOCAL/'jvp_link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-2000:]

def config(c,arm):
    cfg=v.config(c)
    if arm=='legacy':return cfg
    cfg['solver']['region_resolved_interface_assembly']=dict(poisson_charge_node_volume='signed_transport',transport_edge_geometry='element_box')
    cfg['solver']['mobility']['edge_averaging']='element_box'
    geom=cfg.setdefault('mesh_geometry',{});geom['poisson_permittivity_policy']='cell_material'
    if arm=='native_cells':geom['poisson_cell_edge_coefficients']=a.read(prior.LOCAL/(c['device']+'_cell_coefficients.json'))
    elif arm=='generated':geom['cell_box_policy']='delaunay_transfer'
    else:raise ValueError(arm)
    return cfg

def prepare():
    archives=a.read(OUT/'prechange.json')['files']
    for name,expected in a.read(prior.OUT/'final_evidence.json')['input_hashes'].items():
        path=p.REPO/name
        if a.sha(path)!=expected:
            assert name in archives and archives[name]['sha256']==expected
            assert a.sha(p.REPO/archives[name]['archive'])==expected
    cases=prior.cases();files=[Path(__file__).resolve(),prior.OUT/'final_evidence.json',OUT/'prechange.json',v.RUNNER,RUNNER,p.REPO/'build-release/libvela_core.a']
    for c in cases:
        root=LOCAL/c['key'];c['new_jobs']=[];c['new_probes']=[]
        for arm in ('legacy','native_cells','generated'):
            cfg=config(c,arm);state=prior.LOCAL/c['key']/('legacy/state.csv' if arm=='legacy' else 'native_cells_native_mu/state.csv')
            dest=root/'fixed'/arm;paths=v.probes(cfg,dest,state)
            jvp=copy.deepcopy(a.read(dest/'edges.json'));jvp.update(simulation_type='joint_geometry_jvp',output_csv=str(dest/'jvp.csv'));a.write(dest/'jvp.json',jvp);paths.append(dest/'jvp.json')
            if arm!='legacy':
                probe=copy.deepcopy(a.read(dest/'edges.json'));probe.update(simulation_type='element_box_probe',output_csv=str(dest/'cells.csv'));a.write(dest/'cells.json',probe);paths.append(dest/'cells.json')
            c['new_probes'] += [str(x) for x in paths if x.stem in ('functional','edges','mobility','jvp','cells')];files+=paths+[state]
        for arm,seed in [('legacy',v.LOCAL/c['key']/'from_vela/state.csv'),
                         ('native_cells',prior.old.q.LOCAL/c['key']/'replacement/state.csv'),
                         ('generated',v.LOCAL/c['key']/'from_vela/state.csv'),
                         ('generated_native_seed',v.m.previous.LOCAL/c['key']/'native_referenced_joint/state.csv')]:
            cfg=config(c,'generated' if arm.startswith('generated') else arm);dest=root/arm
            cfg.update(state_file=str(seed),output_state_file=str(dest/'state.csv'));a.write(dest/'config.json',cfg);v.post_config(cfg,dest)
            files += [seed,dest/'config.json',dest/'all_row.json',dest/'acceptance_edges.json']+v.probes(cfg,dest/'post',dest/'state.csv')
            c['new_jobs'].append(dict(arm=arm,dest=str(dest)))
        cfg=v.config(c);files += [Path(cfg[n]) for n in ('mesh_file','materials_file','node_doping_file')]
    files += [x for base in ('src','include') for x in (p.REPO/base).rglob('*') if x.suffix in ('.cpp','.h')]
    a.write(OUT/'contract.json',dict(cases=cases,DC=32,gates=a.read(prior.OUT/'contract.json')['gates'],
        fixed='Previous eight-point 300 K Masetti/no-BGN model; signed Si Poisson charge, original SRH volumes and Vela constants. No native mobility table in new candidates.',
        source='Generated local coefficients use signed cotangents and same-region transfer on globally Delaunay interior edges. Unsupported boundary/non-Delaunay cases reject.',
        mobility='Per-node Masetti(total impurity); cell average uses M_i=1/4 sum adjacent g_ij*l_ij^2 with actual sum denominator; edge mu=sum_Si(g*mu_cell)/sum_Si(g).',
        geometry_gate='All native coefficients abs difference <=1e-10; measures abs difference <=1e-14 um2; cell mobility relative <=1e-8. No fitting.',
        equivalence_gate='Conservative edge coefficients and e/h flux relative L2 <=1e-8; fixed-state current relative <=1e-8. Residual differences <=1e-8 times block absolute edge flow norms; legacy bit identity.',
        defaults_changed=False,acceptance_changed=False))
    files += [OUT/'contract.json',LOCAL/'jvp_runner.cpp',LOCAL/'jvp_runner.o',LOCAL/'compile_command.json',LOCAL/'link_command.json']
    d.matrix.freeze(OUT/'freeze.json',files);print('Frozen 32 DC and fixed-state/cell/JVP checks',flush=True)

def execute(path):return v.execute(path,RUNNER if path.stem=='jvp' else v.RUNNER,v.environment())

def preflight():
    a.verify(OUT/'freeze.json');geometry=[];checks=[];jvps=[]
    for c in a.read(OUT/'contract.json')['cases']:
        for path in c['new_probes']:
            s=execute(Path(path));assert s['exit_code']==0,(path,s)
        geo,xy,cells,vol,K,parts,info,modified=prior.old.l.g.geometry(c['device'])
        debug=(native.prior.LOCAL/'native_raw/bundle'/c['device']/'MeasureCoefficients.debug').read_text()
        nativeM=native.geom.box.parse_debug_block(debug,'Measure');vp=info['measure_permutation']
        fields=native.prior.LOCAL/'native_exports'/c['device']/'fields'
        mus={car:{int(x['cell_id']):float(x['component0'])*1e-4 for x in a.rows(fields/(car+'Mobility_region0_cells.csv'))} for car in ('e','h')}
        base_jvp=None
        for arm in ('legacy','native_cells','generated'):
            root=LOCAL/c['key']/'fixed'/arm;ref=prior.LOCAL/c['key']/('legacy/post' if arm=='legacy' else 'native_cells_native_mu/post')
            newR=d.array(d.ordered(root/'residual.csv',geo.count),('psi_residual','phin_residual','phip_residual'))
            refR=d.array(d.ordered(ref/'residual.csv',geo.count),('psi_residual','phin_residual','phip_residual'))
            newE=a.rows(root/'edges.csv');refE=a.rows(ref/'edges.csv');flow=np.zeros((3,geo.count));flux=[]
            for car in ('electron','hole'):
                nr=np.array([float(x[car+'_flux']) for x in newE]);rr=np.array([float(x[car+'_flux']) for x in refE])
                flux.append(float(np.linalg.norm(nr-rr)/max(np.linalg.norm(rr),1e-300)))
            factor=a.read(d.matrix.p.CONTRACT)['frozen_poisson_residual_per_physical_charge']
            for e in refE:
                i,j=int(e['node0']),int(e['node1']);g=-K[i,j];value=abs(float(g)*(float(e['psi0_V'])-float(e['psi1_V'])))
                # Poisson residual is in internal scaled units; use the frozen physical conversion.
                flow[0,i]+=value*factor;flow[0,j]+=value*factor
                for b,car in ((1,'electron'),(2,'hole')):flow[b,i]+=abs(float(e[car+'_flux']));flow[b,j]+=abs(float(e[car+'_flux']))
            residual=[float(np.linalg.norm((newR-refR)[b,geo.free])/max(np.linalg.norm(flow[b,geo.free]),1e-300)) for b in range(3)]
            current=abs(a.read(root/'functional.status.json')['current_A_per_um']/a.read(ref/'functional.status.json')['current_A_per_um']-1)
            bitwise=bool(np.array_equal(newR,refR)) if arm=='legacy' else None
            jrows,groups=prior.old.jc.jvp_metrics(root/'jvp.csv',c['key'],arm,base_jvp if arm=='generated' else None);jvps+=jrows
            if arm=='native_cells':base_jvp=groups
            qualified=max(flux+residual+[current])<=1e-8 and (bitwise is not False) and all(r['qualified'] for r in jrows)
            checks.append(dict(key=c['key'],arm=arm,poisson_error_over_flow=residual[0],electron_error_over_flow=residual[1],hole_error_over_flow=residual[2],electron_flux_relative=flux[0],hole_flux_relative=flux[1],Id_relative=current,legacy_bitwise=bitwise,qualified=qualified))
            if arm!='legacy':
                cg=cm=mu=0.
                for r in a.rows(root/'cells.csv'):
                    cid,k,node=int(r['cell_id']),int(r['local_vertex']),int(r['node_id']);ns=cells[cid]['nodes'];slot=ns.index(node)
                    # C++ local edge ordering is frozen by the actual mesh input.
                    mesh=a.read(Path(config(c,arm)['mesh_file'])) if not hasattr(geo,'new_mesh') else geo.new_mesh;geo.new_mesh=mesh
                    pair=tuple(sorted((node,mesh['triangles'][cid]['node_ids'][(k+1)%3])))
                    expected=next(x['coefficient'] for x in parts[pair] if x['cell']==cid)
                    cg=max(cg,abs(float(r['coefficient_next'])-expected));cm=max(cm,abs(float(r['measure_m2'])*1e12-nativeM[cid]['values'][vp[slot]]))
                    if cells[cid]['material']=='Si':
                        for car,label in (('e','electron'),('h','hole')):mu=max(mu,abs(float(r[label+'_mobility_m2_V_s'])/mus[car][cid]-1))
                geometry.append(dict(key=c['key'],arm=arm,coefficient_max_absolute=cg,measure_max_absolute_um2=cm,cell_mu_max_relative=mu,qualified=cg<=1e-10 and cm<=1e-14 and mu<=1e-8))
        print('Preflight',c['key'],checks[-1]['qualified'],geometry[-1],flush=True)
    a.write_csv(OUT/'preflight.csv',checks);a.write_csv(OUT/'cell_identity.csv',geometry);a.write_csv(OUT/'jvp.csv',jvps)
    d.matrix.freeze(OUT/'preflight_evidence.json',[OUT/'freeze.json',OUT/'preflight.csv',OUT/'cell_identity.csv',OUT/'jvp.csv']+[x for c in prior.cases() for x in (LOCAL/c['key']/'fixed').rglob('*') if x.is_file()])
    assert all(x['qualified'] for x in checks+geometry)

def run():
    a.verify(OUT/'freeze.json');a.verify(OUT/'preflight_evidence.json');assert all(x['qualified']=='True' for x in a.rows(OUT/'preflight.csv')+a.rows(OUT/'cell_identity.csv'))
    def one(c):
        rows=[]
        for j in c['new_jobs']:
            dest=Path(j['dest']);s=execute(dest/'config.json');row=dict(key=c['key'],arm=j['arm'],qualified=False)
            if (dest/'state.csv').exists():
                for name in ('all_row','acceptance_edges'):execute(dest/(name+'.json'))
                row.update(prior.old.qualify(c,dest))
                for name in ('functional','edges','terms','mobility'):execute(dest/'post'/(name+'.json'))
            else:row['failure']=s.get('failure_reason','missing state')
            rows.append(row);print(c['key'],j['arm'],row['qualified'],row.get('max_row_ratio'),flush=True)
        return rows
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(one,a.read(OUT/'contract.json')['cases']) for r in group]
    keys=list(dict.fromkeys(k for r in rows for k in r));a.write_csv(OUT/'dc.csv',[{k:r.get(k,'') for k in keys} for r in rows])

def analyze():
    a.verify(OUT/'freeze.json');dc=a.rows(OUT/'dc.csv');comparison=[];dual=[];states=[];ports=[]
    for c in prior.cases():
        geo,mask=v.m.previous.prior.support(c);rows={x['arm']:x for x in dc if x['key']==c['key']}
        for arm,r in rows.items():
            Id=float(r['current_A_per_um']);comparison.append(dict(key=c['key'],device=c['device'],vg=c['vg'],vd=c['vd'],arm=arm,Id_A_per_um=Id,native_relative_error=Id/c['native_Id_A_per_um']-1,qualified=r['qualified']=='True'))
            s=a.read(LOCAL/c['key']/arm/'post/functional.status.json');e=abs(s['current_A_per_um']/s['contact_current_extractor_A_per_um']-1)
            ports.append(dict(key=c['key'],arm=arm,relative=e,qualified=e<=1e-8))
        get=lambda path:d.ordered(path/'state.csv',geo.count)
        first=get(LOCAL/c['key']/'generated');second=get(LOCAL/c['key']/'generated_native_seed')
        diff=v.m.previous.prior.delta_states(first,second,mask);error=abs(float(rows['generated']['current_A_per_um'])/float(rows['generated_native_seed']['current_A_per_um'])-1)
        gate=lambda x:max(x[k] for k in ('psi_max_V','phin_max_V','phip_max_V'))<=1e-6 and x['density_max_relative']<=1e-4
        dual.append(dict(key=c['key'],**diff,Id_relative=error,qualified=gate(diff) and error<=1e-6))
        oldstate=get(prior.LOCAL/c['key']/'native_cells_native_mu');diff=v.m.previous.prior.delta_states(first,oldstate,mask)
        oldId=a.read(prior.LOCAL/c['key']/'native_cells_native_mu/config.status.json')['contact_currents_A_per_um']['drain'];error=abs(float(rows['generated']['current_A_per_um'])/oldId-1)
        states.append(dict(key=c['key'],**diff,Id_relative=error,qualified=gate(diff) and error<=1e-8))
    for name,data in (('comparison',comparison),('dual',dual),('state_equivalence',states),('ports',ports)):a.write_csv(OUT/(name+'.csv'),data)
    summary=dict(DC=len(dc),qualified_DC=sum(x['qualified']=='True' for x in dc),qualified_dual=sum(x['qualified'] for x in dual),qualified_migration=sum(x['qualified'] for x in states),qualified_ports=sum(x['qualified'] for x in ports),
        maximum_row_ratio=max(float(x['max_row_ratio']) for x in dc),maximum_KCL_over_Id=max(float(x['kcl_over_Id']) for x in dc),maximum_migration_Id_relative=max(x['Id_relative'] for x in states))
    a.write(OUT/'summary.json',summary)
    d.matrix.freeze(OUT/'validation_evidence.json',[OUT/'freeze.json',OUT/'preflight_evidence.json',OUT/'summary.json']+list(OUT.glob('*.csv'))+[x for c in prior.cases() for x in (LOCAL/c['key']).rglob('*') if x.is_file()]);print(summary,flush=True)

if __name__=='__main__':
    arg=argparse.ArgumentParser();arg.add_argument('action',choices=('build','prepare','preflight','run','analyze'));getattr(__import__(__name__),arg.parse_args().action)()
