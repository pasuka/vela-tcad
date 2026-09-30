"""Frozen additive pair-source response at eight SRH-volume hotspots.

The source is fixed from the accepted Vg=0 state, not reevaluated at trial
densities. This qualifies a direction; it does not implement a volume policy.
"""
import argparse
import copy
import json
import math
import shlex
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu
import validate_simplemos_phumob_lowvg_20260909 as prior
import analyze_simplemos_joint_eight_point_20260907 as fields

a, d, V = prior.a, prior.d, prior.q.run.V
REPO = prior.REPO
LOCAL = REPO/'build-release/phumob_local_source_20260909'
OUT = REPO/'reference_tcad/simplemos_sentaurus2022/phumob_local_source_20260909'
RUNNER = LOCAL/'runner.exe'
TAGS = (795, 792, 1114, 794, 791, 1115, 1203, 1188)
AMPLITUDES = (('zero',0.),('plus_full',.001),('minus_full',-.001),('plus_half',.0005),('minus_half',-.0005))
HOOK = r'''
    // Frozen diagnostic pair source. No state dependence and no Jacobian term.
    if (const char* filename = std::getenv("VELA_PAIR_SOURCE_FILE")) {
        static const auto source = [&]() {
            std::ifstream in(filename); std::size_t count, size;
            if (!(in >> size >> count) || size != mesh_.numNodes() || count != 8)
                throw std::runtime_error("Invalid pair-source manifest");
            std::vector<double> values(size, 0.); std::set<Index> seen;
            for (std::size_t k=0; k<count; ++k) {
                Index id; double value;
                if (!(in >> id >> value) || id>=size || contactNodes_[id] || ni_[id]<=0. ||
                    !std::isfinite(value) || !seen.insert(id).second)
                    throw std::runtime_error("Invalid pair-source record");
                values[id]=value;
            }
            std::string extra; if(in >> extra) throw std::runtime_error("Extra pair-source data");
            return values;
        }();
        static const double alpha = []() {
            const char* text=std::getenv("VELA_PAIR_SOURCE_ALPHA");
            if(!text) throw std::runtime_error("Missing pair-source amplitude");
            std::size_t end=0; double value=std::stod(text,&end);
            if(end!=std::string(text).size() || !std::isfinite(value) || std::abs(value)>1.)
                throw std::runtime_error("Invalid pair-source amplitude");
            return value;
        }();
        const Real sf=scaling_.enabled ? scaling_.unitSystem.continuitySourceIntegralFactor() : 1.;
        if(alpha!=0. && source[node]!=0.) rate += alpha*source[node]/(vol_[node]*sf);
    }
'''

def build():
    LOCAL.mkdir(parents=True,exist_ok=False); OUT.mkdir(parents=True,exist_ok=False)
    source=REPO/'src/equation/CoupledDDAssembler.cpp'
    text=source.read_text(encoding='utf-8')
    marker='    return rate;\n}\n\nVectorXd CoupledDDAssembler::electronDensity('
    assert text.count(marker)==1
    target=LOCAL/'CoupledDDAssembler.cpp'
    target.write_text('#include <fstream>\n#include <set>\n#include <cstdlib>\n'+text.replace(marker,HOOK+marker),encoding='utf-8')
    entry=next(x for x in a.read(REPO/'build-release/compile_commands.json') if x['file'].endswith('/CoupledDDAssembler.cpp'))
    # The include/define flags have no spaces; preserve quoted compiler paths.
    cmd=shlex.split(entry['command'].replace('\\','/'))
    cmd[cmd.index('-c')+1]=str(target);cmd[cmd.index('-o')+1]=str(LOCAL/'assembler.o')
    a.write(LOCAL/'compile_command.json',cmd)
    r=subprocess.run(cmd,cwd=entry['directory'],env=V.environment(),capture_output=True,text=True)
    (LOCAL/'compile.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    link=[cmd[0],str(LOCAL/'assembler.o'),str(REPO/'build-release/CMakeFiles/vela_example_runner.dir/src/tools/vela_example_runner.cpp.obj'),str(REPO/'build-release/libvela_core.a')]
    link += ['D:/msys64/ucrt64/lib/lib'+name for name in ('spdlog.dll.a','fmt.a','umfpack.dll.a','spqr.dll.a','cholmod.dll.a')]
    link += ['-o',str(RUNNER)]
    a.write(LOCAL/'link_command.json',link)
    r=subprocess.run(link,env=V.environment(),capture_output=True,text=True)
    (LOCAL/'link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    seal_build()

def seal_build():
    link=a.read(LOCAL/'link_command.json');source=REPO/'src/equation/CoupledDDAssembler.cpp'
    assert RUNNER.exists() and (LOCAL/'compile.log').read_text()=='' and (LOCAL/'link.log').read_text()==''
    a.write(OUT/'external_libraries.json',{x:a.sha(Path(x)) for x in link[1:] if Path(x).is_file() and not Path(x).is_relative_to(REPO)})
    # VERSION is unused in this translation unit; retain the exact resolved command.
    assert 'VELA_VERSION' not in source.read_text()
    a.write(OUT/'build_ledger.json',dict(initial_build='Compile/link succeeded; first manifest failed because external library paths cannot be repo-relative. Dependent prepare stopped before inputs were created.',repair='Hash external libraries in a separate absolute-path registry; no binary or numerical change.'))
    d.matrix.freeze(OUT/'build_evidence.json',[Path(__file__).resolve(),source,REPO/'build-release/compile_commands.json',OUT/'external_libraries.json',OUT/'build_ledger.json']+[Path(x) for x in link[1:] if Path(x).is_file() and Path(x).is_relative_to(REPO)]+list(LOCAL.glob('*')))
    print('Built current-production additive-source overlay',flush=True)

def env(case,alpha):
    e=V.environment();e.update(VELA_PAIR_SOURCE_FILE=case['source_file'],VELA_PAIR_SOURCE_ALPHA=format(alpha,'.17g'));return e

def execute(path,case,alpha):return V.execute(path,RUNNER,env(case,alpha))

def prepare():
    a.verify(OUT/'build_evidence.json');a.verify(prior.OUT/'srh_fields/srh_evidence.json')
    table=a.rows(prior.OUT/'srh_fields/srh_nodes.csv');jobs=[];scope=[];files=[OUT/'build_evidence.json',prior.OUT/'srh_fields/srh_evidence.json']
    geo23,_=fields.support(dict(device='n23'))
    for result in (prior.LOCAL/'candidate/dc/phumob').glob('*/vela/vg_000/attempt_0/result.json'):
        c=a.read(result);assert c['qualified'];c['key']=c['case'];c['baseline']=str(result.parent)
        geo,mask=fields.support(c);base=result.parent;cfg=a.read(base/'config.json')
        edge=next(r for r in a.rows(base/'acceptance_edges.csv') if abs(float(r['electron_flux']))>1e-100)
        conv=float(edge['electron_particle_line_flux_per_m_s'])/float(edge['electron_flux'])
        records=[]
        for tag in TAGS:
            distances=np.linalg.norm(np.asarray([geo.coords[i] for i in range(geo.count)])-np.asarray(geo23.coords[tag]),axis=1)
            node=int(np.argmin(distances));assert distances[node]<1e-8 and mask[node]
            r=next(x for x in table if x['key']==c['key']+'_vg_000' and int(x['node_id'])==node)
            physical=-float(r['volume_contribution_per_m_s']);internal=physical/conv
            records.append((node,internal));scope.append(dict(key=c['key'],device=c['device'],tag=tag,node=node,x_um=geo.coords[node][0],y_um=geo.coords[node][1],integrated_pair_source_per_m_s=physical,internal_source=internal,native_volume_m2=float(r['native_volume_m2']),native_rate_cm3_s=physical/(float(r['native_volume_m2'])*1e6)))
        src=LOCAL/(c['key']+'_source.txt');src.write_text(f'{geo.count} 8\n'+''.join(f'{i} {s:.17g}\n' for i,s in records));c['source_file']=str(src);c['jobs']=[]
        files += [src,result,base/'state.csv',base/'config.json']
        for label,alpha in AMPLITUDES:
            dest=LOCAL/'dc'/c['key']/label;deck=copy.deepcopy(cfg);deck.update(state_file=str(base/'state.csv'),output_state_file=str(dest/'state.csv'))
            a.write(dest/'config.json',deck);V.post_config(deck,dest)
            c['jobs'].append(dict(label=label,alpha=alpha,dest=str(dest)))
            files += list(dest.glob('*.json'))
        for label in ('zero','unit'):
            dest=LOCAL/'fixed'/c['key']/label;files+=V.probes(cfg,dest,base/'state.csv')
            deck=a.read(dest/'functional.json');deck.update(simulation_type='parameter_jacobian',output_csv=str(dest/'jacobian.csv'));a.write(dest/'jacobian.json',deck);files.append(dest/'jacobian.json')
        jobs.append(c)
    assert len(jobs)==4
    a.write_csv(OUT/'scope.csv',scope)
    a.write(OUT/'contract.json',dict(cases=jobs,source='Frozen Vela SRH rate times native-minus-Vela volume on 8 spatially mapped sites, equal integrated electron/hole source in both solvers. State-independent source; zero extra Jacobian. Original SRH/PhuMob/Poisson unchanged.',amplitudes=AMPLITUDES,gates=dict(row=1e-6,kcl=1e-8,source_relative=1e-5,linear_relative=1e-8,prediction=1e-3,two_amplitude=1e-3,even_over_odd=.01,signal_over_drift=100),native_response_qualified=False,finite_volume_replacement=False))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'scope.csv',OUT/'contract.json']);print('Frozen 4 controls / 20 DC / 8 matched sites',flush=True)

def matrix(path,N):
    rows=a.rows(path);return coo_matrix(([float(r['value']) for r in rows],([int(r['row']) for r in rows],[int(r['column']) for r in rows])),shape=(3*N,3*N)).tocsc()

def residual(path,N):return d.array(d.ordered(path/'residual.csv',N),('psi_residual','phin_residual','phip_residual')).reshape(-1)

def preflight():
    a.verify(OUT/'freeze.json');checks=[];tangents=[];scope=a.rows(OUT/'scope.csv')
    for c in a.read(OUT/'contract.json')['cases']:
        geo,mask=fields.support(c);N=geo.count;root=LOCAL/'fixed'/c['key']
        for label,alpha in (('zero',0.),('unit',1.)):
            for name in ('functional','edges','terms','jacobian'):
                s=execute(root/label/(name+'.json'),c,alpha);assert s['exit_code']==0,s
        zero=root/'zero';unit=root/'unit';R0=residual(zero,N);R1=residual(unit,N);S=np.zeros(3*N)
        for r in scope:
            if r['key']==c['key']:
                i=int(r['node']);S[N+i]=S[2*N+i]=float(r['internal_source'])
        error=np.linalg.norm(R1-R0-S)/np.linalg.norm(S)
        J=matrix(zero/'jacobian.csv',N);J1=matrix(unit/'jacobian.csv',N)
        jdiff=J1-J;identity=jdiff.nnz==0 or np.max(abs(jdiff.data))==0
        edge0=a.rows(zero/'edges.csv');edge1=a.rows(unit/'edges.csv')
        transport=edge0==edge1
        scale=1/np.maximum(np.asarray(abs(J).max(axis=1).toarray()).ravel(),1e-300);lu=splu((diags(scale)@J).tocsc());wide=J.tocsr().astype(np.longdouble);dx=lu.solve(-scale*S)
        for _ in range(4):dx-=lu.solve(scale*np.asarray(wide@dx.astype(np.longdouble)+S.astype(np.longdouble),dtype=float))
        linear=np.linalg.norm(np.asarray(wide@dx.astype(np.longdouble)+S.astype(np.longdouble),dtype=float))/np.linalg.norm(S)
        vt=a.read(zero/'jacobian.status.json')['potential_scale_V']
        for i in np.flatnonzero(mask):tangents.append(dict(key=c['key'],node=int(i),psi_V=vt*dx[i],phin_V=vt*dx[N+i],phip_V=vt*dx[2*N+i]))
        baseline=d.ordered(Path(c['baseline'])/'all_row.csv',N);terms=d.ordered(zero/'terms.csv',N)
        same=all(x[k]==y[k] for x,y in zip(baseline,terms) for k in ('electron_residual','hole_residual','electron_recombination','hole_recombination'))
        row=dict(key=c['key'],zero_identity=same,source_relative=float(error),jacobian_identity=bool(identity),transport_identity=transport,linear_relative=float(linear),qualified=bool(same and error<=1e-5 and identity and transport and linear<=1e-8));checks.append(row);print(row,flush=True)
    a.write_csv(OUT/'preflight.csv',checks);a.write_csv(OUT/'tangent.csv',tangents)
    d.matrix.freeze(OUT/'preflight_evidence.json',[OUT/'freeze.json',OUT/'preflight.csv',OUT/'tangent.csv']+[p for p in (LOCAL/'fixed').rglob('*') if p.is_file()])

def run():
    a.verify(OUT/'preflight_evidence.json');assert all(r['qualified']=='True' for r in a.rows(OUT/'preflight.csv'))
    def one(c):
        rows=[]
        for j in c['jobs']:
            dest=Path(j['dest']);s=execute(dest/'config.json',c,j['alpha'])
            for name in ('all_row','acceptance_edges'):execute(dest/(name+'.json'),c,j['alpha'])
            # qualify's cached edge probe is already generated with this source.
            result=prior.q.run.w.old.prior.old.qualify(c,dest)
            r=dict(key=c['key'],label=j['label'],alpha=j['alpha'],elapsed_seconds=s['elapsed_seconds'],**result);a.write(dest/'result.json',r);rows.append(r);print(c['key'],j['label'],result,flush=True)
        return rows
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(one,a.read(OUT/'contract.json')['cases']) for r in group]
    a.write_csv(OUT/'dc.csv',rows)
    d.matrix.freeze(OUT/'dc_evidence.json',[OUT/'preflight_evidence.json',OUT/'dc.csv']+[p for p in (LOCAL/'dc').rglob('*') if p.is_file()])

def analyze():
    a.verify(OUT/'dc_evidence.json');dc=a.rows(OUT/'dc.csv');tr=a.rows(OUT/'tangent.csv');results=[];ports=[]
    for c in a.read(OUT/'contract.json')['cases']:
        geo,mask=fields.support(c);ids=np.flatnonzero(mask);root=LOCAL/'dc'/c['key']
        def state(path):
            rows=d.ordered(path,geo.count);return [[fields.physical(rows[i],k) for i in ids] for k in ('psi','phin','phip')]
        base=state(Path(c['baseline'])/'state.csv');zero=state(root/'zero/state.csv')
        def delta(x,y):return np.array([[float(a-b) for a,b in zip(xx,yy)] for xx,yy in zip(x,y)])
        noise=np.linalg.norm(delta(zero,base));ds={l:delta(state(root/l/'state.csv'),zero) for l,_ in AMPLITUDES if l!='zero'}
        odd={amp:(ds['plus_'+amp]-ds['minus_'+amp])/2 for amp in ('full','half')}
        rows=[r for r in tr if r['key']==c['key']];tangent=np.array([[float(r[k+'_V']) for r in rows] for k in ('psi','phin','phip')]);linearity=np.linalg.norm(odd['full']-2*odd['half'])/np.linalg.norm(odd['full'])
        for amp,alpha in (('full',.001),('half',.0005)):
            pred=np.linalg.norm(odd[amp]-alpha*tangent)/np.linalg.norm(alpha*tangent);even=np.linalg.norm((ds['plus_'+amp]+ds['minus_'+amp])/2)/np.linalg.norm(odd[amp]);snr=np.linalg.norm(odd[amp])/max(noise,1e-300)
            qualified=all(r['qualified']=='True' for r in dc if r['key']==c['key']) and pred<=1e-3 and linearity<=1e-3 and even<=.01 and snr>=100
            results.append(dict(key=c['key'],amplitude=amp,prediction_relative=float(pred),two_amplitude_relative=float(linearity),even_over_odd=float(even),signal_over_drift=float(snr),qualified=bool(qualified)))
            sp=a.read(root/('plus_'+amp)/'config.status.json')['contact_currents_A_per_um'];sm=a.read(root/('minus_'+amp)/'config.status.json')['contact_currents_A_per_um']
            for port in sp:ports.append(dict(key=c['key'],amplitude=amp,port=port,derivative_A_per_um=(sp[port]-sm[port])/(2*alpha)))
    a.write_csv(OUT/'response.csv',results);a.write_csv(OUT/'ports.csv',ports)
    summary=dict(DC=len(dc),qualified_DC=sum(r['qualified']=='True' for r in dc),responses=len(results),qualified_responses=sum(r['qualified'] for r in results),native_response_qualified=False,finite_volume_replacement=False)
    a.write(OUT/'summary.json',summary);d.matrix.freeze(OUT/'evidence.json',[OUT/'dc_evidence.json',OUT/'response.csv',OUT/'ports.csv',OUT/'summary.json']);print(summary,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','seal_build','prepare','preflight','run','analyze'));globals()[p.parse_args().action]()
