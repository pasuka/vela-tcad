"""Frozen read-only Jv / residual FD audit along actual channel DC response."""
import argparse
from decimal import Decimal
import json
import math
from pathlib import Path
import subprocess
import numpy as np
import validate_simplemos_vela_channel_charge as prior
import seal_simplemos_vela_channel_charge as seal_prior

REPO, ROOT = prior.REPO, prior.ROOT
LOCAL = REPO/'build-release/simplemos_actual_direction_jvp'
OUT = ROOT/'actual_direction_jvp'
RUNNER = LOCAL/'direction_runner.exe'
CONTRACT = OUT/'contract.json'
FREEZE = OUT/'freeze.json'
DOC = REPO/'docs/validation/simplemos_actual_direction_jvp_2026-09-05.md'

FUNCTION = r'''
nlohmann::json runActualDirection(const std::string& configFile, const nlohmann::json& cfg)
{
    const auto cfgDir = configDirectory(configFile);
    const auto problem = loadNewtonProblem(configFile, cfg);
    const auto solver = makeNewtonSolver(problem);
    const auto state = readExternalState(cfgDir, cfg, problem.mesh.numNodes());
    const auto x = solver.packArclengthState(state);
    const double scale = solver.evaluateResidual(state).potentialScale;
    const int n = static_cast<int>(problem.mesh.numNodes());
    vela::VectorXd v(x.size());
    const std::vector<std::string> blocks = {"psi", "phin", "phip"};
    for (int b=0;b<3;++b) {
        const auto a = readNodeScalarCsv(cfg.at("direction_files").at(blocks[b]).get<std::string>(), problem.mesh.numNodes());
        for(int i=0;i<n;++i) v(b*n+i)=a[i]/scale;
    }
    const auto system = solver.makeArclengthSystem("drain");
    const double bias = problem.biases.at("drain");
    const auto J = system.jacobian(x,bias);
    auto pcfg=cfg; pcfg["state_file"]=cfg.at("plus_state_file");
    auto mcfg=cfg; mcfg["state_file"]=cfg.at("minus_state_file");
    const auto xp=solver.packArclengthState(readExternalState(cfgDir,pcfg,problem.mesh.numNodes()));
    const auto xm=solver.packArclengthState(readExternalState(cfgDir,mcfg,problem.mesh.numNodes()));
    const vela::VectorXd actualFD=(system.residual(xp,bias)-system.residual(xm,bias))/2.0;
    std::ofstream out(cfg.at("output_csv").get<std::string>());
    if(!out) throw std::runtime_error("Cannot write Jv rows");
    out<<std::setprecision(17)<<"mode,step,row_block,node_id,analytic,fd,actual_endpoint_fd,direction_scaled\n";
    for(int mode=-1;mode<3;++mode) {
        vela::VectorXd d=v;
        if(mode>=0) for(int b=0;b<3;++b) if(b!=mode) d.segment(b*n,n).setZero();
        std::vector<long double> sums(x.size(),0.L);
        for(int col=0;col<J.outerSize();++col)
            for(vela::SparseMatrixd::InnerIterator e(J,col);e;++e)
                sums[e.row()]+=static_cast<long double>(e.value())*static_cast<long double>(d(e.col()));
        for(double h : {1.,0.5,0.25}) {
            const vela::VectorXd fd=(system.residual(x+h*d,bias)-system.residual(x-h*d,bias))/(2*h);
            if(!fd.allFinite()) throw std::runtime_error("Nonfinite actual direction FD");
            for(int b=0;b<3;++b) for(int i=0;i<n;++i) {
                const int k=b*n+i;
                out<<(mode<0?"all":blocks[mode])<<','<<h<<','<<blocks[b]<<','<<i<<','
                   <<static_cast<double>(sums[k])<<','<<fd(k)<<','<<actualFD(k)<<','<<d(k)<<'\n';
            }
        }
    }
    return {{"read_only",true},{"node_count",n},{"potential_scale_V",scale},
        {"direction_max_V",scale*v.lpNorm<Eigen::Infinity>()},{"row_count",36*n},
        {"nonlinear_solves",0}};
}
'''


def build():
    if RUNNER.exists(): raise FileExistsError(RUNNER)
    LOCAL.mkdir(parents=True,exist_ok=True)
    source=(REPO/'src/tools/vela_example_runner.cpp').read_text()
    marker='int main(int argc, char** argv)'
    assert source.count(marker)==1
    source=source.replace(marker,FUNCTION+'\n'+marker)
    marker='        } else if (type == "newton_residual_probe") {'
    assert source.count(marker)==1
    source=source.replace(marker,'        } else if (type == "actual_direction_jvp") {\n            status.update(runActualDirection(configFile, cfg));\n'+marker)
    (LOCAL/'runner.cpp').write_text(source,newline='\n')
    args=prior.audit.read(REPO/'build-release/simplemos_strict_rejection_trace/build_command.json')
    args=[str(LOCAL/'runner.cpp') if a.endswith('runner.cpp') else a for a in args]
    args[-1]=str(RUNNER)
    prior.audit.write(LOCAL/'build_command.json',args)
    p=subprocess.run(args,env=prior.env(),capture_output=True,text=True)
    (LOCAL/'build.log').write_text(p.stdout+p.stderr)
    if p.returncode: raise RuntimeError(p.stderr[-3000:])
    print('Built read-only runner; original solver/core retained',flush=True)


def prepare():
    seal_prior.verify()
    if CONTRACT.exists(): raise FileExistsError(CONTRACT)
    inputs=[Path(__file__),RUNNER,LOCAL/'runner.cpp',LOCAL/'build_command.json',prior.OUT/'evidence.json',prior.FREEZE]
    for a in prior.audit.read(LOCAL/'build_command.json'):
        if Path(a).is_file(): inputs.append(Path(a))
    cases=[]
    for suffix in ('1e13','5e12'):
        dest=LOCAL/suffix
        p=prior.LOCAL/('plus_'+suffix)/'state.csv'; m=prior.LOCAL/('minus_'+suffix)/'state.csv'
        plus=prior.audit.rows(p); minus=prior.audit.rows(m)
        assert len(plus)==len(minus)
        direction={k:[] for k in ('psi','phin','phip')}
        for a,b in zip(plus,minus):
            assert a['node_id']==b['node_id']
            for k,prefix in (('psi',None),('phin','electron'),('phip','hole')):
                if prefix:
                    ref=prefix+'_qf_reference_V'; inc=prefix+'_qf_increment_V'
                    value=(Decimal(a[ref])-Decimal(b[ref])+Decimal(a[inc])-Decimal(b[inc]))/2
                else: value=(Decimal(a[k])-Decimal(b[k]))/2
                direction[k].append({'node_id':int(a['node_id']),'component0':float(value)})
        paths={}
        for k,rows in direction.items():
            path=dest/(k+'.csv');prior.audit.write_csv(path,rows);paths[k]=str(path);inputs.append(path)
        cfg=prior.audit.read(prior.BASE/'config.json');cfg.pop('output_state_file')
        cfg.update(simulation_type='actual_direction_jvp',state_file=str(prior.BASE/'state.csv'),plus_state_file=str(p),minus_state_file=str(m),direction_files=paths,output_csv=str(dest/'jvp.csv'))
        prior.audit.write(dest/'config.json',cfg)
        inputs += [p,m,dest/'config.json',Path(cfg['state_file'])]+[Path(cfg[k]) for k in ('mesh_file','materials_file','node_doping_file')]
        cases.append(suffix)
    inputs += [prior.native.upstream.LOCAL/'m65_n23_vd_0p050000_endpoint/adjoint.csv',prior.OUT/'stationarity_projection.csv',prior.OUT/'actual_direction_gradient.csv',prior.native.OUT/'node_mask.csv']
    prior.audit.write(CONTRACT,{'status':'frozen_before_execution','cases':cases,'base':'Same strict n23 Vd=.05 V Vg=.9 V baseline; complete physical psi/phin/phip actual centered DC response, QF differences formed using reference+increment before conversion to double.',
        'steps':[1,.5,.25],'directions':['all','psi','phin','phip'],
        'operator':'Original unmodified core. Jv accumulated in long double. FD at x0 +/- h*v; separate original +/- endpoint secant. No charge source needed because it is constant in state.',
        'gates':{'weighted_fd_step_relative':.001,'weighted_recovery_of_prior_gap_relative':.001,'gradient_jvp_identity_relative':1e-8},
        'interpretation':'Qualify defect measurement and localization, not the erroneous Jacobian. Tiny hole signals may be noise and must be reported without claiming local relative precision.',
        'new_nonlinear_solves':0,'new_sentaurus_runs':0,'m82_released':False,'m83_released':False})
    inputs.append(CONTRACT)
    prior.audit.write(FREEZE,{'input_hashes':{(prior.audit.rel(p) if p.is_relative_to(REPO) else str(p)):prior.audit.sha(p) for p in sorted(set(inputs))}})


def run():
    prior.audit.verify(FREEZE)
    for case in prior.audit.read(CONTRACT)['cases']:
        dest=LOCAL/case
        p=subprocess.run([str(RUNNER),'--config',str(dest/'config.json'),'--log','off'],env=prior.env(),capture_output=True,text=True)
        (dest/'stdout.txt').write_text(p.stdout);(dest/'stderr.txt').write_text(p.stderr)
        if p.returncode: raise RuntimeError(p.stderr[-3000:]+p.stdout[-1000:])
        status=json.loads(p.stdout.strip().splitlines()[-1]);prior.audit.write(dest/'status.json',status)
        print(case,status,flush=True)


def analyze():
    prior.audit.verify(FREEZE)
    adj={int(r['node_id']):r for r in prior.audit.rows(prior.native.upstream.LOCAL/'m65_n23_vd_0p050000_endpoint/adjoint.csv')}
    mask={int(r['node_id']):r['selected']=='True' for r in prior.audit.rows(prior.native.OUT/'node_mask.csv')}
    weights={'psi':'lambda_poisson','phin':'lambda_electron','phip':'lambda_hole'}
    rows=[];summaries=[];top=[];checks=[]
    cfg=prior.audit.read(CONTRACT)
    for case in cfg['cases']:
        raw=prior.audit.rows(LOCAL/case/'jvp.csv')
        assert len(raw)==36*len(adj)
        for r in raw:
            i=int(r['node_id']);w=float(adj[i][weights[r['row_block']]])
            a,f,e=(float(r[k]) for k in ('analytic','fd','actual_endpoint_fd'))
            assert all(math.isfinite(x) for x in (a,f,e,w))
            rows.append({'case':case,**r,'channel':mask[i],'x_um':float(adj[i]['x_m'])*1e6,'y_um':float(adj[i]['y_m'])*1e6,'lambda':w,'weighted_defect_A_per_um':w*(a-f),'weighted_endpoint_defect_A_per_um':w*(a-e)})
        for mode in cfg['directions']:
            for h in cfg['steps']:
                group=[r for r in rows if r['case']==case and r['mode']==mode and float(r['step'])==h]
                for block in ('all','psi','phin','phip'):
                    for region in ('all','channel','outside_channel'):
                        part=[r for r in group if (block=='all' or r['row_block']==block) and (region=='all' or r['channel']==(region=='channel'))]
                        av=np.array([float(r['analytic']) for r in part]);fv=np.array([float(r['fd']) for r in part])
                        summaries.append({'case':case,'column_mode':mode,'step':h,'row_block':block,'region':region,
                            'analytic_norm':float(np.linalg.norm(av)),'fd_norm':float(np.linalg.norm(fv)),
                            'difference_norm':float(np.linalg.norm(av-fv)),
                            'weighted_defect_A_per_um':math.fsum(r['weighted_defect_A_per_um'] for r in part),
                            'sum_abs_weighted_defect_A_per_um':math.fsum(abs(r['weighted_defect_A_per_um']) for r in part)})
        full=[r for r in rows if r['case']==case and r['mode']=='all' and float(r['step'])==.5]
        top += sorted(full,key=lambda r:abs(r['weighted_defect_A_per_um']),reverse=True)[:30]
        gradient=next(r for r in prior.audit.rows(prior.OUT/'actual_direction_gradient.csv') if float(r['amplitude_cm_3'])==float(case))
        old=next(r for r in prior.audit.rows(prior.OUT/'stationarity_projection.csv') if float(r['amplitude_cm_3'])==float(case))
        analytic=math.fsum(float(r['analytic'])*r['lambda'] for r in full)
        target=float(old['inferred_lambda_dot_J_delta_x_minus_F_delta_A_per_um'])
        endpoint=math.fsum(r['weighted_endpoint_defect_A_per_um'] for r in full)
        bystep=[next(r['weighted_defect_A_per_um'] for r in summaries if r['case']==case and r['column_mode']=='all' and r['step']==h and r['row_block']=='all' and r['region']=='all') for h in cfg['steps']]
        check={'case':case,'analytic_projection_A_per_um':analytic,'endpoint_weighted_defect_A_per_um':endpoint,'prior_gap_A_per_um':target,
            'gradient_jvp_identity_relative':abs(analytic/float(gradient['gradient_dot_actual_state_A_per_um'])-1),
            'weighted_recovery_of_prior_gap_relative':abs(endpoint/target-1),
            'weighted_fd_step_relative':max(abs(v/bystep[1]-1) for v in bystep),
            'baseline_centered_defects_A_per_um':bystep}
        check['pass']=all(check[k]<=v for k,v in cfg['gates'].items())
        checks.append(check)
    prior.audit.write_csv(OUT/'row_ledger.csv',rows);prior.audit.write_csv(OUT/'block_region_ledger.csv',summaries)
    prior.audit.write_csv(OUT/'top_rows.csv',top)
    prior.audit.write(OUT/'result.json',{'defect_measurement_qualified':all(c['pass'] for c in checks),'checks':checks,'m82_released':False,'m83_released':False})
    print(json.dumps(checks,indent=2))


def seal():
    prior.audit.verify(FREEZE)
    paths=[Path(__file__),DOC]+[p for d in (LOCAL,OUT) for p in d.rglob('*') if p.is_file()]
    prior.audit.write(OUT/'evidence.json',{'input_hashes':{prior.audit.rel(p):prior.audit.sha(p) for p in sorted(set(paths))},'new_nonlinear_solves':0,'new_sentaurus_runs':0,'m82_released':False})


def verify():
    prior.audit.verify(FREEZE);prior.audit.verify(OUT/'evidence.json');print('Actual direction Jv evidence verified')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('build','prepare','run','analyze','seal','verify'))
    {'build':build,'prepare':prepare,'run':run,'analyze':analyze,'seal':seal,'verify':verify}[p.parse_args().action]()
