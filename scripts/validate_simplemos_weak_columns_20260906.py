"""Difference physical components before assembling weak Jacobian columns."""
import argparse
from pathlib import Path
import json
import subprocess
import audit_simplemos_surface_derivative_step_20260906 as previous

a=previous.a;d=previous.m.d
LOCAL=d.REPO/'build-release/simplemos_weak_columns_20260906'
OUT=d.ROOT/'weak_columns_20260906';RUNNER=LOCAL/'runner.exe'
FUNCTION=r'''
nlohmann::json runWeakColumns(const std::string& file,const nlohmann::json& cfg) {
    const auto problem=loadNewtonProblem(file,cfg);
    const auto solver=makeNewtonSolver(problem);
    const auto state=readExternalState(configDirectory(file),cfg,problem.mesh.numNodes());
    const auto x=solver.packArclengthState(state);const int N=problem.mesh.numNodes();
    const auto sys=solver.makeArclengthSystem("drain");const double bias=problem.biases.at("drain");
    const double scale=solver.evaluateResidual(state).potentialScale;
    const auto setStep=[](const char* value){_putenv_s("VELA_VALIDATE_TRANSPORT_STEP",value);};
    setStep("1e-6");const auto J0=sys.jacobian(x,bias);
    setStep("1e-7");const auto J1=sys.jacobian(x,bias);
    setStep("1e-8");const auto J2=sys.jacobian(x,bias);setStep("1e-7");
    const auto baseTerms=solver.evaluateCarrierTermDiagnostics(state,true).rows;
    const auto baseEdges=solver.evaluateSgEdgeFluxDiagnostics(state);
    const vela::VectorXd baseR=sys.residual(x,bias);
    std::vector<long double> en(N,0),hp(N,0);
    for(const auto& e:baseEdges){en[e.node0]+=e.electronFlux;en[e.node1]-=e.electronFlux;hp[e.node0]+=e.holeFlux;hp[e.node1]-=e.holeFlux;}
    double reconstruction=0;
    for(int i=0;i<N;++i){const auto& t=baseTerms[i];
        if(t.electronContinuityActive) reconstruction=std::max(reconstruction,static_cast<double>(std::abs(en[i]+t.electronRecombination+t.electronImpact-baseR(N+i))/std::max({(long double)t.electronFluxAbsSum,std::abs((long double)t.electronRecombination),1e-300L})));
        if(t.holeContinuityActive) reconstruction=std::max(reconstruction,static_cast<double>(std::abs(hp[i]+t.holeRecombination+t.holeImpact-baseR(2*N+i))/std::max({(long double)t.holeFluxAbsSum,std::abs((long double)t.holeRecombination),1e-300L})));
    }
    if(reconstruction>1e-10)throw std::runtime_error("Physical component reconstruction failed");
    std::ofstream out(cfg.at("output_csv").get<std::string>());out<<std::setprecision(21);
    out<<"column_node,column_block,row_block,step_V,row_node,jacobian_1e6,jacobian_1e7,jacobian_1e8,component_fd,full_residual_fd,transport_fd,source_fd,charge_fd\n";
    const auto coefficients=cfg.at("poisson_hole_coefficients").get<std::vector<double>>();
    int targets=0;
    for(const auto& target:cfg.at("targets")){
        int node=target.at("node"),cb=target.at("column_block"),rb=target.at("row_block"),col=cb*N+node;
        for(double step:cfg.at("steps_V").get<std::vector<double>>()){
            double h=step/scale;vela::VectorXd xp=x,xm=x;xp(col)+=h;xm(col)-=h;
            const auto sp=solver.unpackArclengthState(xp),sm=solver.unpackArclengthState(xm);
            const auto tp=solver.evaluateCarrierTermDiagnostics(sp,true).rows,tm=solver.evaluateCarrierTermDiagnostics(sm,true).rows;
            const vela::VectorXd rp=sys.residual(xp,bias),rm=sys.residual(xm,bias);
            std::vector<long double> transport(N,0),source(N,0),charge(N,0);
            if(rb==0){
                if(cb!=2)throw std::runtime_error("Poisson component support restricted to phi_p columns");
                for(int i=0;i<N;++i)charge[i]=-static_cast<long double>(coefficients[i])*(static_cast<long double>(tp[i].holeDensity_m3)-tm[i].holeDensity_m3)/(2*h);
            }else{
                const auto ep=solver.evaluateSgEdgeFluxDiagnostics(sp),em=solver.evaluateSgEdgeFluxDiagnostics(sm);
                if(ep.size()!=em.size())throw std::runtime_error("Edge support changed");
                for(std::size_t k=0;k<ep.size();++k){
                    if(ep[k].edgeId!=em[k].edgeId)throw std::runtime_error("Edge order changed");
                    long double df=rb==1?(static_cast<long double>(ep[k].electronFlux)-em[k].electronFlux):(static_cast<long double>(ep[k].holeFlux)-em[k].holeFlux);
                    df/=2*h;transport[ep[k].node0]+=df;transport[ep[k].node1]-=df;
                }
                for(int i=0;i<N;++i){
                    const bool active=rb==1?tp[i].electronContinuityActive:tp[i].holeContinuityActive;
                    source[i]=(rb==1?(static_cast<long double>(tp[i].electronRecombination)-tm[i].electronRecombination):(static_cast<long double>(tp[i].holeRecombination)-tm[i].holeRecombination))/(2*h);
                    if(!active){transport[i]=0;source[i]=(static_cast<long double>(rp(rb*N+i))-rm(rb*N+i))/(2*h);}
                }
            }
            for(int i=0;i<N;++i){
                const long double ref=transport[i]+source[i]+charge[i];
                const long double raw=(static_cast<long double>(rp(rb*N+i))-rm(rb*N+i))/(2*h);
                if(ref==0&&raw==0&&J0.coeff(rb*N+i,col)==0&&J1.coeff(rb*N+i,col)==0&&J2.coeff(rb*N+i,col)==0)continue;
                out<<node<<','<<cb<<','<<rb<<','<<step<<','<<i<<','<<J0.coeff(rb*N+i,col)<<','<<J1.coeff(rb*N+i,col)<<','<<J2.coeff(rb*N+i,col)<<','<<ref<<','<<raw<<','<<transport[i]<<','<<source[i]<<','<<charge[i]<<'\n';
            }
        }++targets;
    }
    return {{"read_only",true},{"targets",targets},{"component_reconstruction_scaled_error",reconstruction},{"nonlinear_solves",0}};
}
'''


def build():
    a.verify(previous.OUT/'freeze.json');LOCAL.mkdir(parents=True,exist_ok=False)
    source=(previous.m.LOCAL/'runner.cpp').read_text()
    source=source.replace('int main(int argc, char** argv)',FUNCTION+'\nint main(int argc, char** argv)')
    marker='        } else if (type == "independent_columns") {'
    assert source.count(marker)==1
    source=source.replace(marker,'        } else if (type == "weak_columns") {\n            status.update(runWeakColumns(configFile,cfg));\n'+marker)
    (LOCAL/'runner.cpp').write_text(source,newline='\n')
    args=[str(LOCAL/'runner.cpp') if Path(s).name=='runner.cpp' else s for s in a.read(previous.LOCAL/'build_command.json')]
    args[-1]=str(RUNNER);a.write(LOCAL/'build_command.json',args)
    result=subprocess.run(args,env=previous.m.env(),capture_output=True,text=True);(LOCAL/'build.log').write_text(result.stdout+result.stderr)
    assert result.returncode==0,result.stderr[-3000:]
    print('Built independent component-difference runner',flush=True)


def prepare():
    rows=[r for r in a.rows(previous.m.OUT/'columns.csv') if r['variant']=='baseline_on' and r['active']=='True' and r['fd_stable']=='False']
    assert len(rows)==33
    files=[Path(__file__).resolve(),RUNNER,LOCAL/'runner.cpp',previous.m.OUT/'columns.csv',previous.OUT/'freeze.json'];cases=[]
    for c in a.read(d.CONTRACT)['cases']:
        selected=[r for r in rows if r['case']==c['case']];geo=d.matrix.spatial.m73.Geometry(c['device'])
        cfg=a.read(previous.m.LOCAL/c['case']/'baseline_on/config.json')
        coef=1e6*d.fixed.Q*c['factor']*geo.volumes['all_cell'];coef[geo.contact_nodes]=0
        cfg.update(simulation_type='weak_columns',targets=[dict(node=int(r['column_node']),column_block=int(r['column_block']),row_block=int(r['row_block'])) for r in selected],
            steps_V=[1e-3,5e-4,2.5e-4,1e-4,5e-5,2.5e-5,1e-5,5e-6,2.5e-6,1e-6,5e-7,2.5e-7,1e-7,5e-8],
            poisson_hole_coefficients=coef.tolist(),output_csv=str(LOCAL/c['case']/'components.csv'))
        path=LOCAL/c['case']/'config.json';a.write(path,cfg);files.append(path);cases.append(dict(case=c['case'],config=str(path)))
        files.extend(Path(cfg[k]) for k in ('state_file','mesh_file','materials_file','node_doping_file'))
    a.write(OUT/'contract.json',dict(status='frozen_before_execution',cases=cases,targets=33,
        method='Poisson: difference hole charge before large dielectric/dopant terms. Continuity: difference corresponding edge fluxes and SRH first, then accumulate in long double. Underlying flux/mobility evaluation remains double.',
        independent_reference='Richardson R(h)=(4D(h/2)-D(h))/3, qualified by comparison with R(h/2); h in 1e-3,1e-4,1e-5,1e-6. Keep all step results. No use of the Jacobian to select the reference.',
        gates={'reference_step_relative':1e-4,'jacobian_relative':1e-3,'component_reconstruction_scaled_error':1e-10},
        jacobian_controls='Same physical residual, vector-chain enabled; original transport step 1e-6 and independent 1e-7/1e-8 controls.',
        limitations='Term-wise cancellation control is not arbitrary-precision mobility; unresolved small signals remain unqualified.',production_changes=False))
    files.append(OUT/'contract.json');d.matrix.freeze(OUT/'freeze.json',files)
    print('Frozen all 33 weak targets',flush=True)


def run():
    a.verify(OUT/'freeze.json')
    for c in a.read(OUT/'contract.json')['cases']:
        path=Path(c['config']);r=subprocess.run([str(RUNNER),'--config',str(path),'--log','off'],env=previous.m.env(),capture_output=True,text=True)
        path.with_suffix('.stdout.txt').write_text(r.stdout);path.with_suffix('.stderr.txt').write_text(r.stderr)
        s=json.loads(r.stdout.strip().splitlines()[-1]);s['exit_code']=r.returncode;a.write(path.with_suffix('.status.json'),s)
        assert r.returncode==0,s
        print(c['case'],s['targets'],'component differences complete',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','prepare','run'));globals()[p.parse_args().action]()
