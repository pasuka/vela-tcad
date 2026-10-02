"""Independent full-row column FD audit, including all local neighboring columns."""
import argparse
import copy
import json
from pathlib import Path
import subprocess
import numpy as np
import validate_simplemos_local_conservative_flux as v

a=v.a;d=v.prior
LOCAL=d.REPO/'build-release/simplemos_jacobian_columns_20260906'
OUT=d.ROOT/'jacobian_columns_20260906'
RUNNER=LOCAL/'column_runner.exe'
FUNCTION=r'''
nlohmann::json runIndependentColumns(const std::string& file, const nlohmann::json& cfg) {
    const auto problem=loadNewtonProblem(file,cfg);
    const auto solver=makeNewtonSolver(problem);
    const auto state=readExternalState(configDirectory(file),cfg,problem.mesh.numNodes());
    const auto x=solver.packArclengthState(state);
    const auto system=solver.makeArclengthSystem("drain");
    const double bias=problem.biases.at("drain");
    const double scale=solver.evaluateResidual(state).potentialScale;
    const auto J=system.jacobian(x,bias);
    const int N=static_cast<int>(problem.mesh.numNodes());
    std::ofstream out(cfg.at("output_csv").get<std::string>());
    out<<std::setprecision(17)<<"column_node,column_block,row_block,step_V,analytic_norm,fd_norm,diff_norm,max_diff,max_diff_node,outside_pattern_norm,fd_step_difference_norm,fd_roundoff_bound_norm\n";
    const vela::VectorXd base=system.residual(x,bias);
    std::ofstream raw(cfg.at("output_csv").get<std::string>()+".residual.csv");
    raw<<std::setprecision(17);
    for(int i=0;i<base.size();++i) raw<<i<<','<<base(i)<<'\n';
    int count=0;
    for(int node:cfg.at("column_nodes").get<std::vector<int>>()) for(int cb=0;cb<3;++cb) {
        const int col=cb*N+node;
        const vela::VectorXd jc=J.col(col);
        vela::VectorXd previous;
        for(double step:cfg.at("steps_V").get<std::vector<double>>()) {
            const double h=step/scale;
            vela::VectorXd xp=x,xm=x;xp(col)+=h;xm(col)-=h;
            const vela::VectorXd rp=system.residual(xp,bias),rm=system.residual(xm,bias);
            const vela::VectorXd fd=(rp-rm)/(2*h);
            for(int rb=0;rb<3;++rb) {
                const vela::VectorXd exact=jc.segment(rb*N,N),ref=fd.segment(rb*N,N);
                const vela::VectorXd diff=exact-ref;
                Eigen::Index worst;const double maximum=diff.cwiseAbs().maxCoeff(&worst);
                long double missing=0,roundoff=0;
                for(int row=rb*N;row<(rb+1)*N;++row) {
                    if(J.coeff(row,col)==0.0) missing+=static_cast<long double>(fd(row))*fd(row);
                    const long double noise=64*std::numeric_limits<double>::epsilon()*(std::abs(rp(row))+std::abs(rm(row)))/(2*h);
                    roundoff+=noise*noise;
                }
                const double stability=previous.size()==0?-1.0:(fd.segment(rb*N,N)-previous.segment(rb*N,N)).norm();
                out<<node<<','<<cb<<','<<rb<<','<<step<<','<<exact.norm()<<','<<ref.norm()<<','<<diff.norm()<<','<<maximum<<','<<worst<<','<<std::sqrt(missing)<<','<<stability<<','<<std::sqrt(roundoff)<<'\n';++count;
            }
            previous=fd;
        }
    }
    return {{"read_only",true},{"nodes",N},{"rows",count},{"nonlinear_solves",0},{"jacobian_nonzeros",J.nonZeros()},{"base_residual_norm",base.norm()}};
}
'''


def env(enabled=True,step=1e-7):
    result=d.matrix.chain.env(enabled)
    result['VELA_LINEAR_SOLVER']='sparselu'
    result['VELA_VALIDATE_VECTOR_CHAIN_STEP']=str(step)
    return result


def build():
    a.verify(v.FREEZE);a.verify(v.OUT/'evidence.json')
    LOCAL.mkdir(parents=True,exist_ok=False)
    source=(d.matrix.chain.LOCAL/'CoupledDDAssembler.cpp').read_text()
    assert source.count('const Real step=1e-7;')==1
    source=source.replace('const Real step=1e-7;','const Real step=std::getenv("VELA_VALIDATE_VECTOR_CHAIN_STEP")?std::stod(std::getenv("VELA_VALIDATE_VECTOR_CHAIN_STEP")):1e-7;')
    (LOCAL/'CoupledDDAssembler.cpp').write_text(source,newline='\n')
    runner=(v.LOCAL/'runner.cpp').read_text()
    runner=runner.replace('int main(int argc, char** argv)',FUNCTION+'\nint main(int argc, char** argv)')
    marker='        } else if (type == "newton_residual_probe") {'
    assert runner.count(marker)==1
    runner=runner.replace(marker,'        } else if (type == "independent_columns") {\n            status.update(runIndependentColumns(configFile,cfg));\n'+marker)
    (LOCAL/'runner.cpp').write_text(runner,newline='\n')
    args=[str(LOCAL/Path(s).name) if Path(s).name in ('runner.cpp','CoupledDDAssembler.cpp') else s for s in a.read(v.LOCAL/'build_command.json') if Path(s).name!='ContactCurrent.cpp']
    args[-1]=str(RUNNER);a.write(LOCAL/'build_command.json',args)
    done=subprocess.run(args,env=env(),capture_output=True,text=True)
    (LOCAL/'build.log').write_text(done.stdout+done.stderr)
    if done.returncode:raise RuntimeError(done.stderr[-3000:])
    print('Built isolated column audit runner',flush=True)


def prepare():
    cases=[];files=[Path(__file__).resolve(),RUNNER,LOCAL/'runner.cpp',LOCAL/'CoupledDDAssembler.cpp',LOCAL/'build_command.json',v.FREEZE]
    for s in a.read(LOCAL/'build_command.json'):
        if Path(s).is_file():files.append(Path(s))
    targetgeo=d.matrix.spatial.m73.Geometry('n23')
    focus=np.array([targetgeo.coords[i] for i in (1091,1092,320,324)])
    for c in a.read(d.CONTRACT)['cases']:
        geo=d.matrix.spatial.m73.Geometry(c['device']);xy=np.array([geo.coords[i] for i in range(geo.count)])
        selected={int(np.argmin(np.linalg.norm(xy-f,axis=1))) for f in focus}
        mesh=a.read(d.matrix.spatial.m73.M8/'vela'/c['device']/'mesh.json')
        seed=set(selected)
        for tri in mesh['triangles']:
            if seed.intersection(tri['node_ids']):selected.update(tri['node_ids'])
        variants=[('baseline_off',False,{},1e-7),('baseline_on',True,{},1e-7)]
        if c['device']=='n23' and float(c['vd'])==1.:
            variants += [('chain_step_1e6',True,{},1e-6),('chain_step_1e8',True,{},1e-8),
                ('vector_no_surface',True,{'model':'phumob_field'},1e-7),
                ('edge_projection',True,{'high_field_gradient_discretization':'edge_projection'},1e-7),
                ('constant_vector',True,{'model':'constant_field'},1e-7),
                ('field_frozen_off',False,{'jacobian_field_derivatives':False},1e-7),
                ('field_frozen_on',True,{'jacobian_field_derivatives':False},1e-7)]
        for label,on,changes,step in variants:
            cfg=a.read(d.matrix.spatial.precision.LOCAL/c['case']/'config.json');cfg.pop('output_state_file')
            cfg.update(simulation_type='independent_columns',state_file=str(d.matrix.spatial.precision.LOCAL/c['case']/'state.csv'),
                column_nodes=sorted(selected),steps_V=[1e-5,1e-6,1e-7],output_csv=str(LOCAL/c['case']/label/'columns.csv'))
            cfg['solver']['mobility'].update(changes)
            path=LOCAL/c['case']/label/'config.json';a.write(path,cfg);files.append(path)
            files += [Path(cfg[k]) for k in ('state_file','mesh_file','materials_file','node_doping_file')]
            cases.append(dict(case=c['case'],variant=label,enabled=on,chain_step=step,config=str(path),columns=sorted(selected)))
    a.write(OUT/'contract.json',dict(status='frozen_before_execution',cases=cases,
        scope='Full residual rows for all three potential columns on every node in the full triangle neighborhoods of the two focal edges; not a global all-column certificate.',
        gates={'active_block_relative_floor':1e-10,'stable_fd_relative':1e-4,'jacobian_relative':1e-3,'roundoff_multiplier':64},
        frozen_field_cases='Intentional quasi-Newton control: compare switch on/off identity; full-residual derivative agreement not required.',
        branch_cases='Read-only derivative tests at the same imported state, not self-consistent physical model A/B.',production_changes=False))
    files.append(OUT/'contract.json')
    a.write(OUT/'freeze.json',{'input_hashes':{a.rel(p) if p.is_relative_to(d.REPO) else str(p):a.sha(p) for p in sorted(set(files))}})
    print('Frozen',len(cases),'column audit configurations',flush=True)


def run():
    a.verify(OUT/'freeze.json')
    for c in a.read(OUT/'contract.json')['cases']:
        path=Path(c['config']);status=path.with_suffix('.status.json')
        if status.exists():continue
        done=subprocess.run([str(RUNNER),'--config',str(path),'--log','off'],env=env(c['enabled'],c['chain_step']),capture_output=True,text=True)
        path.with_suffix('.stdout.txt').write_text(done.stdout);path.with_suffix('.stderr.txt').write_text(done.stderr)
        s=json.loads(done.stdout.strip().splitlines()[-1]);s['exit_code']=done.returncode;a.write(status,s)
        print(c['case'],c['variant'],done.returncode,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('build','prepare','run'))
    globals()[parser.parse_args().action]()
