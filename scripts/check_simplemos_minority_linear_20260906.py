"""Static same-matrix linear accuracy control; no production or DC acceptance change."""
import argparse
from pathlib import Path
import json
import subprocess
import audit_simplemos_minority_residual_20260906 as v

a=v.a;d=v.d;LOCAL=v.LOCAL/'linear';OUT=v.OUT/'linear';RUNNER=LOCAL/'runner.exe'
FUNCTION=r'''
nlohmann::json runMinorityLinear(const std::string& file,const nlohmann::json& cfg) {
    using MP=boost::multiprecision::cpp_dec_float_100;
    const auto problem=loadNewtonProblem(file,cfg);const auto solver=makeNewtonSolver(problem);
    const auto state=readExternalState(configDirectory(file),cfg,problem.mesh.numNodes());
    const auto x=solver.packArclengthState(state);const auto sys=solver.makeArclengthSystem("drain");const double bias=problem.biases.at("drain");
    const auto J=sys.jacobian(x,bias);const vela::VectorXd F=sys.residual(x,bias);const int N=problem.mesh.numNodes();
    const auto terms=solver.evaluateCarrierTermDiagnostics(state,true).rows;
    const double potentialScale=solver.evaluateResidual(state).potentialScale;
    const auto highResidual=[&](const vela::VectorXd& dx) {
        std::vector<MP> r(F.size());for(int i=0;i<F.size();++i)r[i]=MP(F(i));
        for(int col=0;col<J.outerSize();++col)for(vela::SparseMatrixd::InnerIterator it(J,col);it;++it)r[it.row()]+=MP(it.value())*MP(dx(it.col()));
        vela::VectorXd answer(F.size());for(int i=0;i<F.size();++i)answer(i)=static_cast<double>(r[i]);return answer;
    };
    std::ofstream out(cfg.at("output_csv").get<std::string>());out<<std::setprecision(17)<<"profile,refinement,row,carrier,node,original_residual,row_scale,step_V,linear_residual_hp,linear_ratio,trial_residual,trial_ratio,row_weight,column_weight\n";
    for(const std::string profile:{"legacy","equilibrated"}) {
        vela::VectorXd rw=vela::VectorXd::Ones(F.size()),cw=rw;
        if(profile=="legacy") {
            for(int i=0;i<N;++i)for(int b=1;b<=2;++b) {
                const auto& t=terms[i];const bool active=b==1?t.electronContinuityActive:t.holeContinuityActive;if(!active)continue;
                double source=std::abs(b==1?t.electronRecombination:t.holeRecombination),flux=b==1?t.electronFluxAbsSum:t.holeFluxAbsSum;
                if(source>=1e-18&&source>=.001*flux)rw(b*N+i)=std::clamp(1/std::max(source,1e-30),1e-12,1e12);
            }
        }else {
            vela::VectorXd maximum=vela::VectorXd::Zero(F.size());
            for(int col=0;col<J.outerSize();++col)for(vela::SparseMatrixd::InnerIterator it(J,col);it;++it)maximum(it.row())=std::max(maximum(it.row()),std::abs(it.value()));
            for(int i=0;i<F.size();++i)if(maximum(i)>0)rw(i)=std::ldexp(1.,-std::ilogb(maximum(i)));
            maximum.setZero();
            for(int col=0;col<J.outerSize();++col)for(vela::SparseMatrixd::InnerIterator it(J,col);it;++it)maximum(it.col())=std::max(maximum(it.col()),std::abs(it.value()*rw(it.row())));
            for(int i=0;i<F.size();++i)if(maximum(i)>0)cw(i)=std::ldexp(1.,-std::ilogb(maximum(i)));
        }
        vela::SparseMatrixd scaled=J;
        for(int col=0;col<scaled.outerSize();++col)for(vela::SparseMatrixd::InnerIterator it(scaled,col);it;++it)it.valueRef()*=rw(it.row())*cw(it.col());
        vela::LinearSolver linear;vela::VectorXd dx=linear.solve(scaled,-F.cwiseProduct(rw)).cwiseProduct(cw);
        for(int k=0;k<=4;++k) {
            const vela::VectorXd defect=highResidual(dx),trial=sys.residual(x+dx,bias);
            for(int b=1;b<=2;++b)for(int i=0;i<N;++i) {
                const auto& t=terms[i];if(!(b==1?t.electronContinuityActive:t.holeContinuityActive))continue;
                double scale=std::max(b==1?t.electronFluxAbsSum:t.holeFluxAbsSum,std::abs(b==1?t.electronRecombination:t.holeRecombination));
                if(scale<=0)continue;int r=b*N+i;
                out<<profile<<','<<k<<','<<r<<','<<(b==1?"electron":"hole")<<','<<i<<','<<F(r)<<','<<scale<<','<<dx(r)*potentialScale<<','<<defect(r)<<','<<std::abs(defect(r))/scale<<','<<trial(r)<<','<<std::abs(trial(r))/scale<<','<<rw(r)<<','<<cw(r)<<'\n';
            }
            if(k<4)dx+=linear.solve(scaled,-defect.cwiseProduct(rw)).cwiseProduct(cw);
        }
    }
    return {{"read_only",true},{"nonlinear_solves",0},{"matrix_rows",F.size()},{"note","Uncapped static trial, fixed original row scales; not DC qualification."}};
}
'''


def build():
    a.verify(v.OUT/'freeze.json');LOCAL.mkdir(parents=True,exist_ok=False)
    src=(v.LOCAL/'runner.cpp').read_text();src='#include "vela/solver/LinearSolver.h"\n#include <boost/multiprecision/cpp_dec_float.hpp>\n'+src
    src=src.replace('int main(int argc, char** argv)',FUNCTION+'\nint main(int argc, char** argv)')
    marker='        } else if (type == "independent_columns") {';assert src.count(marker)==1
    src=src.replace(marker,'        } else if (type == "minority_linear") {\n            status.update(runMinorityLinear(configFile,cfg));\n'+marker)
    (LOCAL/'runner.cpp').write_text(src,newline='\n')
    args=[str(LOCAL/'runner.cpp') if Path(s).name=='runner.cpp' else s for s in a.read(v.LOCAL/'build_command.json')];args[-1]=str(RUNNER)
    a.write(LOCAL/'build_command.json',args);p=subprocess.run(args,env=v.m.env(),capture_output=True,text=True);(LOCAL/'build.log').write_text(p.stdout+p.stderr)
    assert p.returncode==0,p.stderr[-3000:];print('Built static linear accuracy control',flush=True)


def retry_build():
    assert not RUNNER.exists() and not (OUT/'freeze.json').exists()
    src=(LOCAL/'runner.cpp').read_text();assert '#include "vela/solver/LinearSolver.h"' not in src
    (LOCAL/'runner_initial.cpp').write_text(src,newline='\n')
    (LOCAL/'runner.cpp').write_text('#include "vela/solver/LinearSolver.h"\n'+src,newline='\n')
    p=subprocess.run(a.read(LOCAL/'build_command.json'),env=v.m.env(),capture_output=True,text=True)
    (LOCAL/'build_retry.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr[-3000:]
    print('Built linear control after adding its explicit solver header',flush=True)


def prepare():
    jobs=[];files=[Path(__file__).resolve(),RUNNER,LOCAL/'runner.cpp',v.OUT/'freeze.json']
    for j in a.read(v.OUT/'contract.json')['probes']:
        chosen=(j['device']=='n19' and j['vd']==.05 and j['vg']==.35 and j['initialization']=='vela') or (j['device']=='n23' and j['initialization']=='vela' and ((j['vd']==.05 and j['vg']==1) or (j['vd']==1 and j['vg']==.35)))
        if not chosen:continue
        cfg=a.read(Path(j['config']));cfg.update(simulation_type='minority_linear');dest=LOCAL/Path(j['tag']);cfg['output_csv']=str(dest/'linear.csv')
        # Match the legacy row-weight rule exactly; no inferred defaults in this control.
        assert cfg['solver']['continuity_row_scaling']==dict(enabled=True,flux_fraction=.001,scale_floor=1e-30,min_source_scale=1e-18,min_weight=1e-12,max_weight=1e12)
        a.write(dest/'config.json',cfg);jobs.append(dict(tag=j['tag'],config=str(dest/'config.json')));files += [dest/'config.json',Path(cfg['state_file'])]
    assert len(jobs)==3
    a.write(OUT/'contract.json',dict(status='frozen_before_execution',jobs=jobs,
        trigger='Exact n19 rejection replay: worst hole node 41 raw step=0 and local Jdx+F approximately F, despite globally small linear norm; original row weight=1.',
        axes='Same complete J,F and Eigen SparseLU. Compare legacy row weighting with powers-of-two row/column equilibration; each gets zero through four refinement corrections using 100-digit Jdx+F.',
        scope='Three fixed states, no nonlinear solves. Uncapped x+dx trial uses frozen original row scales only, and is never a DC qualification. Same physical residual and Jacobian model.',
        old_record='Six exact trajectory replays and existing source/gates stay immutable.'))
    files.append(OUT/'contract.json');d.matrix.freeze(OUT/'freeze.json',files)


def run():
    a.verify(OUT/'freeze.json')
    for j in a.read(OUT/'contract.json')['jobs']:
        path=Path(j['config']);p=subprocess.run([str(RUNNER),'--config',str(path),'--log','off'],env=v.m.env(),capture_output=True,text=True)
        path.with_suffix('.stdout.txt').write_text(p.stdout);path.with_suffix('.stderr.txt').write_text(p.stderr)
        s=json.loads(p.stdout.strip().splitlines()[-1]);s['exit_code']=p.returncode;a.write(path.with_suffix('.status.json'),s)
        assert s['exit_code']==0,s;print(j['tag'],'linear control complete',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('build','retry_build','prepare','run'));globals()[parser.parse_args().action]()
