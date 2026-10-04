"""Build a fully recompiled, isolated constants/SRH diagnostic executable.

All 47 vela_core translation units use the same runtime constants header.
Production constants and binaries are untouched. The local volume hook applies
only after the independent Poisson volume arrays have been constructed.
"""
import copy
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import complete_simplemos_generated_box_mobility_20260907 as done

s=done.s;a=s.a;d=s.d;v=s.v;p=s.p
LOCAL=p.REPO/'build-release/simplemos_constants_srh_calibration_20260907'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/constants_srh_calibration_20260907'
RUNNER=LOCAL/'runner.exe'
HEADER=r'''
#include <cstdlib>
#include <cmath>
#include <stdexcept>
namespace vela_parameter_audit {
inline double relative(const char* name) {
    const char* text=std::getenv(name);const double value=text?std::stod(text):0.;
    if(!std::isfinite(value)||std::abs(value)>.01)throw std::runtime_error("Invalid diagnostic relative parameter");
    return value;
}
}
'''
SRH=r'''
    if(const char* nodeText=std::getenv("VELA_CANDIDATE_SRH_NODE")) {
        const auto node=static_cast<Index>(std::stoul(nodeText));
        if(node>=mesh_.numNodes() || contactNodes_[node] || ni_[node]<=0.)
            throw std::runtime_error("SRH diagnostic node must be free semiconductor");
        // In the frozen SRH-only/no-avalanche/no-carrier-floor profile this
        // affects R, its six state derivatives, and solved source diagnostics.
        vol_[node]*=1.+vela_parameter_audit::relative("VELA_CANDIDATE_SRH_RELATIVE");
    }
'''
JACOBIAN=r'''
nlohmann::json runParameterJacobian(const std::string& configFile,const nlohmann::json& cfg) {
    const auto problem=loadNewtonProblem(configFile,cfg);
    const auto solver=makeNewtonSolver(problem);
    const auto state=readExternalState(configDirectory(configFile),cfg,problem.mesh.numNodes());
    const auto x=solver.packArclengthState(state);
    const auto system=solver.makeArclengthSystem("drain");
    const auto J=system.jacobian(x,problem.biases.at("drain"));
    std::ofstream out(cfg.at("output_csv").get<std::string>());
    if(!out)throw std::runtime_error("Cannot write parameter Jacobian");
    out<<std::setprecision(17)<<"row,column,value\n";
    for(int k=0;k<J.outerSize();++k)for(vela::SparseMatrixd::InnerIterator it(J,k);it;++it)
        out<<it.row()<<','<<it.col()<<','<<it.value()<<'\n';
    return {{"potential_scale_V",solver.evaluateResidual(state).potentialScale},{"read_only",true},{"rows",J.rows()}};
}
'''

def build():
    assert not (OUT/'build_evidence.json').exists()
    LOCAL.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    header=(p.REPO/'include/vela/core/PhysicalConstants.h').read_text()
    header=header.replace('namespace vela::constants {',HEADER+'\nnamespace vela::constants {')
    for key,value in [('q','1.602176634e-19'),('kb','1.380649e-23'),('eps0','8.8541878128e-12')]:
        import re
        header,count=re.subn(r'constexpr double '+key+r'\s*=\s*'+re.escape(value)+r';',
            'inline const double '+key+' = '+value+' * (1.+vela_parameter_audit::relative("VELA_CANDIDATE_'+key.upper()+'_RELATIVE"));',header)
        assert count==1
    header=header.replace('constexpr double Vt_300','inline const double Vt_300')
    target=LOCAL/'include/vela/core/PhysicalConstants.h';target.parent.mkdir(parents=True,exist_ok=True);target.write_text(header)
    path=p.REPO/'include/vela/equation/ElementEdgeGssLauxAD.inl'
    target=LOCAL/path.relative_to(p.REPO);target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(path.read_text().replace('constexpr Real kBoltzmannEvPerK','const Real kBoltzmannEvPerK'))
    src=(p.REPO/'src/equation/CoupledDDAssembler.cpp').read_text()
    marker='    detail::validateElementBoxMobilityContext(mesh_,mobilityConfig_,scaling_.regionResolvedInterfaceAssembly);'
    assert src.count(marker)==1;src=src.replace(marker,SRH+'\n'+marker)
    (LOCAL/'CoupledDDAssembler.cpp').write_text(src)
    src=(p.REPO/'src/physics/ImpactIonizationModel.cpp').read_text().replace('constexpr Real kBoltzmann_eV_per_K','const Real kBoltzmann_eV_per_K')
    (LOCAL/'ImpactIonizationModel.cpp').write_text(src)
    runner=(s.LOCAL/'jvp_runner.cpp').read_text()
    marker='int main(int argc, char** argv)';assert runner.count(marker)==1
    runner=runner.replace(marker,JACOBIAN+'\n'+marker)
    marker='        } else if (type == "newton_residual_probe") {';assert runner.count(marker)==1
    runner=runner.replace(marker,'        } else if (type == "parameter_jacobian") {\n            status.update(runParameterJacobian(configFile,cfg));\n'+marker)
    marker='std::cout << status.dump()';assert runner.count(marker)==1
    runner=runner.replace(marker,'status["diagnostic_constants"]={{"q",vela::constants::q},{"kb",vela::constants::kb},{"eps0",vela::constants::eps0},{"Vt_300",vela::constants::Vt_300}};\n        '+marker)
    (LOCAL/'runner.cpp').write_text(runner)
    template=list(a.read(s.prior.OVERLAY/'compile_commands.json')[0]);template[1]='-I'+str(LOCAL/'include')
    jobs=[]
    entries=[x for x in a.read(p.REPO/'build-release/compile_commands.json') if '/vela_core.dir/' in x['output']]
    assert len(entries)==47
    for entry in entries:
        path=Path(entry['file']);cmd=list(template)
        src=LOCAL/path.name if path.name in ('CoupledDDAssembler.cpp','ImpactIonizationModel.cpp') else path
        obj=LOCAL/'objects'/Path(entry['output']).relative_to(p.REPO/'build-release/CMakeFiles/vela_core.dir').with_suffix('.o')
        obj.parent.mkdir(parents=True,exist_ok=True);cmd[cmd.index('-c')+1]=str(src);cmd[cmd.index('-o')+1]=str(obj);jobs.append(cmd)
    cmd=list(template);cmd[cmd.index('-c')+1]=str(LOCAL/'runner.cpp');cmd[cmd.index('-o')+1]=str(LOCAL/'runner.o');jobs.append(cmd)
    if (LOCAL/'compile_commands.json').exists():assert a.read(LOCAL/'compile_commands.json')==jobs
    else:a.write(LOCAL/'compile_commands.json',jobs)
    def compile(cmd):
        obj=Path(cmd[cmd.index('-o')+1]);log=Path(str(obj)+'.log')
        if obj.exists() and log.exists() and not log.read_text():return
        r=subprocess.run(cmd,cwd=p.REPO/'build-release',env=v.environment(),capture_output=True,text=True)
        Path(cmd[cmd.index('-o')+1]+'.log').write_text(r.stdout+r.stderr)
        assert r.returncode==0,r.stderr[-3500:]
        print('Compiled',Path(cmd[cmd.index('-c')+1]).name,flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(compile,jobs))
    libs=[x for x in a.read(s.prior.OVERLAY/'link_command.json') if x.endswith('.a') and not x.endswith('libvela_core.a')]
    cmd=['D:/msys64/ucrt64/bin/g++.exe']+[x[x.index('-o')+1] for x in jobs]+libs+['-o',str(RUNNER)]
    a.write(LOCAL/'link_command.json',cmd);r=subprocess.run(cmd,env=v.environment(),capture_output=True,text=True)
    (LOCAL/'link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    d.matrix.freeze(OUT/'build_evidence.json',[Path(__file__).resolve()]+[x for x in LOCAL.rglob('*') if x.is_file()]+[Path(x['file']) for x in entries])

if __name__=='__main__':build()
