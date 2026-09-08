"""Isolated exact native Si transport/Poisson charge-volume coefficient A/B."""
from pathlib import Path
import subprocess
import prepare_simplemos_masetti_runtime_20260907 as p
import build_simplemos_distributed_flux_20260907 as previous

a=p.a;d=p.d;LOCAL=p.REPO/'build-release/simplemos_joint_geometry_20260907'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/joint_geometry_20260907'
RUNNER=LOCAL/'runner.exe'
HEADER=p.REPO/'scripts/diagnostics/simplemos_native_geometry.hpp'

FUNCTION=r'''
nlohmann::json runJointGeometryJvp(const std::string& configFile,const nlohmann::json& cfg) {
    const auto problem=loadNewtonProblem(configFile,cfg);
    const auto solver=makeNewtonSolver(problem);
    const auto state=readExternalState(configDirectory(configFile),cfg,problem.mesh.numNodes());
    const auto x=solver.packArclengthState(state);
    const double potential=solver.evaluateResidual(state).potentialScale;
    const int n=static_cast<int>(problem.mesh.numNodes());
    const auto system=solver.makeArclengthSystem("drain");
    const auto J=system.jacobian(x,problem.biases.at("drain"));
    const auto contacts=contactNodeMask(problem.mesh);
    std::ofstream out(cfg.at("output_csv").get<std::string>());
    if(!out)throw std::runtime_error("Cannot write joint geometry JVP");
    out<<std::setprecision(17)<<"direction,step_V,block,node_id,analytic,fd\n";
    for(int mode=0;mode<3;++mode) {
        vela::VectorXd v=vela::VectorXd::Zero(3*n);
        for(int i=0;i<n;++i)if(!contacts[i])v(mode*n+i)=std::sin(.37*(i+1));
        const vela::VectorXd analytic=J*v;
        for(double hV:{1e-5,5e-6,2.5e-6}) {
            const double h=hV/potential;
            const vela::VectorXd fd=(system.residual(x+h*v,problem.biases.at("drain"))-system.residual(x-h*v,problem.biases.at("drain")))/(2*h);
            for(int block=0;block<3;++block)for(int i=0;i<n;++i)
                out<<mode<<','<<hV<<','<<block<<','<<i<<','<<analytic(block*n+i)<<','<<fd(block*n+i)<<'\n';
        }
    }
    return {{"read_only",true},{"node_count",n},{"potential_scale_V",potential},{"rows",27*n}};
}
'''

def env(root,ratios=None,T=0.,P=0.):
    e=p.prior.prev.env('masetti_refined',root)
    for k in list(e):
        if k.startswith(('VELA_DIAGNOSTIC_DISTRIBUTED_','VELA_NATIVE_GEOMETRY_')):e.pop(k)
    e['VELA_NATIVE_GEOMETRY_T']=format(T,'.17g');e['VELA_NATIVE_GEOMETRY_P']=format(P,'.17g')
    if ratios:
        e['VELA_NATIVE_GEOMETRY_EDGE_RATIOS']=str(ratios/'edges.txt')
        e['VELA_NATIVE_GEOMETRY_NODE_RATIOS']=str(ratios/'nodes.txt')
    return e

def main():
    a.verify(p.OUT/'validation_evidence.json');LOCAL.mkdir(parents=True,exist_ok=False)
    args=a.read(previous.LOCAL/'build_command.json')
    assembler=Path(next(x for x in args if Path(x).name=='CoupledDDAssembler.cpp'))
    contact=Path(next(x for x in args if Path(x).name=='ContactCurrent.cpp'))
    for path in (assembler,contact):
        text=path.read_text();assert text.count(previous.HELPER)==1
        text=text.replace(previous.HELPER,'')
        text=text.replace('            nFlux += diagnosticConstantElectronEdgeFlux(e);\n','')
        # Remove the preceding source hook from the contact copy as well.
        if path==contact:
            import validate_simplemos_local_conservative_flux as old
            original=(p.REPO/'src/post/ContactCurrent.cpp').read_text()
            text=original
            marker='''{
    if (bandgapNarrowingConfig.equalNiFluxEvaluation !='''
            patch='''{
    simplemos_native_geometry::transport(couple_);
    if (bandgapNarrowingConfig.equalNiFluxEvaluation !='''
        else:
            assert 'diagnosticConstantElectronEdgeFlux' not in text
            marker='''{
    if (bandgapNarrowingConfig_.equalNiFluxEvaluation !='''
            patch='''{
    simplemos_native_geometry::transport(couple_);
    simplemos_native_geometry::poissonVolume(poissonElectronVol_);
    simplemos_native_geometry::poissonVolume(poissonHoleVol_);
    simplemos_native_geometry::poissonVolume(poissonDopantVol_);
    if (bandgapNarrowingConfig_.equalNiFluxEvaluation !='''
        assert text.count(marker)==1
        (LOCAL/path.name).write_text('#include "simplemos_native_geometry.hpp"\n'+text.replace(marker,patch),newline='\n')
    source=Path(next(x for x in args if Path(x).name=='runner.cpp'));text=source.read_text()
    marker='int main(int argc, char** argv)';assert text.count(marker)==1;text=text.replace(marker,FUNCTION+'\n'+marker)
    marker='        } else if (type == "newton_residual_probe") {';assert text.count(marker)==1
    text=text.replace(marker,'        } else if (type == "joint_geometry_jvp") {\n            status.update(runJointGeometryJvp(configFile,cfg));\n'+marker)
    (LOCAL/'runner.cpp').write_text(text,newline='\n')
    args=[str(LOCAL/Path(x).name) if x in (str(assembler),str(contact),str(source)) else x for x in args];args[-1]=str(RUNNER)
    a.write(LOCAL/'build_command.json',args)
    result=subprocess.run(args,env=env(LOCAL),capture_output=True,text=True)
    (LOCAL/'build.log').write_text(result.stdout+result.stderr);assert result.returncode==0,result.stderr[-4000:]
    print('Built isolated native geometry runner with unchanged solver/physics',flush=True)

if __name__=='__main__':main()
