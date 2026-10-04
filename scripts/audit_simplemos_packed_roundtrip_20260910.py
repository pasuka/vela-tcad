"""Output-only replay isolates free/contact and psi/electron/hole packing loss."""
import argparse
import shlex
import subprocess
from pathlib import Path
import validate_simplemos_production_weighted_20260910 as v
a,d=v.a,v.d
LOCAL=v.REPO/'build-release/pr_audit';OUT=v.OUT.parent/'packed_roundtrip'
RUNNER=LOCAL/'runner.exe'
HOOK=r'''
    if(const char* directory=std::getenv("VELA_PACK_AUDIT")) {
        static int sequence=0;
        const std::string prefix=std::string(directory)+"/state_"+std::to_string(sequence++);
        const auto bc=buildBoundaryConditions(assembler);
        VectorXd y=packReferencedSolution(assembler,sol,bc);
        for(const auto& [i,z]:bc.psi)y(i)=(z*potentialScale)/potentialScale;
        for(const auto& [i,z]:bc.phin)y(N+i)=(z*potentialScale)/potentialScale-assembler.electronQuasiFermiReferenceAt(i)/potentialScale;
        for(const auto& [i,z]:bc.phip)y(2*N+i)=(z*potentialScale)/potentialScale-assembler.holeQuasiFermiReferenceAt(i)/potentialScale;
        const VectorXd original=assembler.residual(x,bc);
        std::ofstream states(prefix+".csv");states<<std::setprecision(17)<<"row,block,node,contact,x,repacked,R,potential_scale\n";
        for(int j=0;j<3*N;++j) {
            const int block=j/N,node=j%N;const bool contact=block==0?bc.psi.contains(node):block==1?bc.phin.contains(node):bc.phip.contains(node);
            states<<j<<','<<block<<','<<node<<','<<contact<<','<<x(j)<<','<<y(j)<<','<<original(j)<<','<<potentialScale<<'\n';
        }
        std::ofstream log(prefix+"_norms.csv");log<<std::setprecision(17)<<"mask,changed,psi_norm,electron_norm,hole_norm,total_norm,delta_norm\n";
        for(int mask=0;mask<9;++mask) {
            VectorXd trial=x;int changed=0;
            for(int j=0;j<3*N;++j) {
                const int block=j/N,node=j%N;const bool contact=block==0?bc.psi.contains(node):block==1?bc.phin.contains(node):bc.phip.contains(node);
                if(mask==7 || (mask>=1 && mask<=6 && 2*block+int(contact)==mask-1))trial(j)=y(j);
            }
            if(mask==8) {
                trial=y;
                for(const auto& [i,z]:bc.psi)trial(i)=z;
                for(const auto& [i,z]:bc.phin)trial(N+i)=z-assembler.electronQuasiFermiReferenceAt(i)/potentialScale;
                for(const auto& [i,z]:bc.phip)trial(2*N+i)=z-assembler.holeQuasiFermiReferenceAt(i)/potentialScale;
            }
            for(int j=0;j<3*N;++j)changed+=trial(j)!=x(j);
            const VectorXd r=assembler.residual(trial,bc);
            log<<mask<<','<<changed<<','<<r.head(N).norm()<<','<<r.segment(N,N).norm()<<','<<r.tail(N).norm()<<','<<r.norm()<<','<<(r-original).norm()<<'\n';
        }
    }
'''

def build():
    a.verify(v.OUT/'completion_evidence_amended.json');LOCAL.mkdir(parents=True,exist_ok=False)
    src=v.REPO/'src/solver/NewtonSolver.cpp';text=src.read_text(encoding='utf-8')
    start=text.index('DDSolution NewtonSolver::makeSolution(');end=text.index('std::shared_ptr<CoupledDDAssembler> NewtonSolver::makeArclengthAssembler()',start)
    part=text[start:end];assert part.count('    return sol;')==1
    target=LOCAL/'NewtonSolver.cpp';target.write_text('#include <cstdlib>\n'+text[:start]+part.replace('    return sol;',HOOK+'\n    return sol;')+text[end:],encoding='utf-8')
    entry=next(e for e in a.read(v.REPO/'build-release/compile_commands.json') if Path(e['file'])==src)
    cmd=shlex.split(entry['command'].replace('\\','/'));cmd=[('-DVELA_VERSION="0.1.0"' if x.startswith('-DVELA_VERSION=') else x) for x in cmd]
    cmd[cmd.index('-c')+1]=str(target);cmd[cmd.index('-o')+1]=str(LOCAL/'newton.o');a.write(LOCAL/'compile.json',cmd)
    p=subprocess.run(cmd,cwd=entry['directory'],env=v.run.V.environment(),capture_output=True,text=True);(LOCAL/'compile.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr[-3000:]
    link=[cmd[0],str(LOCAL/'newton.o'),str(v.REPO/'build-release/CMakeFiles/vela_example_runner.dir/src/tools/vela_example_runner.cpp.obj'),str(v.REPO/'build-release/libvela_core.a')]
    link += ['D:/msys64/ucrt64/lib/lib'+n for n in ('spdlog.dll.a','fmt.a','umfpack.dll.a','spqr.dll.a','cholmod.dll.a')]+['-o',str(RUNNER)]
    a.write(LOCAL/'link.json',link);p=subprocess.run(link,env=v.run.V.environment(),capture_output=True,text=True);(LOCAL/'link.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
    d.matrix.freeze(OUT/'build_evidence.json',[Path(__file__).resolve(),src,v.OUT/'completion_evidence_amended.json']+[p for p in LOCAL.iterdir() if p.is_file()])

def replay():
    a.verify(OUT/'build_evidence.json')
    base=v.LOCAL/'dc/phumob/m65_n23_vd_1p000000_endpoint/native/vg_000/attempt_1'
    dest=LOCAL/'replay';cfg=a.read(base/'config.json');cfg['output_state_file']=str(dest/'state.csv');a.write(dest/'config.json',cfg)
    env=v.run.V.environment();env['VELA_PACK_AUDIT']=str(dest)
    status=v.run.V.execute(dest/'config.json',RUNNER,env);old=a.read(base/'config.status.json')
    same=a.sha(dest/'state.csv')==a.sha(base/'state.csv') and all(status[k]==old[k] for k in ('iterations','exit_code','final_residual','failure_reason'))
    a.write(OUT/'identity.json',dict(state_status_identical=same));assert same
    paths=sorted(dest.glob('state_*_norms.csv'));assert paths
    for p in paths:print(p.name,a.rows(p),flush=True)
    d.matrix.freeze(OUT/'replay_evidence.json',[OUT/'build_evidence.json',OUT/'identity.json',base/'config.json',base/'state.csv']+[p for p in dest.iterdir() if p.is_file()])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','replay'));globals()[p.parse_args().action]()
