"""Isolated fixed edge-mobility interpolation on the current production core."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess
import validate_simplemos_production_migration_20260907 as v
import build_simplemos_joint_geometry_20260907 as old
a=v.a;p=v.p;LOCAL=p.REPO/'build-release/simplemos_native_mobility_candidate_20260907'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/native_mobility_candidate_20260907';RUNNER=LOCAL/'runner.exe'

HOOK=r'''
namespace vela_candidate {
inline double factor(std::size_t edge, std::size_t size, bool electron) {
    const char* value=std::getenv("VELA_CANDIDATE_MOBILITY_ALPHA");
    if (!value || std::stod(value)==0.) return 1.;
    const double alpha=std::stod(value);
    if (!std::isfinite(alpha)) throw std::runtime_error("Invalid mobility amplitude");
    static const auto ratios=[&] {
        const char* path=std::getenv("VELA_CANDIDATE_MOBILITY_RATIOS");
        if(!path)throw std::runtime_error("Missing mobility ratios");
        std::ifstream in(path);std::vector<std::array<double,2>> rows(size);
        for(std::size_t k=0;k<size;++k) {
            std::size_t id;
            if(!(in>>id>>rows[k][0]>>rows[k][1]) || id!=k ||
                !std::isfinite(rows[k][0]) || !std::isfinite(rows[k][1]) ||
                rows[k][0]<=0 || rows[k][1]<=0)throw std::runtime_error("Invalid mobility ratio row");
        }
        std::string extra;if(in>>extra)throw std::runtime_error("Extra mobility ratio rows");
        return rows;
    }();
    if(size!=ratios.size())throw std::runtime_error("Mobility ratio mesh mismatch");
    const double multiplier=1.+alpha*(ratios.at(edge)[electron?0:1]-1.);
    if(!(multiplier>0.) || !std::isfinite(multiplier))throw std::runtime_error("Invalid mobility multiplier");
    return multiplier;
}
}
'''

def main():
    LOCAL.mkdir(parents=True,exist_ok=False)
    header=(p.REPO/'include/vela/equation/AssemblerUtils.h').read_text();begin=header.index('inline Real edgeMobility(');end=header.index('/// Return average model mobility',begin)
    body=header[begin:end];marker='    return sum / static_cast<Real>(contributingCells);';assert body.count(marker)==1
    body=body.replace(marker,'    return (sum / static_cast<Real>(contributingCells)) * vela_candidate::factor(edgeId, mesh.edges().size(), carrier == CarrierType::Electron);')
    header=header[:begin]+body+header[end:];include=LOCAL/'include/vela/equation/AssemblerUtils.h';include.parent.mkdir(parents=True)
    include.write_text('#include <fstream>\n#include <cstdlib>\n#include <array>\n#include <vector>\n#include <cmath>\n#include <stdexcept>\n'+HOOK+header,newline='\n')
    runner=(p.REPO/'src/tools/vela_example_runner.cpp').read_text();marker='int main(int argc, char** argv)';assert runner.count(marker)==1
    runner=runner.replace(marker,old.FUNCTION+'\n'+marker)
    marker='        } else if (type == "newton_residual_probe") {';assert runner.count(marker)==1
    runner=runner.replace(marker,'        } else if (type == "joint_geometry_jvp") {\n            status.update(runJointGeometryJvp(configFile,cfg));\n'+marker)
    (LOCAL/'runner.cpp').write_text(runner,newline='\n')
    # Use current Release flags, including actual available HDF5/SuiteSparse defines.
    database=a.read(p.REPO/'build-release/compile_commands.json')
    import shlex
    jobs=[]
    for source in ('src/equation/CoupledDDAssembler.cpp','src/post/ContactCurrent.cpp','src/tools/vela_example_runner.cpp'):
        entry=next(x for x in database if Path(x['file']).as_posix().endswith(source))
        args=shlex.split(entry['command'],posix=False)
        # CMake on Windows quotes only paths; remove quotes without touching define text.
        args=[(x[1:-1] if x.startswith('"') and x.endswith('"') else x).replace('\\"','"') for x in args]
        oi=args.index('-o');args[oi+1]=str(LOCAL/(Path(source).stem+'.o'))
        ci=args.index('-c');args[ci+1]=str(LOCAL/'runner.cpp') if 'vela_example_runner' in source else str(p.REPO/source)
        args.insert(1,'-I'+str(LOCAL/'include'));jobs.append(args)
    a.write(LOCAL/'compile_commands.json',jobs)
    def compile(cmd):
        r=subprocess.run(cmd,cwd=p.REPO/'build-release',env=v.environment(),capture_output=True,text=True)
        (LOCAL/(Path(cmd[cmd.index('-o')+1]).stem+'.build.log')).write_text(r.stdout+r.stderr)
        assert r.returncode==0,r.stderr[-4000:];print('Candidate compiled',Path(cmd[cmd.index('-o')+1]).name,flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(compile,jobs))
    # Obtain the actual CMake linker command and replace only this executable's objects.
    command=subprocess.check_output(['D:/msys64/ucrt64/bin/ninja.exe','-t','commands','vela_example_runner'],cwd=p.REPO/'build-release',text=True).splitlines()[-1]
    args=command.replace('"','').split()
    # Ninja Windows command starts/ends with cmd /C wrappers. Libraries are explicit paths.
    libs=[x for x in args if x.endswith('.a') and not x.startswith('-')]
    libs=[str(p.REPO/'build-release'/x) if not Path(x).is_absolute() else x for x in libs]
    cmd=['D:/msys64/ucrt64/bin/g++.exe']+[j[j.index('-o')+1] for j in jobs]+libs+['-o',str(RUNNER)]
    a.write(LOCAL/'link_command.json',cmd);r=subprocess.run(cmd,env=v.environment(),capture_output=True,text=True)
    (LOCAL/'link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-4000:]
    print('Isolated native mobility candidate built',flush=True)

if __name__=='__main__':main()
