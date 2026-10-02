"""Isolated multi-node SRH volume interpolation on coherent constant objects."""
import subprocess
from pathlib import Path
import validate_simplemos_constants_finite_20260908 as c

a=c.a;d=c.d;v=c.v;p=c.p;old=c.old
LOCAL=p.REPO/'build-release/simplemos_constants_srh_finite_20260908/srh'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh'
RUNNER=LOCAL/'runner.exe'
HOOK=r'''
    if(const char* path=std::getenv("VELA_CANDIDATE_SRH_VOLUMES")) {
        if(std::getenv("VELA_CANDIDATE_SRH_NODE"))
            throw std::runtime_error("Conflicting SRH diagnostic source selectors");
        const char* text=std::getenv("VELA_CANDIDATE_SRH_ALPHA");
        if(!text)throw std::runtime_error("Missing finite SRH amplitude");
        const double alpha=std::stod(text);
        if(!std::isfinite(alpha)||alpha<-.001||alpha>1.)
            throw std::runtime_error("SRH amplitude outside frozen interpolation range");
        std::ifstream in(path);std::size_t size,count;
        if(!(in>>size>>count)||size!=mesh_.numNodes()||count==0||count>2)
            throw std::runtime_error("Invalid finite SRH node manifest");
        std::set<Index> seen;
        for(std::size_t k=0;k<count;++k) {
            Index node;double ratio;
            if(!(in>>node>>ratio)||node>=size||contactNodes_[node]||ni_[node]<=0.||!seen.insert(node).second||!std::isfinite(ratio)||ratio<=0.)
                throw std::runtime_error("Invalid finite SRH volume record");
            const double factor=1.+alpha*(ratio-1.);
            if(!std::isfinite(factor)||factor<=0.)throw std::runtime_error("Nonpositive finite SRH volume");
            // Poisson volume copies were already constructed. In the frozen
            // SRH-only profile, vol_ is used for recombination and its Jacobian.
            vol_[node]*=factor;
        }
        std::string extra;if(in>>extra)throw std::runtime_error("Extra finite SRH data");
    }
'''

def build():
    a.verify(old.OUT/'build_evidence.json')
    LOCAL.mkdir(parents=True,exist_ok=False);OUT.mkdir(parents=True,exist_ok=False)
    source=(old.LOCAL/'CoupledDDAssembler.cpp').read_text(encoding='utf-8')
    assert source.count(old.b.SRH)==1
    source='#include <fstream>\n#include <set>\n'+source.replace(old.b.SRH,HOOK)
    target=LOCAL/'CoupledDDAssembler.cpp';target.write_text(source,encoding='utf-8')
    commands=a.read(old.LOCAL/'compile_commands.json');cmd=list(next(x for x in commands if Path(x[x.index('-c')+1]).name=='CoupledDDAssembler.cpp'))
    original_object=cmd[cmd.index('-o')+1];cmd[cmd.index('-c')+1]=str(target);cmd[cmd.index('-o')+1]=str(LOCAL/'CoupledDDAssembler.o')
    a.write(LOCAL/'compile_command.json',cmd);r=subprocess.run(cmd,cwd=p.REPO/'build-release',env=v.environment(),capture_output=True,text=True)
    (LOCAL/'build.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    link=[str(LOCAL/'CoupledDDAssembler.o') if x==original_object else x for x in a.read(old.LOCAL/'link_command.json')];link[link.index('-o')+1]=str(RUNNER)
    a.write(LOCAL/'link_command.json',link);r=subprocess.run(link,env=v.environment(),capture_output=True,text=True)
    (LOCAL/'link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    d.matrix.freeze(OUT/'build_evidence.json',[Path(__file__).resolve(),old.OUT/'build_evidence.json']+[x for x in LOCAL.rglob('*') if x.is_file()]+[Path(x) for x in link[1:] if x.endswith('.o')])
    print('Built finite SRH overlay; coherent constants objects retained',flush=True)

if __name__=='__main__':build()
