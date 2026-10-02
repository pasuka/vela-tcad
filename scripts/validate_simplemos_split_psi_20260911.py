"""Isolate split-psi arithmetic at the six frozen rejected Newton directions.

No new DC acceptance or model qualification is inferred from this mixed-path probe.
"""
import argparse,csv,hashlib,json,os,shlex,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parents[1]
LOCAL=ROOT/'build-release/psi_split_0911'
OUT=ROOT/'reference_tcad/simplemos_sentaurus2022/phumob_split_psi_20260911'
OLD=ROOT/'reference_tcad/simplemos_sentaurus2022/phumob_packed_restart_20260910'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def rows(p):
    with Path(p).open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
def write(p,value):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',encoding='utf-8') as f:json.dump(value,f,indent=2)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        while b:=f.read(1024*1024):h.update(b)
    return h.hexdigest()
def freeze(p,files):write(p,dict(input_hashes={str(Path(f).relative_to(ROOT)).replace('\\','/'):sha(f) for f in files}))
def verify(p):
    for file,digest in read(p)['input_hashes'].items():
        assert sha(ROOT/file)==digest,file
def env():
    e=os.environ.copy();e['PATH']='D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+e.get('PATH','');return e
def execute(cfg,runner):
    p=subprocess.run([str(runner),'--config',str(cfg),'--log','off'],env=env(),capture_output=True,text=True,cwd=ROOT)
    cfg.with_suffix('.log').write_text(p.stdout+p.stderr,encoding='utf-8');return p.returncode
def csvout(p,data):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('x',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)

DUMP=r'''
        if(const char* path=std::getenv("VELA_SPLIT_INPUT")) {
            static bool written=false;
            if(!written) {
                written=true;nlohmann::json j;
                for(const auto& a:pn)j["nodes"].push_back({a.psi,a.n,a.p,a.xpsi,a.xn,a.xp,a.potentialScale,a.eref,a.href,a.ni,a.vt,a.doping,a.voln,a.volp,a.vold,a.q,a.area,a.scale,a.interfaceRhs});
                for(const auto& e:pe)j["edges"].push_back({e.i,e.j,e.g});
                j["bcpsi"]=nlohmann::json::array();
                for(const auto& [node,value]:bcs.psi)j["bcpsi"].push_back({node,value});
                std::ofstream(path)<<std::setprecision(17)<<j.dump(2);
            }
        }
'''

CPP=r'''
#include "vela/equation/ExtendedPoissonResidual.h"
#include "vela/numerics/SplitCoordinate.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iomanip>
int main(int argc,char**argv) {
    if(argc!=3)return 2;
    namespace pp=vela::poisson_precision;
    nlohmann::json j;std::ifstream(argv[1])>>j;
    std::vector<pp::Node> nodes;std::vector<pp::Edge> edges;
    for(const auto& a:j["nodes"])nodes.push_back({a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18]});
    for(const auto& e:j["edges"])edges.push_back({e[0],e[1],e[2]});
    const int N=nodes.size();std::ofstream out(argv[2]);out<<std::setprecision(40)<<"trial,mode,node,hi,lo,R\n";
    for(const auto& t:j["trials"]) {
        for(const std::string mode:{"rounded","split"}) {
            auto nn=nodes;std::vector<double> low(N,0.);
            for(int i=0;i<N;++i) {
                const auto pair=vela::numerics::SplitCoordinate::sum(t["x"][i],double(t["alpha"])*double(t["step"][i]));
                nn[i].xpsi=t["candidate"][i];nn[i].xn=t["candidate"][N+i];nn[i].xp=t["candidate"][2*N+i];
                if(mode=="split") {nn[i].xpsi=pair.hi;low[i]=pair.lo;}
            }
            auto r=pp::evaluate(nn,edges,true,low);
            for(const auto& bc:j["bcpsi"]) {int i=bc[0];r[i]=pp::Q(nn[i].xpsi)+pp::Q(low[i])-pp::Q(double(bc[1]));}
            for(int i=0;i<N;++i)out<<t["label"].get<std::string>()<<','<<mode<<','<<i<<','<<nn[i].xpsi<<','<<low[i]<<','<<r[i]<<'\n';
        }
    }
}
'''

def build():
    verify(OLD/'stalls/replay_evidence.json')
    dest=LOCAL/'probe';dest.mkdir(parents=True,exist_ok=True)
    assert not (OUT/'build_evidence.json').exists()
    src=ROOT/'src/equation/CoupledDDAssembler.cpp';s=src.read_text(encoding='utf-8')
    marker='        const auto pr=pp::evaluate(pn,pe,true);';assert s.count(marker)==1
    target=dest/'CoupledDDAssembler.cpp';target.write_text('#include <nlohmann/json.hpp>\n#include <cstdlib>\n#include <fstream>\n#include <iomanip>\n'+s.replace(marker,DUMP+'\n'+marker),encoding='utf-8')
    entry=next(e for e in read(ROOT/'build-release/compile_commands.json') if Path(e['file'])==src)
    cmd=shlex.split(entry['command'].replace('\\','/'));cmd=[('-DVELA_VERSION="0.1.0"' if x.startswith('-DVELA_VERSION=') else x) for x in cmd]
    cmd[cmd.index('-c')+1]=str(target);cmd[cmd.index('-o')+1]=str(dest/'assembler.o')
    if not (dest/'compile.json').exists():write(dest/'compile.json',cmd)
    p=subprocess.run(cmd,cwd=entry['directory'],env=env(),capture_output=True,text=True);(dest/'compile.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr[-4000:]
    obj=ROOT/'build-release/CMakeFiles/vela_example_runner.dir/src/tools/vela_example_runner.cpp.obj';lib=ROOT/'build-release/libvela_core.a'
    link=[cmd[0],str(dest/'assembler.o'),str(obj),str(lib)]+['D:/msys64/ucrt64/lib/lib'+n for n in ('spdlog.dll.a','fmt.a','umfpack.dll.a','spqr.dll.a','cholmod.dll.a')]+['-o',str(dest/'runner.exe')]
    write(dest/'link.json',link);p=subprocess.run(link,env=env(),capture_output=True,text=True);(dest/'link.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
    (dest/'kernel.cpp').write_text(CPP);cmd=[cmd[0],'-std=c++20','-O2','-I'+str(ROOT/'include'),str(dest/'kernel.cpp'),'-o',str(dest/'kernel.exe')]
    write(dest/'kernel_compile.json',cmd);p=subprocess.run(cmd,env=env(),capture_output=True,text=True);(dest/'kernel.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr
    freeze(OUT/'build_evidence.json',[Path(__file__).resolve(),src,obj,lib,ROOT/'build-release/vela_example_runner.exe',ROOT/'include/vela/numerics/SplitCoordinate.h',ROOT/'include/vela/equation/ExtendedPoissonResidual.h']+[p for p in dest.iterdir() if p.is_file()])

def run():
    verify(OUT/'build_evidence.json');verify(OLD/'stalls/replay_evidence.json');jobs=[];files=[OUT/'build_evidence.json',OLD/'stalls/replay_evidence.json']
    selected=rows(OLD/'cohort_v2/selected_states.csv')
    for c in rows(OLD/'stalls/identity.csv'):
        old=Path(c['dest']);dest=LOCAL/'states'/c['device']/c['vd']/c['index'];dest.mkdir(parents=True,exist_ok=False)
        cfg=read(old/'config.json');state=Path(cfg['state_file'])
        cfg.update(simulation_type='terminal_current_functional_probe',state_file=str(state),contact='drain',residual_output_csv=str(dest/'residual.csv'),contact_edge_output_csv=str(dest/'contact.csv'))
        cfg.pop('output_state_file',None);cfg['solver']['stable_merit_comparison']=False;cfg['solver'].pop('local_update_diagnostics',None)
        write(dest/'config.json',cfg)
        e=env();e['VELA_SPLIT_INPUT']=str(dest/'geometry.json')
        p=subprocess.run([str(LOCAL/'probe/runner.exe'),'--config',str(dest/'config.json'),'--log','off'],env=e,cwd=ROOT,capture_output=True,text=True);(dest/'probe.log').write_text(p.stdout+p.stderr);assert p.returncode==0,p.stderr[-1000:]
        data=read(dest/'geometry.json');base=rows(old/'trial_0.csv');n=len(base)//3
        actual=rows(dest/'residual.csv');names=('psi_residual','phin_residual','phip_residual')
        assert all(float(r['R'])==float(actual[i%n][names[i//n]]) for i,r in enumerate(base))
        for label,alpha,dr in [('base',0.,base)]+[(f'trial_{k}',2.**(-k),rows(old/f'trial_{k}.csv')) for k in range(13)]:
            data.setdefault('trials',[]).append(dict(label=label,alpha=alpha,x=[float(r['x']) for r in dr],step=[float(r['capped_step']) for r in dr],candidate=[float(r['x'] if label=='base' else r['candidate']) for r in dr]))
        write(dest/'inputs.json',data)
        p=subprocess.run([str(LOCAL/'probe/kernel.exe'),str(dest/'inputs.json'),str(dest/'kernel.csv')],env=env(),capture_output=True,text=True);assert p.returncode==0,p.stderr
        jobs.append(dict(c,new_dest=str(dest),coordinates=3*n,residual_identity=True));print(jobs[-1],flush=True)
        files += [state,old/'config.json']+[p for p in dest.iterdir() if p.is_file()]
    csvout(OUT/'identity.csv',jobs);freeze(OUT/'replay_evidence.json',files+[OUT/'identity.csv'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','run'));globals()[p.parse_args().action]()
