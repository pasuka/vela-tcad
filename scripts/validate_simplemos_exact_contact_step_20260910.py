"""Isolate exact identity-row Newton updates from linear solver roundoff."""
import argparse
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import validate_simplemos_poisson_precision_20260910 as c
s=c.s;LOCAL=c.LOCAL/'contact_step';OUT=c.OUT/'contact_step';RUNNER=LOCAL/'runner.exe'
HOOK=r'''
        if(const char* exactText=std::getenv("VELA_EXACT_CONTACT_STEP")) {
            if(std::string(exactText)!="1")throw std::runtime_error("Invalid exact contact step selector");
            const auto projectIdentityRow=[&](int row) {
                if(const char* path=std::getenv("VELA_EXACT_CONTACT_TRACE")) {
                    const bool header=!std::filesystem::exists(path);std::ofstream out(path,std::ios::app);
                    if(header)out<<"iteration,row,original_step,identity_step,x,R\n";
                    out<<std::setprecision(17)<<iter<<','<<row<<','<<step(row)<<','<<-r(row)<<','<<x(row)<<','<<r(row)<<'\n';
                }
                step(row)=-r(row);
            };
            for(const auto& [node,value]:bcs.psi)projectIdentityRow(static_cast<int>(node));
            for(const auto& [node,value]:bcs.phin)projectIdentityRow(N+static_cast<int>(node));
            for(const auto& [node,value]:bcs.phip)projectIdentityRow(2*N+static_cast<int>(node));
        }
'''

def build():
    s.a.verify(c.OUT/'build_evidence.json');LOCAL.mkdir(parents=True,exist_ok=False)
    text=(c.m.LOCAL/'NewtonSolver.cpp').read_text(encoding='utf-8');marker='        const VectorXd rawStep = step;';assert text.count(marker)==1
    (LOCAL/'NewtonSolver.cpp').write_text(text.replace(marker,HOOK+'\n'+marker),encoding='utf-8')
    cmd=s.a.read(c.m.LOCAL/'compile_command.json');cmd[cmd.index('-c')+1]=str(LOCAL/'NewtonSolver.cpp');cmd[cmd.index('-o')+1]=str(LOCAL/'newton.o');s.a.write(LOCAL/'compile_command.json',cmd)
    r=subprocess.run(cmd,cwd=s.REPO/'build-release',env=s.V.environment(),capture_output=True,text=True);(LOCAL/'compile.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    link=s.a.read(c.LOCAL/'link_command.json');assert link[1]==str(c.m.LOCAL/'newton.o');link[1]=str(LOCAL/'newton.o');link[link.index('-o')+1]=str(RUNNER);s.a.write(LOCAL/'link_command.json',link)
    r=subprocess.run(link,env=s.V.environment(),capture_output=True,text=True);(LOCAL/'link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    s.d.matrix.freeze(OUT/'build_evidence.json',[Path(__file__).resolve(),c.OUT/'build_evidence.json']+[p for p in LOCAL.iterdir() if p.is_file()])

def run():
    s.a.verify(OUT/'build_evidence.json');jobs=[];files=[OUT/'build_evidence.json',c.OUT/'freeze.json']
    for case in s.a.read(c.OUT/'contract.json')['jobs']:
        if case['mode']!='packed' or (case['device'],case['label']) not in (('n19','minus_full'),('n23','zero'),('n23','plus_full')):continue
        base=Path(case['dest']);dest=LOCAL/case['key']/case['label'];cfg=s.a.read(base/'config.json');cfg['output_state_file']=str(dest/'state.csv');s.a.write(dest/'config.json',cfg);s.V.post_config(cfg,dest)
        jobs.append(dict(case,base=str(base),dest=str(dest)));files+=list(dest.glob('*.json'))+[Path(cfg['state_file'])]
    s.a.write(OUT/'contract.json',dict(jobs=jobs,scope='After linear refinement, set the Newton step of each exact Dirichlet identity row to minus its residual. Leave all free-row updates, caps, merit, source, packed Poisson and 200-iteration/gating settings unchanged. Three isolated controls, not uniform 16-point qualification.'))
    s.d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json'])
    def one(j):
        dest=Path(j['dest']);env=s.env(j,j['alpha']);env.update(VELA_STABLE_MERIT='1',VELA_POISSON_PRECISION='packed',VELA_EXACT_CONTACT_STEP='1',VELA_EXACT_CONTACT_TRACE=str(dest/'contact.csv'),VELA_STABLE_MERIT_TRACE=str(dest/'merit.csv'))
        s.V.execute(dest/'config.json',RUNNER,env)
        for name in ('all_row','acceptance_edges'):s.V.execute(dest/(name+'.json'),RUNNER,env)
        row=dict(**j,**s.prior.q.run.w.old.prior.old.qualify(j,dest));print(row,flush=True);return row
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(one,jobs))
    s.a.write_csv(OUT/'dc.csv',rows);s.d.matrix.freeze(OUT/'dc_evidence.json',[OUT/'freeze.json',OUT/'dc.csv']+[p for j in jobs for p in Path(j['dest']).rglob('*') if p.is_file()])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','run'));globals()[p.parse_args().action]()
