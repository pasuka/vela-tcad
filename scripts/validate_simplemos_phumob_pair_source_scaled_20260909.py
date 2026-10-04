"""Correct the diagnostic source's normalized-to-assembly unit conversion.

The first preflight's 4/4 source failures remain in the parent evidence root.
No failed experiment is promoted or used for a DC launch.
"""
import argparse
import subprocess
from pathlib import Path
import calibrate_simplemos_phumob_local_source_20260909 as s
import complete_simplemos_phumob_local_source_20260909 as j

BASE_LOCAL,BASE_OUT=s.LOCAL,s.OUT
s.LOCAL=BASE_LOCAL/'scaled_source';s.OUT=BASE_OUT/'scaled_source';s.RUNNER=s.LOCAL/'runner.exe'
j.ROOT=s.LOCAL/'jacobian_adapter';j.RUNNER=j.ROOT/'runner.exe'

def build():
    s.LOCAL.mkdir(parents=True,exist_ok=False);j.ROOT.mkdir(parents=True,exist_ok=False)
    text=(BASE_LOCAL/'CoupledDDAssembler.cpp').read_text()
    old='rate += alpha*source[node]/(vol_[node]*sf);'
    new='rate += alpha*source[node]*(scaling_.enabled ? scaling_.C0*scaling_.D0 : 1.)/(vol_[node]*sf);'
    assert text.count(old)==1
    src=s.LOCAL/'CoupledDDAssembler.cpp';src.write_text(text.replace(old,new))
    cmd=s.a.read(BASE_LOCAL/'compile_command.json');cmd=[('-DVELA_VERSION="0.1.0"' if x.startswith('-DVELA_VERSION=') else x) for x in cmd]
    cmd[cmd.index('-c')+1]=str(src);cmd[cmd.index('-o')+1]=str(s.LOCAL/'assembler.o')
    s.a.write(s.LOCAL/'compile_command.json',cmd)
    r=subprocess.run(cmd,cwd=s.REPO/'build-release',env=s.V.environment(),capture_output=True,text=True);(s.LOCAL/'compile.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    inputs=[Path(__file__).resolve(),Path(s.__file__),Path(j.__file__),BASE_OUT/'preflight_evidence.json',BASE_OUT/'build_evidence.json',BASE_OUT/'jacobian_build_evidence.json']
    for oldroot,newroot,exe in ((BASE_LOCAL,s.LOCAL,s.RUNNER),(BASE_LOCAL/'jacobian_adapter',j.ROOT,j.RUNNER)):
        link=s.a.read(oldroot/'link_command.json');link=[str(s.LOCAL/'assembler.o') if x==str(BASE_LOCAL/'assembler.o') else x for x in link];link[link.index('-o')+1]=str(exe)
        s.a.write(newroot/'link_command.json',link);r=subprocess.run(link,env=s.V.environment(),capture_output=True,text=True);(newroot/'link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
        inputs += [Path(x) for x in link[1:] if Path(x).is_file() and Path(x).is_relative_to(s.REPO)]
    s.a.write(s.OUT/'unit_amendment.json',dict(prior='Four zero/source preflights had relative insertion error 1; DC launch was blocked.',cause='Scope source is in normalized continuity-row units. nodeRecombinationRate feeds an unscaled assembly; continuity residual is subsequently divided by C0*D0. The initial injection omitted this factor.',repair='Multiply only the diagnostic additive source by C0*D0 before source assembly; preserve source definition, physical pair-source scope, alpha and all gates. No production source change.'))
    s.d.matrix.freeze(s.OUT/'build_evidence.json',inputs+[s.OUT/'unit_amendment.json']+[p for p in s.LOCAL.rglob('*') if p.is_file()])
    s.d.matrix.freeze(s.OUT/'jacobian_build_evidence.json',[s.OUT/'build_evidence.json',j.RUNNER])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','prepare','preflight','run','analyze'));x=p.parse_args().action
    build() if x=='build' else j.preflight() if x=='preflight' else getattr(s,x)()
