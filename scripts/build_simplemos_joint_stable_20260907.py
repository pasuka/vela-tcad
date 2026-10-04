"""Separate numeric prerequisite after the retained baseline hole/psi JVP failure."""
from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor
import build_simplemos_joint_geometry_20260907 as old

p=old.p;a=p.a;d=p.d
LOCAL=p.REPO/'build-release/simplemos_joint_stable_20260907'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/joint_stable_20260907'
RUNNER=LOCAL/'runner.exe';HEADER=p.REPO/'scripts/diagnostics/simplemos_stable_sg_psi_derivative.hpp'
env=old.env

def main():
    a.verify(old.OUT/'freeze.json');LOCAL.mkdir(parents=True,exist_ok=False)
    source=(old.LOCAL/'CoupledDDAssembler.cpp').read_text()
    for carrier,qf,offset in [('Electron','phin','electron'),('Hole','phip','hole')]:
        marker=f'                add({qf}Offset() + i, psiOffset() + i, dF_dpsi_i);'
        assert source.count(marker)==1
        psi0='psi_i - electronQuantumPotential_V_(i)' if carrier=='Electron' else 'psi_i'
        psi1='psi_j - electronQuantumPotential_V_(j)' if carrier=='Electron' else 'psi_j'
        patch=f'''                if (compensatedEqualNiFlux_ && !bgnEnabled_ && niI == niJ) {{
                    const Real stableFlux = sg{carrier}BoltzmannContinuityFlux(
                        niI, niJ,
                        {psi0} - {offset}QuasiFermiReferenceAt(idxI),
                        {psi1} - {offset}QuasiFermiReferenceAt(idxI),
                        {qf}_i, {qf}_j_from_i, Vt_, coef,
                        SGBoltzmannFluxPolicy{{false,true}});
                    const auto stableDerivative = simplemos_stable_sg::psiDerivative(
                        stableFlux, {'Bminus, dBminusArg' if carrier=='Electron' else 'Bplus, dBplus'}, Vt_, {'true' if carrier=='Electron' else 'false'});
                    dF_dpsi_i=stableDerivative[0];dF_dpsi_j=stableDerivative[1];
                }}
'''
        source=source.replace(marker,patch+marker)
    (LOCAL/'CoupledDDAssembler.cpp').write_text('#include "simplemos_stable_sg_psi_derivative.hpp"\n'+source,newline='\n')
    original=a.read(old.LOCAL/'build_command.json');sources=[x for x in original if x.endswith('.cpp')]
    sources=[str(LOCAL/'CoupledDDAssembler.cpp') if Path(x).name=='CoupledDDAssembler.cpp' else x for x in sources]
    flags=[x for x in original[1:] if x.startswith(('-I','-D','-O','-std='))]
    jobs=[]
    for src in sources:
        obj=LOCAL/(Path(src).stem+'.o');jobs.append([original[0]]+flags+['-c',src,'-o',str(obj)])
    def compile(command):
        r=subprocess.run(command,env=old.env(LOCAL),capture_output=True,text=True)
        (LOCAL/(Path(command[-1]).stem+'.build.log')).write_text(r.stdout+r.stderr)
        assert r.returncode==0,r.stderr[-3000:];print('Compiled',Path(command[-1]).name,flush=True)
    a.write(LOCAL/'compile_commands.json',jobs)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(compile,jobs))
    libraries=[x for x in original if x.endswith('.a')]
    args=[original[0]]+[x[-1] for x in jobs]+libraries+['-o',str(RUNNER)];a.write(LOCAL/'build_command.json',args)
    r=subprocess.run(args,env=old.env(LOCAL),capture_output=True,text=True);(LOCAL/'build.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    print('Built separate stable SG psi-derivative runner',flush=True)

if __name__=='__main__':main()
