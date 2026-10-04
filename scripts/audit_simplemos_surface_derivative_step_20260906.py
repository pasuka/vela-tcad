"""Separately frozen step control for stable remaining surface-mobility defects."""
import argparse
from pathlib import Path
import subprocess
import audit_simplemos_jacobian_columns_20260906 as m

a=m.a;LOCAL=m.LOCAL/'surface_step';OUT=m.OUT/'surface_step';RUNNER=LOCAL/'runner.exe'


def build():
    a.verify(m.OUT/'freeze.json');LOCAL.mkdir(parents=True,exist_ok=False)
    source=(m.LOCAL/'CoupledDDAssembler.cpp').read_text()
    start=source.index('            if (transportMobilityDerivative || usesFermiDirac_) {')
    end=source.index('        // Isolated diagnostic: add only the missing field-drive chain term.',start)
    fragment=source[start:end];assert fragment.count('1.0e-6')==4
    fragment=fragment.replace('1.0e-6','(std::getenv("VELA_VALIDATE_TRANSPORT_STEP") ? std::stod(std::getenv("VELA_VALIDATE_TRANSPORT_STEP")) : 1.0e-6)')
    (LOCAL/'CoupledDDAssembler.cpp').write_text(source[:start]+fragment+source[end:],newline='\n')
    args=[str(LOCAL/'CoupledDDAssembler.cpp') if Path(s).name=='CoupledDDAssembler.cpp' else s for s in a.read(m.LOCAL/'build_command.json')]
    args[-1]=str(RUNNER);a.write(LOCAL/'build_command.json',args)
    p=subprocess.run(args,env=m.env(),capture_output=True,text=True);(LOCAL/'build.log').write_text(p.stdout+p.stderr)
    assert p.returncode==0,p.stderr[-2000:]
    print('Built isolated surface derivative step control',flush=True)


def prepare():
    files=[Path(__file__).resolve(),RUNNER,LOCAL/'CoupledDDAssembler.cpp',LOCAL/'build_command.json',m.OUT/'freeze.json']
    jobs=[]
    for c in a.read(m.OUT/'contract.json')['cases']:
        if c['variant']!='baseline_on' or '1p000000' not in c['case']:continue
        for step in (1e-7,1e-8):
            cfg=a.read(Path(c['config']));dest=LOCAL/c['case']/str(step)
            cfg['output_csv']=str(dest/'columns.csv');a.write(dest/'config.json',cfg);files.append(dest/'config.json')
            jobs.append(dict(case=c['case'],step=step,config=str(dest/'config.json')))
    a.write(OUT/'contract.json',dict(status='frozen_before_execution',jobs=jobs,
        trigger='Four stable psi-column carrier-block errors >.1% remain with vector-chain correction and disappear without Lombardi.',
        comparison_axis='Change only local transport/surface finite-difference step coefficient from 1e-6 to 1e-7 or 1e-8; vector chain step remains 1e-7. Physical residual unchanged.',
        gates=a.read(m.OUT/'contract.json')['gates'],production_changes=False))
    files.append(OUT/'contract.json');m.d.matrix.freeze(OUT/'freeze.json',files)


def run():
    a.verify(OUT/'freeze.json')
    for j in a.read(OUT/'contract.json')['jobs']:
        path=Path(j['config']);e=m.env();e['VELA_VALIDATE_TRANSPORT_STEP']=str(j['step'])
        p=subprocess.run([str(RUNNER),'--config',str(path),'--log','off'],env=e,capture_output=True,text=True)
        path.with_suffix('.stdout.txt').write_text(p.stdout);path.with_suffix('.stderr.txt').write_text(p.stderr)
        s=__import__('json').loads(p.stdout.strip().splitlines()[-1]);s['exit_code']=p.returncode;a.write(path.with_suffix('.status.json'),s)
        assert p.returncode==0,s
        print(j['case'],j['step'],'complete',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('build','prepare','run'))
    globals()[parser.parse_args().action]()
