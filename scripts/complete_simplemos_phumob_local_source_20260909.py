"""Add the missing read-only full-J export without changing the DC executable."""
import argparse
import inspect
import shlex
import subprocess
from pathlib import Path
import calibrate_simplemos_phumob_local_source_20260909 as s
import build_simplemos_constants_srh_calibration_20260907 as old

a,d=s.a,s.d
ROOT=s.LOCAL/'jacobian_adapter'
RUNNER=ROOT/'runner.exe'

def build():
    ROOT.mkdir(parents=True,exist_ok=False)
    src=s.REPO/'src/tools/vela_example_runner.cpp';text=src.read_text(encoding='utf-8')
    marker='int main(int argc, char** argv)';assert text.count(marker)==1
    text=text.replace(marker,old.JACOBIAN+'\n'+marker)
    marker='        } else if (type == "newton_residual_probe") {';assert text.count(marker)==1
    text=text.replace(marker,'        } else if (type == "parameter_jacobian") {\n            status.update(runParameterJacobian(configFile,cfg));\n'+marker)
    target=ROOT/'runner.cpp';target.write_text(text,encoding='utf-8')
    entry=next(x for x in a.read(s.REPO/'build-release/compile_commands.json') if x['file'].endswith('/vela_example_runner.cpp'))
    cmd=shlex.split(entry['command'].replace('\\','/'));cmd=[('-DVELA_VERSION="0.1.0"' if x.startswith('-DVELA_VERSION=') else x) for x in cmd]
    cmd[cmd.index('-c')+1]=str(target);cmd[cmd.index('-o')+1]=str(ROOT/'runner.o')
    a.write(ROOT/'compile_command.json',cmd)
    r=subprocess.run(cmd,cwd=entry['directory'],env=s.V.environment(),capture_output=True,text=True);(ROOT/'compile.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    link=a.read(s.LOCAL/'link_command.json');link=[str(ROOT/'runner.o') if 'vela_example_runner.cpp.obj' in x else x for x in link];link[link.index('-o')+1]=str(RUNNER)
    a.write(ROOT/'link_command.json',link);r=subprocess.run(link,env=s.V.environment(),capture_output=True,text=True);(ROOT/'link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    a.write(s.OUT/'jacobian_adapter.json',dict(reason='Production runner does not expose parameter_jacobian. Initial rejected probe retained; add a read-only runner adapter using the current solver API. DC executable and frozen solver unchanged.',source=str(src)))
    d.matrix.freeze(s.OUT/'jacobian_build_evidence.json',[Path(__file__).resolve(),src,Path(old.__file__),s.OUT/'jacobian_adapter.json',s.OUT/'build_evidence.json']+list(ROOT.glob('*')))

def preflight():
    a.verify(s.OUT/'jacobian_build_evidence.json');inputs=[]
    for c in a.read(s.OUT/'contract.json')['cases']:
        for label in ('zero','unit'):
            src=s.LOCAL/'fixed'/c['key']/label/'jacobian.json';cfg=a.read(src)
            cfg['output_csv']=str(src.with_name('jacobian_supported.csv'))
            dest=src.with_name('jacobian_supported.json');a.write(dest,cfg);inputs.append(dest)
    d.matrix.freeze(s.OUT/'jacobian_probe_freeze.json',inputs+[s.OUT/'freeze.json',s.OUT/'jacobian_build_evidence.json'])
    original=s.execute
    def execute(path,c,alpha):
        if path.stem=='jacobian_supported':return s.V.execute(path,RUNNER,s.env(c,alpha))
        return original(path,c,alpha)
    s.execute=execute
    text=inspect.getsource(s.preflight).replace("'jacobian'","'jacobian_supported'").replace("'jacobian.csv'","'jacobian_supported.csv'").replace("'jacobian.status.json'","'jacobian_supported.status.json'")
    exec(text,s.__dict__);s.preflight()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','preflight'));globals()[p.parse_args().action]()
