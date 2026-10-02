"""Replay the qualified pair source with production Newton and Poisson numerics."""
import argparse
import shlex
import subprocess
from pathlib import Path
import validate_simplemos_combined_source_20260910 as prior
import validate_simplemos_production_unit_block_20260910 as production
import build_simplemos_constants_srh_calibration_20260907 as adapter

s=prior.s
LOCAL=production.LOCAL/'source';OUT=production.OUT/'source'
RUNNER=LOCAL/'runner.exe';JRUNNER=LOCAL/'jacobian_runner.exe'


def build():
    s.a.verify(production.OUT/'validation_freeze.json')
    LOCAL.mkdir(parents=True,exist_ok=False)
    source=s.REPO/'src/equation/CoupledDDAssembler.cpp';text=source.read_text(encoding='utf-8')
    hook=s.HOOK.replace('rate += alpha*source[node]/(vol_[node]*sf);',
        'rate += alpha*source[node]*(scaling_.enabled ? scaling_.C0*scaling_.D0 : 1.)/(vol_[node]*sf);')
    marker='    return rate;\n}\n\nVectorXd CoupledDDAssembler::electronDensity('
    assert text.count(marker)==1
    target=LOCAL/'CoupledDDAssembler.cpp';target.write_text('#include <fstream>\n#include <set>\n#include <cstdlib>\n'+text.replace(marker,hook+marker),encoding='utf-8')
    src=s.REPO/'src/tools/vela_example_runner.cpp';text=src.read_text(encoding='utf-8')
    marker='int main(int argc, char** argv)';assert text.count(marker)==1;text=text.replace(marker,adapter.JACOBIAN+'\n'+marker)
    marker='        } else if (type == "newton_residual_probe") {';assert text.count(marker)==1
    text=text.replace(marker,'        } else if (type == "parameter_jacobian") {\n            status.update(runParameterJacobian(configFile,cfg));\n'+marker)
    (LOCAL/'runner.cpp').write_text(text,encoding='utf-8')
    for original,target,obj in ((source,target,'assembler.o'),(src,LOCAL/'runner.cpp','runner.o')):
        entry=next(e for e in s.a.read(s.REPO/'build-release/compile_commands.json') if Path(e['file'])==original)
        cmd=shlex.split(entry['command'].replace('\\','/'))
        cmd=[('-DVELA_VERSION="0.1.0"' if x.startswith('-DVELA_VERSION=') else x) for x in cmd]
        cmd[cmd.index('-c')+1]=str(target);cmd[cmd.index('-o')+1]=str(LOCAL/obj)
        s.a.write(LOCAL/(obj+'_compile.json'),cmd)
        r=subprocess.run(cmd,cwd=entry['directory'],env=s.V.environment(),capture_output=True,text=True)
        (LOCAL/(obj+'_compile.log')).write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    for obj,exe in ((s.REPO/'build-release/CMakeFiles/vela_example_runner.dir/src/tools/vela_example_runner.cpp.obj',RUNNER),(LOCAL/'runner.o',JRUNNER)):
        link=[cmd[0],str(LOCAL/'assembler.o'),str(obj),str(s.REPO/'build-release/libvela_core.a')]
        link += ['D:/msys64/ucrt64/lib/lib'+name for name in ('spdlog.dll.a','fmt.a','umfpack.dll.a','spqr.dll.a','cholmod.dll.a')]
        link += ['-o',str(exe)];s.a.write(LOCAL/(exe.stem+'_link.json'),link)
        r=subprocess.run(link,env=s.V.environment(),capture_output=True,text=True)
        (LOCAL/(exe.stem+'_link.log')).write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    s.d.matrix.freeze(OUT/'build_evidence.json',[Path(__file__).resolve(),Path(adapter.__file__),Path(s.__file__),production.OUT/'validation_freeze.json',source,src]+[p for p in LOCAL.iterdir() if p.is_file()])


def bind():
    def execute(path,case,alpha):return s.V.execute(path,JRUNNER if path.stem=='jacobian' else RUNNER,s.env(case,alpha))
    s.LOCAL=LOCAL;s.OUT=OUT;s.RUNNER=RUNNER;s.execute=execute


def prepare():
    s.a.verify(OUT/'build_evidence.json');s.a.verify(prior.OUT/'combined_evidence.json')
    contract=s.a.read(prior.OUT/'contract.json');files=[OUT/'build_evidence.json',prior.OUT/'combined_evidence.json']
    for case in contract['cases']:
        base=Path(case['baseline']);files += [base/'state.csv',base/'all_row.csv',Path(case['source_file'])]
        for job in case['jobs']:
            cfg=s.a.read(Path(job['dest'])/'config.json');dest=LOCAL/'dc'/case['key']/job['label']
            cfg['solver'].update(production.OPTIONS);cfg['output_state_file']=str(dest/'state.csv')
            s.a.write(dest/'config.json',cfg);s.V.post_config(cfg,dest);job['dest']=str(dest);files+=list(dest.glob('*.json'))
        cfg=s.a.read(Path(case['jobs'][0]['dest'])/'config.json')
        for label in ('zero','unit'):
            dest=LOCAL/'fixed'/case['key']/label;files+=s.V.probes(cfg,dest,base/'state.csv')
            deck=s.a.read(dest/'functional.json');deck.update(simulation_type='parameter_jacobian',output_csv=str(dest/'jacobian.csv'))
            s.a.write(dest/'jacobian.json',deck);files.append(dest/'jacobian.json')
    contract['scope']='Same four baselines and twenty pair-source perturbations, current production numerics with explicit options; diagnostic overlay adds only the unchanged fixed source and read-only full-J export.'
    s.a.write(OUT/'contract.json',contract);s.a.write_csv(OUT/'scope.csv',s.a.rows(prior.OUT/'scope.csv'))
    s.d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json',OUT/'scope.csv'])
    bind();s.preflight()


def run():bind();s.run()


def analyze():
    bind();s.analyze()
    dc=s.a.rows(OUT/'dc.csv');native=s.a.rows(prior.ORIGINAL_NATIVE/'terminal_derivatives.csv');ports=[];identity=[]
    for r in s.a.rows(OUT/'ports.csv'):
        if r['port'] not in ('drain','substrate'):continue
        nr=next(n for n in native if all(n[k]==r[k] for k in ('key','amplitude','port')))
        error=abs(float(r['derivative_A_per_um'])/float(nr['derivative_A_per_um'])-1)
        ports.append(dict(key=r['key'],amplitude=r['amplitude'],port=r['port'],relative=error,qualified=error<=1e-3 and all(x['qualified']=='True' for x in dc if x['key']==r['key'])))
    for case in s.a.read(OUT/'contract.json')['cases']:
        for job in case['jobs']:
            old=prior.LOCAL/'dc'/case['key']/job['label'];new=Path(job['dest'])
            astate=(old/'state.csv').read_bytes();bstate=(new/'state.csv').read_bytes()
            oid=s.a.read(old/'config.status.json')['contact_currents_A_per_um'];nid=s.a.read(new/'config.status.json')['contact_currents_A_per_um']
            identity.append(dict(key=case['key'],label=job['label'],state_bitwise_identical=astate==bstate,
                max_port_relative=max(abs(nid[k]/oid[k]-1) for k in oid if oid[k]!=0),ports_identical=oid==nid))
    s.a.write_csv(OUT/'native_ports.csv',ports);s.a.write_csv(OUT/'diagnostic_identity.csv',identity)
    summary=s.a.read(OUT/'summary.json');summary.update(qualified_native_ports=sum(r['qualified'] for r in ports),native_ports=len(ports),
        max_native_relative=max(r['relative'] for r in ports),identical_states=sum(r['state_bitwise_identical'] for r in identity),
        identical_ports=sum(r['ports_identical'] for r in identity),production_numerics=True)
    s.a.write(OUT/'production_summary.json',summary)
    s.d.matrix.freeze(OUT/'production_evidence.json',[OUT/'evidence.json',prior.OUT/'combined_evidence.json',OUT/'production_summary.json',OUT/'native_ports.csv',OUT/'diagnostic_identity.csv'])
    print(summary,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','prepare','run','analyze'));globals()[p.parse_args().action]()
