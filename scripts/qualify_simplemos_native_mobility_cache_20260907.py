"""Retain rejected candidate preflight; cover its separate analytic-J cache."""
import argparse,copy,subprocess
from pathlib import Path
import validate_simplemos_native_mobility_candidate_20260907 as f
import complete_simplemos_production_migration_20260907 as adapter

oldLOCAL=f.LOCAL;oldOUT=f.OUT;oldRUNNER=f.b.RUNNER
LOCAL=oldLOCAL/'cache_consistent';OUT=oldOUT/'cache_consistent';RUNNER=LOCAL/'runner.exe'
a=f.a;p=f.p

def build():
    a.verify(oldOUT/'preflight_evidence.json');LOCAL.mkdir(parents=True,exist_ok=False)
    source=(p.REPO/'src/equation/CoupledDDAssembler.cpp').read_text()
    marker='    return sum / static_cast<Real>(lowFieldMobilities.size());';assert source.count(marker)==1
    source=source.replace(marker,'    return (sum / static_cast<Real>(lowFieldMobilities.size())) * vela_candidate::factor(edgeId, mesh_.numEdges(), carrier == CarrierType::Electron);')
    (LOCAL/'CoupledDDAssembler.cpp').write_text(source,newline='\n')
    cmd=copy.deepcopy(a.read(oldLOCAL/'compile_commands.json')[0]);cmd[cmd.index('-o')+1]=str(LOCAL/'CoupledDDAssembler.o');cmd[cmd.index('-c')+1]=str(LOCAL/'CoupledDDAssembler.cpp')
    a.write(LOCAL/'compile_command.json',cmd);r=subprocess.run(cmd,cwd=p.REPO/'build-release',env=f.v.environment(),capture_output=True,text=True)
    (LOCAL/'build.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    cmd=a.read(oldLOCAL/'link_command.json');cmd=[str(LOCAL/'CoupledDDAssembler.o') if x==str(oldLOCAL/'CoupledDDAssembler.o') else x for x in cmd];cmd[-1]=str(RUNNER)
    a.write(LOCAL/'link_command.json',cmd);r=subprocess.run(cmd,env=f.v.environment(),capture_output=True,text=True)
    (LOCAL/'link.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    a.write(OUT/'cache_addendum.json',dict(status='fixed_before_new_preflight',retained_failure=str(oldOUT/'preflight.csv'),
        reason='Candidate-only edgeMobility override missed CoupledDDAssembler.cachedEdgeMobility used by analytic Jacobian. Zero-amplitude production blocks passed; 64 finite-candidate checks failed.',
        change='Multiply the cached fixed-mobility return by exactly the same carrier/edge factor. Residual, contact, thresholds and candidate direction unchanged.',production_source_changed=False))
    print('Cache-consistent candidate built',flush=True)

def route():
    f.LOCAL=LOCAL;f.OUT=OUT;f.b.RUNNER=RUNNER
    original_execute=adapter.execute
    def execute(path,runner=adapter.RUNNER,env=None):
        if path.is_relative_to(LOCAL):
            runner=RUNNER
            if env is not None:
                env=dict(env);env['VELA_CANDIDATE_MOBILITY_RATIOS']=str(LOCAL/'ratios'/(Path(env['VELA_CANDIDATE_MOBILITY_RATIOS']).name))
        return original_execute(path,runner,env)
    adapter.execute=execute

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('build','prepare','preflight','small','response','finite','analyze'));action=parser.parse_args().action
    if action=='build':build()
    else:
        route()
        if action=='prepare':
            f.prepare()
            f.d.matrix.freeze(OUT/'cache_addendum_freeze.json',[Path(__file__).resolve(),OUT/'cache_addendum.json',OUT/'freeze.json',LOCAL/'CoupledDDAssembler.cpp',LOCAL/'compile_command.json',LOCAL/'link_command.json',oldOUT/'preflight_evidence.json'])
        else:
            a.verify(OUT/'cache_addendum_freeze.json')
            if action in ('small','finite'):f.run(action=='finite')
            else:getattr(f,action)()
