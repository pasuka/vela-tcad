"""Build an isolated opt-in mapping runner against the frozen VM core.

Run from a fresh campaign with payload/; scientific artifacts remain there.
"""
from pathlib import Path
import datetime, json, shutil, sys

R=Path.cwd(); B=R.parent
sys.path.insert(0,str(B/'s21p1_parallel'))
import common as C
C.END=datetime.datetime.fromisoformat('2026-10-22T19:16:00+08:00').timestamp()
P=B/'s22p1_complete_r2/build'; D=R/'build'
S=dict(completed=False,phase='prepare',deadline_epoch=C.END)
def save(): C.put(R/'build_status.json',S)
def replace(s,a,b):
    assert s.count(a)==1,(a,s.count(a))
    return s.replace(a,b)
def compile_one(name,path,extra=()):
    d=D/name;d.mkdir(parents=True)
    obj=d/(path.name+'.o')
    a=C.get(B/'s20p1_units/contract.json')['compile_command'][:]
    a[a.index('-c')+1]=str(path);a[a.index('-o')+1]=str(obj)
    a[a.index('-MF')+1]=str(obj)+'.d';a[a.index('-MT')+1]=str(obj)
    a[1:1]=['-I'+str(D/'include')]
    a+=['-I'+str(path.parent),'-I'+str(C.SOURCE/'src/tools'),*extra]
    C.run(d,'compile',a,cwd=C.SOURCE/'build',cap=600)
    return obj
try:
    (R/'build.lock').mkdir();save()
    D.mkdir(); (D/'include/vela/solver').mkdir(parents=True)
    shutil.copy2(R/'payload/SlotPopulationNumerics.h',D/'include/vela/solver/SlotPopulationNumerics.h')
    shutil.copytree(P/'tools',D/'tools')
    p=D/'tools/SlotPopulationMemory.inl';s=p.read_text()
    s=replace(s,' const bool stableFiniteUpdate=opt.value("stable_finite_update",false);',
              ' const bool stableFiniteUpdate=opt.value("stable_finite_update",false);\n const bool diluteLogInverse=opt.value("dilute_log_inverse",false);')
    s=replace(s,' result["stable_finite_update"]=stableFiniteUpdate;',
              ' result["stable_finite_update"]=stableFiniteUpdate;\n result["dilute_log_inverse"]=diluteLogInverse;')
    s=replace(s,'  int stableFiniteCount=0;double finiteUpdateChange=0.;',
              '  int stableFiniteCount=0;double finiteUpdateChange=0.;\n  int diluteLogSelected=0,diluteLogApplied=0;')
    s=replace(s,'    const double de=slot_population::relativeInverse(density[b][i]/dos[b][i],rel);',
              '    const bool logSelected=diluteLogInverse && slot_population::usesDiluteLogInverse(density[b][i]/dos[b][i],rel);\n'
              '    if(logSelected){++diluteLogSelected;++diluteLogApplied;}\n'
              '    const double de=diluteLogInverse ? slot_population::relativeInverseDiluteLog(density[b][i]/dos[b][i],rel) : slot_population::relativeInverse(density[b][i]/dos[b][i],rel);')
    s=replace(s,'      ++stableFiniteCount;','      ++stableFiniteCount;\n      if(logSelected)--diluteLogApplied;')
    token='  h["finite_update_arithmetic"]={{"stable_nodes",stableFiniteCount},{"max_qf_change_from_legacy_V",finiteUpdateChange},{"max_psi_writeback_error_V",psiWritebackError},{"consistent_psi_writeback",consistentPsiWriteback}};'
    s=replace(s,token,token+'\n  h["finite_update_arithmetic"]["dilute_log_selected"]=diluteLogSelected;\n  h["finite_update_arithmetic"]["dilute_log_applied"]=diluteLogApplied;')
    # P1 hypothetical finite maps must match the selected production map too.
    for expr in ['slot_population::relativeInverse(density[b][i]/dos[b][i],rel):0.',
                 'counterfactual=slot_population::relativeInverse(density[b][i]/dos[b][i],o.relativeBeforeHold[i]);']:
        if expr.startswith('counterfactual'):
            new='counterfactual=diluteLogInverse ? slot_population::relativeInverseDiluteLog(density[b][i]/dos[b][i],o.relativeBeforeHold[i]) : slot_population::relativeInverse(density[b][i]/dos[b][i],o.relativeBeforeHold[i]);'
        else:
            new='(diluteLogInverse ? slot_population::relativeInverseDiluteLog(density[b][i]/dos[b][i],rel) : slot_population::relativeInverse(density[b][i]/dos[b][i],rel)):0.'
        s=replace(s,expr,new)
    p.write_text(s)
    S['phase']='regression';save()
    objects=[compile_one('test_new',R/'payload/test_slot_log_inverse.cpp'),
             compile_one('test_statistics',R/'payload/test_carrier_statistics.cpp')]
    C.link(D,P/'libvela_core.a',D/'statistics',test=True,objects=objects)
    C.run(D,'statistics',[D/'statistics'],cap=240)
    S['phase']='compile_runner';save()
    obj=compile_one('runner_compile',D/'tools/vela_example_runner.cpp')
    C.link(D,P/'libvela_core.a',D/'runner',objects=[obj])
    C.put(R/'manifest.json',{str(p):C.sha(p) for p in [
        D/'runner',P/'runner',P/'libvela_core.a',P/'tools/SlotPopulationMemory.inl',
        D/'tools/SlotPopulationMemory.inl',D/'include/vela/solver/SlotPopulationNumerics.h',
        R/'payload/test_slot_log_inverse.cpp',R/'payload/test_carrier_statistics.cpp']})
    S.update(completed=True,phase='ready',passed=True);save()
except Exception as e:
    S.update(completed=True,passed=False,error=repr(e));save();raise
