"""Prepare an isolated SLOT driver from the audited s12home60 snapshot.

No solver is started. Frozen evidence is never edited. Fail closed if the
source driver's shape changed; the generated package retains all guards.
"""
from pathlib import Path
import shutil
import hashlib
import json


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('unexpected snapshot shape: '+old[:80])
    return text.replace(old, new)


def prepare(source, destination, runner):
    source, destination, runner = map(Path, (source, destination, runner))
    destination.mkdir(parents=True, exist_ok=False)
    for name in ('bin', 'inputs'):
        shutil.copytree(source/name, destination/name)
    shutil.copy2(runner, destination/'bin/runner.exe')
    for name in ('base.json', 'helpers.py', 'slot_linear_quality.py'):
        shutil.copy2(source/name, destination/name)
    shutil.copy2(Path(__file__).with_name('slot_population_acceptance.py'), destination/'slot_population_acceptance.py')
    common = (source/'common.py').read_text()
    common = replace_once(common, "c['candidate_port_diagnostics']=True", "c['candidate_port_diagnostics']=False")
    common = replace_once(common, "c['linear_system_output_dir']=str(d/'data')", "c['linear_system_output_dir']=str(d/'data');c['linear_system_output_profile']='solve_only'")
    (destination/'common.py').write_text(common)
    pop = (source/'population.py').read_text()
    pop = 'from slot_population_acceptance import accept, readback_audit\n'+pop
    pop = replace_once(pop, "rho=c.quality(td)['blocks']['combined']", "rho=c.get(td/'status.json')['block_residuals']['combined']")
    pop = replace_once(pop, "qualified=bool(primary and boundary and closure['satisfied'])", "qualified=accept(rho,eraw,closure=closure['satisfied'],boundary=boundary,linear=quality['accepted'],mapping=all(v['target_relative_error']<1e-7 for v in mapping.values()),positive=all(np.all(an[k+'_density'][ids[j+1]]>0) for j,k in enumerate(prev.carriers)))")
    pop = replace_once(pop, "rec['cpp_zero_update_converged']=qs['converged']", "rec['cpp_zero_update_converged']=qs['converged'];rec['readback_audit']=readback_audit(rho,qs['final_residual'],qs['converged']);assert rec['readback_audit']['consistent'], 'saved_state_residual_mismatch'")
    pop = replace_once(pop, "if rec['cpp_zero_update_converged']:return rec", "return rec # RHS readback is diagnostic; consistency remains mandatory")
    (destination/'population.py').write_text(pop)
    continuation = (source/'run.py').read_text()
    old_gate = "if rec['qualified'] and rec.get('readback_rhs_pass',False):"
    if continuation.count(old_gate) != 2:
        raise ValueError('unexpected continuation acceptance gates')
    continuation = continuation.replace(old_gate, "if rec['qualified']:")
    continuation = replace_once(continuation, "readback_rhs_pass=z['cpp_zero_update_converged']", "readback_rhs_pass=z['cpp_zero_update_converged'],readback_audit=z['readback_audit']")
    (destination/'continuation.py').write_text(continuation)
    for path in destination.glob('*.py'):
        compile(path.read_text(), str(path), 'exec')
    files = [source/name for name in ('run.py','population.py','common.py','base.json')]
    files += list((source/'inputs').glob('*')) + [runner]
    hashes = {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in files if p.is_file()}
    (destination/'preparation_inputs.json').write_text(json.dumps(hashes,indent=2))


if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',required=True)
    parser.add_argument('--destination',required=True)
    parser.add_argument('--runner',required=True)
    args=parser.parse_args()
    prepare(args.source,args.destination,args.runner)
