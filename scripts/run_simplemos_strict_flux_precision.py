"""Single-axis stable equal-ni flux evaluation under unchanged strict gates."""
import argparse
import json
from pathlib import Path
import run_simplemos_convergence_acceptance_isolation as previous
import run_simplemos_convergence_audit as audit
import trace_simplemos_strict_rejection as trace

REPO,ROOT=previous.REPO,previous.ROOT
LOCAL=REPO/'build-release/simplemos_strict_flux_precision'
OUT=ROOT/'strict_flux_precision'
CONTRACT=OUT/'contract.json'
FREEZE=OUT/'freeze.json'
SCRIPT=Path(__file__).resolve()


def prepare():
    if CONTRACT.exists():raise FileExistsError(CONTRACT)
    audit.verify(previous.FREEZE)
    original=audit.read(previous.CONTRACT)
    paths=[SCRIPT,previous.CONTRACT,previous.FREEZE,previous.OUT/'case_ledger.csv',previous.first.PRODUCTION,
           REPO/'src/discretization/ScharfetterGummel.cpp',REPO/'src/equation/CoupledDDAssembler.cpp',REPO/'src/post/ContactCurrent.cpp',trace.CONTRACT]
    for c in original['cases']:
        src=previous.LOCAL/c['case'];dest=LOCAL/c['case']
        deck=audit.read(src/'config.json')
        deck['solver']['bandgap_narrowing']['equal_ni_flux_evaluation']='compensated_log_expm1'
        deck['output_state_file']=str(dest/'state.csv')
        # Only precision-evaluation option changes inside the solver.
        audit.write(dest/'config.json',deck)
        paths += [src/'config.json',Path(deck['state_file']),dest/'config.json']
    audit.write(CONTRACT,{'status':'frozen_before_execution','cases':original['cases'],'additional_nonlinear_reclosures':4,
        'single_axis':'solver.bandgap_narrowing.equal_ni_flux_evaluation=compensated_log_expm1, including model=none. Same analytic equations, ni, Vt, mobility, original states, row gates and external KCL.',
        'trigger':'Strict rejection trace reproduced prior states; failed step carrier caps absent and raw linear solve accurate. Existing legacy equal-ni path subtracts two nearly equal exp values, whereas alternate policy uses a log-ratio/expm1 evaluation.',
        'gates':original['gates'],'global_profile':original['global_profile'],
        'interpretation':'Precision experiment, not a physical current correction. Even 4/4 strict success requires separate assessment of Id and high/low paired errors. No automatic M82/M83 release.'})
    paths.append(CONTRACT)
    audit.write(FREEZE,{'input_hashes':{audit.rel(p):audit.sha(p) for p in sorted(set(paths))}})


def run():
    audit.verify(FREEZE)
    previous.LOCAL=LOCAL
    results=[]
    for c in audit.read(CONTRACT)['cases']:
        if (LOCAL/c['case']/'result.json').exists():raise FileExistsError(LOCAL/c['case'])
        results.append(previous.run_case(c))
    audit.write_csv(OUT/'case_ledger.csv',results)
    audit.write(OUT/'summary.json',{'cases':4,'qualified':sum(r['qualified'] for r in results),
        'new_nonlinear_solves':4,'m82_released':False,'production_changes':False})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','run','verify'))
    {'prepare':prepare,'run':run,'verify':lambda:audit.verify(FREEZE)}[p.parse_args().action]()
