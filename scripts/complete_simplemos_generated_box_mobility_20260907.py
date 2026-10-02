"""Retain the failed terminal-extraction attempt and rerun the same gates.

Only Newton-from-state terminal extraction now forwards the complete assembly
configuration. This script uses new output directories, seeds, and frozen inputs;
the original 32-job attempt remains immutable.
"""
import argparse
import shutil
from pathlib import Path
import validate_simplemos_generated_box_mobility_20260907 as s

INITIAL_LOCAL=s.LOCAL
INITIAL_OUT=s.OUT
s.LOCAL=INITIAL_LOCAL/'current_extraction_fix'
s.OUT=INITIAL_OUT/'current_extraction_fix'
s.RUNNER=s.LOCAL/'jvp_runner.exe'

def prepare():
    s.LOCAL.mkdir(parents=True,exist_ok=True)
    s.OUT.mkdir(parents=True,exist_ok=True)
    shutil.copy2(INITIAL_OUT/'prechange.json',s.OUT/'prechange.json')
    s.prepare()
    s.a.write(s.OUT/'addendum.json',dict(
        failure='First attempt: all 24 new-mode solves exported state but terminal extraction rejected missing assembly configuration; all original statuses retained.',
        change='Forward complete regionResolvedInterfaceAssembly to Newton-from-state ContactCurrent. No equation, Jacobian, solver gate, seed or model change.',
        rerun='Same 32 DC jobs and all fixed-state/cell/JVP checks in new directories.',
        initial_evidence=s.a.rel(INITIAL_OUT/'failed_attempt_evidence.json')))
    s.d.matrix.freeze(s.OUT/'addendum_freeze.json',[Path(__file__).resolve(),s.OUT/'freeze.json',s.OUT/'addendum.json'])

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('build','prepare','preflight','run','analyze'))
    action=parser.parse_args().action
    s.LOCAL.mkdir(parents=True,exist_ok=True)
    if action=='prepare':prepare()
    else:
        if action!='build':s.a.verify(s.OUT/'addendum_freeze.json')
        getattr(s,action)()
