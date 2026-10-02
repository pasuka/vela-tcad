"""Path-type adapter for the frozen precision analyzer; equations unchanged."""
import argparse
from pathlib import Path
import analyze_simplemos_linesearch_precision_20260908 as core

original_out=core.OUT
core.OUT=original_out/'path_compatibility'
original_rows=core.a.rows
core.a.rows=lambda path: original_rows(Path(path))


def prepare():
    core.a.verify(original_out/'freeze.json')
    core.prepare()
    core.a.write(core.OUT/'adapter.json',dict(reason='Original analyzer stopped before calculations because the shared CSV helper requires pathlib.Path; normalize only path argument types.',numerical_contract_changed=False,original_freeze=core.a.rel(original_out/'freeze.json')))
    core.d.matrix.freeze(core.OUT/'runtime_freeze.json',[Path(__file__).resolve(),original_out/'freeze.json',core.OUT/'freeze.json',core.OUT/'adapter.json'])


def run():
    core.a.verify(core.OUT/'runtime_freeze.json');core.run()
    core.d.matrix.freeze(core.OUT/'runtime_evidence.json',[core.OUT/'runtime_freeze.json',core.OUT/'evidence.json'])


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run'));globals()[parser.parse_args().action]()
