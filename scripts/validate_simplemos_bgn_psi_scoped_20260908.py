"""Repeat the unchanged validation matrix after limiting the repair's scope.

The original runner functions receive an explicit, separate experiment root
in this process. V1 inputs, failed regressions and outputs remain immutable.
"""
import argparse
from pathlib import Path
import validate_simplemos_bgn_psi_derivative_20260908 as run

ORIGINAL_LOCAL,ORIGINAL_OUT=run.LOCAL,run.OUT
LOCAL,OUT=ORIGINAL_LOCAL/'state_independent',ORIGINAL_OUT/'state_independent'
run.LOCAL,run.OUT=LOCAL,OUT
a,d=run.a,run.d


def prepare():
    a.verify(ORIGINAL_OUT/'candidate_v1_evidence.json')
    previous=a.read(ORIGINAL_OUT/'validation_contract.json')
    identities=a.read(ORIGINAL_OUT/'validation_freeze.json')['input_hashes']
    files=[]
    for job in previous['jobs']:
        cfg_path=Path(job['original_config'])
        cfg=a.read(cfg_path)
        paths=[cfg_path,Path(job['seed'])]+[Path(cfg[k]) for k in ('mesh_file','materials_file','node_doping_file')]
        for path in paths:
            key=path.resolve().relative_to(run.w.REPO).as_posix()
            assert a.sha(path)==identities[key],path
            files.append(path)
    contract=dict(previous)
    contract.update(change='Stable BGN psi partials only for state-independent mobility; retain field/surface approximate Jacobians when chain derivatives were deliberately disabled.',
        previous_candidate=str(ORIGINAL_OUT/'candidate_v1_evidence.json'),
        regression_repair='DCSweep test directories are atomically reserved to prevent cross-process path collisions. No numerical acceptance or frozen transition expectations changed.',
        input_matrix_changed=False)
    a.write(OUT/'validation_contract.json',contract)
    files += [Path(__file__).resolve(),Path(run.__file__),run.RUNNER,run.w.REPO/'build-release/libvela_core.a',
              ORIGINAL_OUT/'candidate_v1_evidence.json',OUT/'validation_contract.json']
    files += [run.w.REPO/name for name in ('src/equation/CoupledDDAssembler.cpp','include/vela/discretization/StableSGDerivative.h',
                                          'tests/test_production_numerics.cpp','tests/test_sg_flux.cpp','tests/test_dc_sweep.cpp','docs/config_schema.md')]
    d.matrix.freeze(OUT/'validation_freeze.json',files)
    print('Scoped repair frozen; all original seeds, models and gates unchanged.',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','solve','jvp','analyze'))
    action=parser.parse_args().action
    prepare() if action=='prepare' else getattr(run,action)()
