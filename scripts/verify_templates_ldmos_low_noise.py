"""Final evidence audit for the diagnostics-disabled D5 repeat experiment."""
import argparse
from pathlib import Path
from run_templates_ldmos_linked_d5 import read,write,digest
from analyze_templates_ldmos_isothermal_backends import summarize

def verify(root):
    summary=read(root/'summary.json');analysis=summarize(root)
    assert summary['status']=='pass' and summary['points']==31 and summary['rounds']==2
    assert summary['profiles']==['D5'] and summary['factor_statistics']=='off'
    assert len(analysis['cases'])==8 and len(analysis['joint_qualifications'])==4
    for gate in analysis['joint_qualifications'].values():
        assert gate['engineering']['status']==gate['final']['status']=='pass'
    assert all(x['exact_state_hashes_identical'] and x['same_target_and_update_history'] for x in analysis['repeats'])
    for gate in [4,8]:
        for backend in ['sparselu','umfpack']:
            cases=[c for c in analysis['cases'] if f'vg{gate}_' in c['name'] and c['name'].endswith(backend)]
            assert len(cases)==2 and cases[0]['counters']==cases[1]['counters']
    assert sum(c['audit']['rollbacks'] for c in analysis['cases'])==0
    previous=root.parent/'templates_ldmos_isothermal_umfpack_full_20260918'
    comparisons=[]
    for case in analysis['cases']:
        assert case['stages_seconds']['linear.factor_statistics']==0.
        assert not case['factor_statistics']['observations']
        ledger=read(root/case['name']/'fixed/ledger.json')
        old=read(previous/case['name']/'fixed/ledger.json')
        same=[(x['bias_V'],x['sha256']) for x in ledger['exact_points']]==[(x['bias_V'],x['sha256']) for x in old['exact_points']]
        assert same,case['name']+' differs from prior diagnostics-enabled state'
        assert [(r['target_V'],r['Newton_updates']) for r in ledger['runs']]==[(r['target_V'],r['Newton_updates']) for r in old['runs']]
        comparisons.append(dict(case=case['name'],matches_prior_diagnostics_enabled_states_and_history=True))
    checked={}
    for backend in ['sparselu','umfpack']:
        manifest=read(root/'binary'/f'{backend}.json')
        assert manifest['factor_statistics']=='off'
        for group in ['frozen_sources','runtime_sha256']:
            for path,sha in manifest[group].items():assert digest(Path(path))==sha;checked[path]=sha
        assert digest(Path(manifest['runner']))==manifest['runner_sha256']
    for path,sha in read(root/'matrix_inputs.json').items():assert digest(Path(path))==sha;checked[path]=sha
    result=dict(status='pass',exact_output_points=248,full_curves=8,joint_qualifications=4,
        newton_updates=sum(c['audit']['Newton_updates'] for c in analysis['cases']),
        rollbacks=0,repeated_linear_and_assembly_counters_identical=True,
        diagnostics_off_confirmed=True,prior_curve_comparisons=comparisons,
        hash_checks=len(checked),runner_sha256=summary['runner_sha256'],
        max_cross_backend_potential_difference_V=max(x['max_potential_difference_V'] for x in analysis['pairs']),
        max_cross_backend_current_relative_difference=max(x['max_current_relative_difference'] for x in analysis['pairs']))
    write(root/'analysis.json',analysis);write(root/'verification.json',result);return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);a=p.parse_args();print(verify(a.root))
