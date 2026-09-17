"""Actual first8 pause/resume against an independently completed combined curve."""
import argparse
import json
import os
from pathlib import Path
import subprocess

from analyze_templates_ldmos_neutral_cache import signature
from run_templates_ldmos_contact_sweep import read, sha
from run_templates_ldmos_electrothermal_curve import state_gate


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--reference',type=Path,help='Explicit completed curve; default is matrix/r0_vg8_combined')
    a=p.parse_args();matrix=a.matrix.resolve();out=a.output.resolve()
    summary=read(matrix/'summary.json');assert summary['status']=='complete'
    runner=matrix/'vela_example_runner';assert sha(runner)==summary['runner_sha256']
    reference=a.reference.resolve() if a.reference else matrix/'r0_vg8_combined';old=read(reference/'results/ledger.json')
    assert old['status']=='complete' and len(old['exact_points'])==31
    out.mkdir(parents=True,exist_ok=False)
    cfg,deck=read(reference/'input.json'),read(reference/'deck.json')
    deck.update(input_file=str(out/'input.json'),output_directory=str(out/'results'),pause_after_attempts=3)
    deck['sweep']['bias_points_V']=deck['sweep']['bias_points_V'][:8]
    (out/'input.json').write_text(json.dumps(cfg));(out/'deck.json').write_text(json.dumps(deck))
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    def execute(name,expected):
        with (out/(name+'.log')).open('x') as log:
            code=subprocess.call([str(runner),'--config',str(out/(name+'.json'))],cwd=str(out),env=env,stdout=log,stderr=subprocess.STDOUT)
        assert code==expected,(name,code)
    execute('deck',1)
    ledger=read(out/'results/ledger.json');assert ledger['status']=='stopped_at_checkpoint'
    accepted=Path(ledger['accepted_result']);accepted_hash=sha(accepted)
    (out/'paused_ledger.json').write_text(json.dumps(ledger,indent=2))
    deck.pop('pause_after_attempts');deck['resume']=True
    (out/'resume.json').write_text(json.dumps(deck));execute('resume',0)
    assert sha(accepted)==accepted_hash
    ledger=read(out/'results/ledger.json');assert ledger['status']=='complete' and len(ledger['exact_points'])==8
    hits=[];counts={}
    for stage in ('initialization_runs','runs'):
        rows=ledger[stage];counts[stage]=len(rows)
        if stage=='initialization_runs':assert len(rows)==len(old[stage])
        for row,prior in zip(rows,old[stage]):
            path=Path(row['result']) if stage=='initialization_runs' else Path(row['directory'])/'output.json'
            ref=Path(prior['result']) if stage=='initialization_runs' else Path(prior['directory'])/'output.json'
            result=read(path);assert signature(result)==signature(read(ref))
            hits.append(result['performance']['static_preparation_reused'])
            if stage=='runs':
                for key in ('bias_V','parent_bias_V','prediction'):assert row[key]==prior[key]
    assert len(ledger['runs'])<=len(old['runs'])
    for point,prior in zip(ledger['exact_points'],old['exact_points']):
        assert point['bias_V']==prior['bias_V']
        assert state_gate(read(point['result']),point['bias_V'])['pass_gate']
    assert sum(hits)==len(hits)-2,hits
    report=dict(status='pass',runner_sha256=sha(runner),exact_to_uninterrupted=True,
                accepted_checkpoint_preserved=True,cache_hits=sum(hits),point_services=len(hits),counts=counts)
    (out/'audit.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))


if __name__=='__main__':main()
