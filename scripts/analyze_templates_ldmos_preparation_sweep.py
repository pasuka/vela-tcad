"""Verify both full repeat rounds on original VM paths, including failed attempts."""
import argparse
import json
from pathlib import Path
import statistics

from analyze_templates_ldmos_neutral_cache import inspect
from run_templates_ldmos_contact_sweep import read,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('matrix','reference','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();summary=read(a.matrix/'summary.json')
    assert summary['status']=='complete' and summary['repeats']>=2
    assert sha(a.matrix/'vela_example_runner')==summary['runner_sha256']
    rows=summary['runs'];repeats=range(summary['repeats'])
    expected={(r,g,v) for r in repeats for g in (4,8) for v in ('native','baseline','combined')}
    assert len(rows)==len(expected) and {(r['repeat'],r['gate_V'],r['variant']) for r in rows}==expected
    report=dict(status='pass',runner_sha256=summary['runner_sha256'],repeated_numerics_exact=True,
                faster_than_baseline_every_pair=True,native_field_audit_required=True,gates={})
    for gate in (4,8):
        prior,_=inspect(a.reference/('r0_vg%d_contact'%gate))
        pairs=[]
        for repeat in repeats:
            group={r['variant']:r for r in rows if r['gate_V']==gate and r['repeat']==repeat}
            pair=dict(repeat=repeat)
            for variant in ('baseline','combined'):
                row=group[variant];directory=Path(row['directory'])
                assert sha(directory/'input.json')==row['input_sha256']
                assert sha(directory/'results/ledger.json')==row['ledger_sha256']
                actual,costs=inspect(directory);assert actual==prior
                assert costs==row['costs']
                ratio=row['external_wall_seconds']/group['native']['external_wall_seconds']
                pair[variant]=dict(wall_seconds=row['external_wall_seconds'],cpu_seconds=row['child_cpu_seconds'],
                    wall_ratio_to_native=ratio,cpu_ratio_to_native=row['child_cpu_seconds']/group['native']['child_cpu_seconds'],
                    costs=costs,performance=row['performance'])
                if ratio>1.5:report['status']='performance_gate_failed'
            pair['native']={k:group['native'][k] for k in ('external_wall_seconds','child_cpu_seconds')}
            pair['wall_reduction_fraction']=1-pair['combined']['wall_seconds']/pair['baseline']['wall_seconds']
            pair['cpu_reduction_fraction']=1-pair['combined']['cpu_seconds']/pair['baseline']['cpu_seconds']
            if pair['wall_reduction_fraction']<=0:report['faster_than_baseline_every_pair']=False
            pairs.append(pair)
        report['gates'][str(gate)]=dict(pairs=pairs,
            median_wall_reduction_fraction=statistics.median(p['wall_reduction_fraction'] for p in pairs))
    a.output.write_text(json.dumps(report,indent=2));print(json.dumps(report))


if __name__=='__main__':main()
