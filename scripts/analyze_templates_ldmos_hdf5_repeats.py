"""Strict full-curve repeated state-format audit and paired timing summary."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics
from analyze_templates_ldmos_binary_state import analyze as paired_analysis


def validate_batch(batch):
    if batch.get('status')!='completed' or batch.get('points')!=31:
        raise ValueError('Completed full curves required')
    if batch.get('control')!='binary' or batch.get('candidate')!='hdf5':
        raise ValueError('Expected VDS1/HDF5 pairing')
    expected={(r,m,g) for r in range(3) for m in ('binary','hdf5') for g in (4,8)}
    cases=batch['cases'];keys=[(c['repeat'],c['mode'],c['gate']) for c in cases]
    if len(keys)!=12 or set(keys)!=expected or batch['planned_curves']!=12:
        raise ValueError('Missing, duplicate or extra curves')
    joint=batch['joint'];jkeys=[(j['repeat'],j['mode']) for j in joint]
    if len(jkeys)!=6 or set(jkeys)!={(r,m) for r in range(3) for m in ('binary','hdf5')}:
        raise ValueError('Incomplete joint gate checks')
    if any(not j['checks'] or not all(j['checks'].values()) for j in joint):
        raise ValueError('Joint gates failed')
    for c in cases:
        if c['status']!='completed' or not c['trajectory_exact'] or not c['audit']['integrity_pass']:
            raise ValueError('Unqualified curve')
        if c['audit']['exact_points']!=31 or c['bias_V']!=40.:
            raise ValueError('Wrong curve range')
        if not c['verdicts'] or not all(v['pass_all'] for v in c['verdicts'].values()):
            raise ValueError('Original curve gates failed')
    return cases


def analyze(root):
    import dd_state_binary as binary
    import dd_state_public as public
    batch=json.loads((root/'summary.json').read_text());cases=validate_batch(batch)
    result=paired_analysis(root)
    if any(p['counter_differences'] for p in result['pairs']):
        raise ValueError('Paired calculation counters changed')
    hashes={};frame_checks=[];sizes=[];trajectory={}
    for c in cases:
        key=(c['repeat'],c['mode'],c['gate']);directory=Path(c['directory'])
        ledger=json.loads((directory/'fixed/ledger.json').read_text())
        if ledger['status']!='completed':raise ValueError('Incomplete ledger')
        events=ledger['frame_events']
        if len(events)!=1 or any(e['status']!='qualified' or not all(e['checks'].values()) for e in events):
            raise ValueError('Frame transition not qualified')
        frame_checks.append(dict(repeat=c['repeat'],mode=c['mode'],gate=c['gate'],checks=events[0]['checks']))
        trajectory[key]=[(r['parent_V'],r['target_V'],r['cap_V'],r['stage'],r['frame_V'],r['Newton_updates']) for r in ledger['runs']]
        values=[]
        for point in ledger['exact_points']:
            path=Path(point['state'])
            if path.suffix=='.vds':rows=binary.read(path)
            elif path.suffix=='.h5':rows=public.read(path)
            else:
                with path.open(newline='') as f:rows=list(csv.DictReader(f))
            values.append((point['bias_V'],hashlib.sha256(binary.encode(rows)).hexdigest()))
        if len(values)!=31 or values[0][0]!=0. or values[-1][0]!=40.:raise ValueError('Incomplete exact states')
        hashes[key]=values
        files=[p for p in directory.rglob('*') if p.is_file() and p.suffix in ('.vds','.h5')]
        sizes.append(dict(repeat=c['repeat'],mode=c['mode'],gate=c['gate'],files=len(files),bytes=sum(p.stat().st_size for p in files)))
    for repeat in range(3):
        for gate in (4,8):
            if hashes[repeat,'binary',gate]!=hashes[repeat,'hdf5',gate]:raise ValueError('Cross-format saved fields differ')
    for repeat in (1,2):
        for mode in ('binary','hdf5'):
            for gate in (4,8):
                if hashes[0,mode,gate]!=hashes[repeat,mode,gate] or trajectory[0,mode,gate]!=trajectory[repeat,mode,gate]:
                    raise ValueError('Repeated states or trajectory differ')
                a=next(c for c in cases if (c['repeat'],c['mode'],c['gate'])==(0,mode,gate))
                b=next(c for c in cases if (c['repeat'],c['mode'],c['gate'])==(repeat,mode,gate))
                if a['performance']['counters']!=b['performance']['counters']:raise ValueError('Repeated counters changed')
    timing=[]
    for mode in ('binary','hdf5'):
        for gate in (4,8):
            v=[next(c['wall_seconds'] for c in cases if (c['repeat'],c['mode'],c['gate'])==(r,mode,gate)) for r in range(3)]
            timing.append(dict(mode=mode,gate=gate,seconds=v,median=statistics.median(v),minimum=min(v),maximum=max(v),
                mean=statistics.mean(v),sample_stdev=statistics.stdev(v)))
    result.update(numerical_integrity_pass=True,full_curves=12,exact_points=372,joint_passes=6,
        cross_format_exact_point_pairs=186,repeat_exact_point_pairs=248,all_saved_fields_exact=True,
        frame_checks=frame_checks,state_sizes=sizes,timing=timing)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--batch',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();result=analyze(a.batch)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('runtime_sha256','frame_checks','state_sizes','pairs')},indent=2))
