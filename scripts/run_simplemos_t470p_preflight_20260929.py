"""Cross-platform replay and reclosure of frozen complete split checkpoints.

This is a build qualification, not a second independent initialization.
"""
import argparse
import copy
import json
import math
from pathlib import Path
import re
import sys

import simplemos_hfs_cloud_20260926 as h
from check_simplemos_cold_overlap_20260927 import load, potential
from decimal import localcontext


def run(base,runner=None):
    base=base.resolve();runner=Path(runner).resolve() if runner else base/'source/build-release/vela_example_runner.exe'
    root=base/'preflight';root.mkdir(exist_ok=True)
    seal=dict(binary_sha256=h.sha(runner),driver_sha256=h.sha(__file__),
              inputs_sha256=h.sha(base/'data/inputs/hashes.json'),
              seed_sha256={p.name:h.sha(p) for p in sorted((base/'data/seeds').glob('*.h5'))})
    if (root/'seal.json').exists():
        if h.read(root/'seal.json')!=seal:raise ValueError('Preflight identity changed')
    else:h.write(root/'seal.json',seal)
    executor=object.__new__(h.Run);executor.runner=runner
    frozen={(r['case'],int(r['index'])):r for r in h.rows(base/'data/frozen_points.csv')}
    results=[]
    for seed in sorted((base/'data/seeds').glob('*.h5')):
        match=re.fullmatch(r'(n\d+)_vd_(0p05|1)_vg_(\d+)\.h5',seed.name)
        if not match:raise ValueError(seed.name)
        device,vdtag,index=match.groups();index=int(index);vg=index*.05;vd=float(vdtag.replace('p','.'))
        case=device+'_vd_'+vdtag;dest=root/seed.stem;dest.mkdir(exist_ok=True)
        cfg=h.read(base/'data/inputs'/device/'template.json')
        for name in ('mesh_file','node_doping_file','materials_file'):cfg[name]=str(base/'data/inputs'/cfg[name])
        for contact in cfg['contacts']:contact['bias']=vd if contact['name']=='drain' else vg if contact['name']=='gate' else 0.
        probe=copy.deepcopy(cfg);probe.update(simulation_type='terminal_current_functional_probe',state_format='hdf5',state_file=str(seed),contact='drain')
        h.write(dest/'replay.json',probe);replayed=executor.execute(dest/'replay.json')
        expected=float(frozen[case,index]['current_A_per_um'])
        measured=replayed.get('contact_current_extractor_A_per_um',math.nan)
        relative=abs(measured/expected-1)
        record=dict(case=case,vg=vg,replay_Id_relative=relative,
                    replay_Id_absolute=abs(measured-expected),replay_passed=replayed['exit_code']==0 and math.isfinite(relative) and relative<=1e-12)
        if not record['replay_passed']:
            h.write(dest/'result.json',record);raise ValueError(f'Frozen split replay failed: {record}')
        cfg.update(simulation_type='newton_solve_from_state',state_format='hdf5',state_file=str(seed),output_state_file=str(dest/'state.h5'))
        h.write(dest/'config.json',cfg);status=executor.execute(dest/'config.json')
        if status['exit_code']!=0 or not status.get('converged'):raise ValueError(f'Reclosure failed: {case}')
        x,y=load(seed),load(dest/'state.h5');geo=h.read(base/'data/inputs'/device/'geometry.json')
        if x[1]['mesh_sha256']!=y[1]['mesh_sha256']:raise ValueError('Mesh identity mismatch')
        with localcontext() as context:
            context.prec=100
            diff={f+'_max_V':max(float(abs(potential(x,f,i)-potential(y,f,i))) for i in geo['free_si']) for f in ('psi','phin','phip')}
        density=max(abs(float(x[0][f][i])/float(y[0][f][i])-1) for f in ('electrons_m3','holes_m3') for i in geo['free_si'])
        currents=status['contact_currents_A_per_um'];current=currents['drain']
        kcl=abs(math.fsum(currents.values()))/abs(current)
        row=status['carrier_row_convergence']
        record.update(**diff,density_max_relative=density,Id_relative=abs(current/expected-1),kcl_over_Id=kcl,
                      convergence_reason=status.get('convergence_reason'),final_residual=status.get('final_residual'),carrier_row=row)
        # The full external all-row/port qualification is performed by the
        # matrix harness; this preflight preserves its corresponding gates.
        record['passed']=all(math.isfinite(v) and v<=1e-6 for v in diff.values()) and density<=1e-4 and record['Id_relative']<=1e-6 and kcl<=1e-8 and row.get('satisfied',False)
        h.write(dest/'result.json',record);results.append(record)
        h.write(root/'progress.json',dict(completed=len(results),last=record))
        print(case,vg,'passed',record['passed'],flush=True)
        if not record['passed']:raise ValueError(record)
    report=dict(points=len(results),passed=len(results)==6 and all(r['passed'] for r in results),
                scope='Six existing low-Vd outlier checkpoints: exact split replay and same-seed reclosure, not independent dual initialization',results=results)
    h.write(root/'summary.json',report)
    if not report['passed']:raise ValueError('Incomplete preflight')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--base',type=Path,required=True);parser.add_argument('--runner',type=Path);args=parser.parse_args()
    try:run(args.base,args.runner)
    except Exception as exc:
        h.write(args.base/'preflight/failure.json',dict(error=repr(exc)));raise
