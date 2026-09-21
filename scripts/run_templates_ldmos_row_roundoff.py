"""Replay six original Linux R7 prefixes and observe frozen rows without updates."""
import argparse
import json
import os
from pathlib import Path
import shutil
from run_templates_ldmos_poisson_initialization import run_point,restored_config
from run_templates_ldmos_density_projection import read,sha
from analyze_templates_ldmos_predictor_study import state_delta


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('baseline','probe','input','reference','output'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
    for name,path in (('original_probe',a.baseline),('audit_probe',a.probe)):shutil.copy2(str(path),str(root/name))
    cfg=read(a.input);full=read(a.reference)
    s=dict(schema='vela.hole_row_roundoff.v1',status='running',baseline_sha256=sha(a.baseline),audit_sha256=sha(a.probe),
           source_input_sha256=sha(a.input),source_reference_sha256=sha(a.reference),runs=[])
    def save():
        tmp=root/'summary.tmp';tmp.write_text(json.dumps(s,indent=2,allow_nan=False));os.replace(str(tmp),str(root/'summary.json'))
    save()
    for count in range(9,15):
        row,state=run_point(root/'original_probe',root/('prefix_%02d'%count),dict(cfg,diagnostic_newton_max_iterations=count))
        if row['exit_code'] or state is None:raise RuntimeError('Original prefix failed')
        if state['history']!=full['history'][:count] or state['newton_updates']!=count:raise ValueError('Frozen original prefix changed')
        frozen=restored_config(cfg,state)
        frozen.update(diagnostic_newton_max_iterations=0,diagnostic_hole_row_audit_nodes=[4634,4635,4636])
        audit,data=run_point(root/'audit_probe',root/('audit_%02d'%count),frozen)
        if audit['exit_code'] or data is None:raise RuntimeError('Read-only audit failed')
        delta=state_delta(state,data)
        if any(delta) or data['newton_updates']:raise ValueError('Observer changed physical state')
        if state['residual']!=data['residual']:raise ValueError('Observer residual differs from frozen original')
        # Deliberately separate the pre-existing load-time contact projection
        # from the raw frozen-state audit; it is NOT a Newton candidate.
        pcfg=restored_config(cfg,state);pcfg['diagnostic_newton_max_iterations']=0
        prep,projected=run_point(root/'original_probe',root/('projected_%02d'%count),pcfg)
        if prep['exit_code'] or projected is None:raise RuntimeError('Projection control failed')
        pd=state_delta(state,projected)
        if any(pd[1:]) or projected['newton_updates']:raise ValueError('Projection changed QF/T or performed updates')
        pcfg=restored_config(cfg,projected);pcfg.update(diagnostic_newton_max_iterations=0,diagnostic_hole_row_audit_nodes=[4634,4635,4636])
        pr,pa=run_point(root/'audit_probe',root/('audit_projected_%02d'%count),pcfg)
        if pr['exit_code'] or pa is None or pa['residual']!=projected['residual']:raise ValueError('Projected observer mismatch')
        s['runs'].append(dict(updates=count,prefix=row,audit=audit,prefix_exact=True,state_delta=delta,residual_exact=True,
            original_row_gate=state['carrier_row_gate'],audited_row_gate=data['carrier_row_gate'],
            projection=prep,projected_audit=pr,projection_state_delta=pd,projected_row_gate=pa['carrier_row_gate'],
            projected_block_gates=pa['electrical_block_gates']))
        save()
    s['status']='completed_frozen_audit';save()


if __name__=='__main__':main()
