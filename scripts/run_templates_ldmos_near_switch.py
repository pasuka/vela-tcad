"""Full frozen-seed trajectories for the opt-in near-steady QF switch."""
import argparse
import copy
import json
import os
from pathlib import Path
import shutil

from run_templates_ldmos_density_projection import read,sha
from run_templates_ldmos_poisson_initialization import run_point,restored_config,preparation_gate,total_cost,variant_config
from analyze_templates_ldmos_predictor_study import state_delta
from run_templates_ldmos_electrothermal_curve import state_gate


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--control',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    control=read(args.control/'summary.json')
    if control['status']!='completed_diagnostic_matrix' or sha(args.manifest)!=control['manifest_sha256']:
        raise ValueError('Control or manifest invalid')
    cases=[c for c in read(args.manifest) if c['variant']=='linear_baseline'
           and (c['gate_V'],c['baseline_index']) in ((4,11),(8,9))]
    if len(cases)!=2:raise ValueError('Expected two original frozen seeds')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    binary=out/args.probe.name;shutil.copy2(args.probe,binary)
    for path in args.probe.parent.glob('*.dll'):shutil.copy2(path,out/path.name)
    if os.name=='nt':
        import ctypes
        ctypes.windll.kernel32.SetErrorMode(0x0001|0x0002|0x8000)
    summary=dict(schema='vela.near_steady_switch.v1',status='running',binary_sha256=sha(binary),
        runtime_sha256={p.name:sha(p) for p in out.glob('*.dll')},manifest_sha256=sha(args.manifest),
        control_directory=str(args.control.resolve()),control_summary_sha256=sha(args.control/'summary.json'),
        preparations=[],runs=[])
    def save():
        tmp=out/'summary.tmp';tmp.write_text(json.dumps(summary,indent=2,allow_nan=False),encoding='utf-8');os.replace(tmp,out/'summary.json')
    save()
    for c in cases:
        vg=c['gate_V'];name='vg%d_p%d'%(vg,c['baseline_index']);original=read(Path(c['input']));reference=read(Path(c['expected']))
        old=next(r for r in control['runs'] if r['variant']=='baseline' and r['vg']==vg)
        if sha(Path(c['input']))!=old['source_input_sha256'] or sha(Path(c['expected']))!=old['source_output_sha256']:
            raise ValueError('Frozen state changed')
        def run(variant,cfg,prior=()):
            row,data=run_point(binary,out/(variant+'_'+name),cfg)
            row.update(variant=variant,vg=vg,source_input_sha256=sha(Path(c['input'])),
                       reference_path=c['expected'],reference_sha256=sha(Path(c['expected'])),qualified=False)
            if data is not None:
                gate=state_gate(data,c['bias_V']);delta=state_delta(data,reference)
                row.update(gate=gate,max_state_delta=delta,switch=data.get('near_steady_qf_switch'),
                    qualified=row['exit_code']==0 and gate['pass_gate'] and all(d<=t for d,t in zip(delta,(1e-8,1e-8,1e-8,1e-7))),
                    charged_cost=total_cost(*prior,row))
            summary['runs'].append(row);save()
            if data is None or row['exit_code']:raise RuntimeError('Probe failed')
            return row,data
        baseline,base=run('baseline',original)
        if not baseline['qualified']:raise ValueError('Original baseline failed')
        cfg=variant_config(original,'ptc_defect_model_s1');cfg['diagnostic_near_steady_qf_switch']=True
        run('switch_raw',cfg)
        prep,prepared=run_point(binary,out/('prepare_'+name),dict(original,solve_mode='poisson'))
        if prep['exit_code'] or prepared is None:raise RuntimeError('Preparation execution failed')
        restored=restored_config(original,prepared)
        audit,audited=run_point(binary,out/('audit_'+name),dict(restored,diagnostic_newton_max_iterations=0))
        if audit['exit_code'] or audited is None:raise RuntimeError('Preparation audit failed')
        gate=preparation_gate(original,prepared,audited)
        summary['preparations'].append(dict(vg=vg,preparation=prep,audit=audit,gate=gate));save()
        if not gate['pass_gate']:continue
        cfg=variant_config(restored,'ptc_defect_model_s1');cfg['diagnostic_near_steady_qf_switch']=True
        run('switch_prepared',cfg,(prep,audit))
    summary['status']='completed_diagnostic_matrix';save()


if __name__=='__main__':main()
