"""Re-solve LDMOS equilibrium after a statistics/material change, with original gates.

All carrier contacts are at zero voltage at these drain-zero gate biases.
Constant zero QFs express equilibrium; Poisson updates the charge/potential.
Only the subsequent full Newton check can qualify a seed for a drain sweep.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import types
import run_templates_ldmos_linked_d5 as linked


def equilibrium_guess(rows):
    result=[dict(row) for row in rows]
    for row in result:
        for key in ('phin','phip','electron_qf_increment_V','hole_qf_increment_V',
                    'electron_qf_reference_V','hole_qf_reference_V'):
            if key in row:row[key]='0'
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--source-bundle',type=Path,default=Path('reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_units_inputs.json'))
    parser.add_argument('--template',type=Path,default=Path('reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_fermi_accurate_config.json'))
    args=parser.parse_args();root=Path(__file__).resolve().parents[1]
    bundle=linked.read(args.source_bundle);manifest=linked.read(args.manifest)
    runner=Path(manifest['runner']);assert linked.digest(runner)==manifest['runner_sha256']
    for path,expected in bundle['files'].items():
        if linked.digest(root/path)!=expected:raise ValueError(f'Input changed: {path}')
    base=linked.read(args.template);linked.PHYSICS_PROFILE='D5';linked.assert_linked_physics_contract(base)
    for c in base['contacts']:
        if c['name']!='gate' and c.get('bias',0.)!=0.:raise ValueError('Requires zero carrier-contact biases')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    linked.HERE=out;linked.ROOT=root;linked.BASE=base;linked.RUNNER=runner
    linked.EXPECTED=manifest['runner_sha256'];linked.WORKER=None
    linked.ENV=dict(os.environ,VELA_LINEAR_SOLVER='sparselu')
    if os.name=='nt':
        linked.ENV['PATH']=r'D:\msys64\ucrt64\bin'+os.pathsep+linked.ENV.get('PATH','')
    linked.write(out/'inputs.json',dict(bundle=bundle,template=base,manifest=manifest,
        preparation='Constant zero QF guess, Poisson, then full original Newton gates'))
    summary={}
    for gate in (4,8):
        linked.args=types.SimpleNamespace(gate=gate,worker=False)
        source=root/bundle['seeds'][str(gate)]
        seed=out/f'vg{gate}_equilibrium_guess.csv'
        linked.csv_out(seed,equilibrium_guess(linked.rows(source)))
        dest=out/f'vg{gate}_poisson';dest.mkdir()
        cfg=linked.prepare_config(base,dest,seed,0.,.1,gate,0.)
        cfg['solver']['method']='poisson_only';linked.write(dest/'control.json',cfg)
        with (dest/'stdout.log').open('w') as log:
            code=subprocess.run([str(runner),'--config',str(dest/'control.json')],
                cwd=dest,env=linked.ENV,stdout=log,stderr=subprocess.STDOUT,timeout=600).returncode
        if code:raise RuntimeError(f'Vg={gate} Poisson preparation failed: {code}')
        result=linked.run_child(f'vg{gate}_qualified',dest/'state.csv',0.,0.,.1,0.)
        qualified=linked.good(result,out/f'vg{gate}_qualified',0.)
        state=out/f'vg{gate}_qualified/state.csv'
        summary[str(gate)]=dict(qualified=qualified,source_sha256=linked.digest(source),
            state=str(state),state_sha256=linked.digest(state),
            poisson_attempts=linked.rows(dest/'newton_attempts.csv'),
            qualification_attempts=result['attempts'])
        linked.write(out/'summary.json',summary)
        if not qualified:raise RuntimeError(f'Vg={gate} failed full Newton qualification')
        print(f'Vg={gate}: equilibrium seed qualified',flush=True)


if __name__=='__main__':main()
