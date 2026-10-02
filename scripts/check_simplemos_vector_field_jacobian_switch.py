"""Read-only control of the field-Jacobian switch on the frozen vector-drive case."""
import json
import subprocess
from pathlib import Path
import check_simplemos_actual_direction_jvp as audit

LOCAL=audit.LOCAL/'field_switch_off'
CONTRACT=audit.OUT/'field_switch_contract.json'


def main():
    audit.prior.audit.verify(audit.FREEZE)
    cfg=audit.prior.audit.read(audit.LOCAL/'1e13/config.json')
    cfg['solver']['mobility']['jacobian_field_derivatives']=False
    cfg['output_csv']=str(LOCAL/'jvp.csv')
    audit.prior.audit.write(LOCAL/'config.json',cfg)
    paths=[Path(__file__),LOCAL/'config.json',audit.FREEZE,audit.RUNNER,
        audit.REPO/'include/vela/physics/MobilityModel.h',audit.REPO/'src/physics/MobilityModel.cpp',
        audit.REPO/'src/equation/CoupledDDAssembler.cpp']
    audit.prior.audit.write(CONTRACT,{'status':'frozen_before_read_only_control',
        'single_axis':'Explicit jacobian_field_derivatives=false against omitted default true; transport_cell_vector unchanged.',
        'scope':'No solve, no model selection change. Compare all residual FD rows and Jv rows; this is not a candidate correction.',
        'input_hashes':{audit.prior.audit.rel(p):audit.prior.audit.sha(p) for p in paths}})
    p=subprocess.run([str(audit.RUNNER),'--config',str(LOCAL/'config.json'),'--log','off'],env=audit.prior.env(),capture_output=True,text=True)
    (LOCAL/'stdout.txt').write_text(p.stdout);(LOCAL/'stderr.txt').write_text(p.stderr)
    if p.returncode:raise RuntimeError(p.stderr[-1000:])
    a=audit.prior.audit.rows(audit.LOCAL/'1e13/jvp.csv');b=audit.prior.audit.rows(LOCAL/'jvp.csv')
    assert len(a)==len(b)
    samefd=all(x['fd']==y['fd'] and x['actual_endpoint_fd']==y['actual_endpoint_fd'] for x,y in zip(a,b))
    samej=all(x['analytic']==y['analytic'] for x,y in zip(a,b))
    result={'rows':len(a),'residual_differences_identical':samefd,'Jv_identical':samej,
        'interpretation':'Switch alone cannot repair observed vector-drive defect if both operators are identical. Current default is true; do not diagnose missing user configuration.',
        'nonlinear_solves':0,'m82_released':False}
    audit.prior.audit.write(audit.OUT/'field_switch_result.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
