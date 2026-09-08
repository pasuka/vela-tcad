"""Retain failed original preflight; apply the unchanged experiment to stable J."""
import argparse
from pathlib import Path
import validate_simplemos_joint_geometry_20260907 as v
import check_simplemos_joint_geometry_20260907 as check
import build_simplemos_joint_stable_20260907 as b

v.b=b;v.LOCAL=b.LOCAL;v.OUT=b.OUT
check.LOCAL=b.LOCAL;check.OUT=b.OUT


def prepare():
    v.prepare()
    files=[Path(__file__).resolve(),Path(check.__file__).resolve(),Path(b.__file__).resolve(),b.HEADER,b.old.HEADER,
        b.old.OUT/'freeze.json',b.old.OUT/'preflight_checks.csv',b.old.OUT/'preflight_jvp.csv',b.old.OUT/'preflight_ports.csv',
        b.LOCAL/'compile_commands.json']
    for cmd in b.a.read(b.LOCAL/'compile_commands.json'):
        files += [Path(x) for x in cmd if x.endswith('.cpp')]
    b.a.write(b.OUT/'numerical_addendum.json',dict(status='frozen_before_stable_execution',
        description='Only the analytic equal-ni/no-BGN Boltzmann SG psi derivative is evaluated via stable flux times logarithmic Bernoulli derivative. Original physical residual, geometry experiment and acceptance gates unchanged.',
        retained_failure='Original baseline and three variants fail hole-row/psi-column JVP. No original full replacement DC has run.',
        scope='Masetti, fixed 300 K and no field-dependent mobility; no production change or claim about clipping-active states.',
        input_hashes={b.a.rel(f):b.a.sha(f) for f in files}))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','preflight','check','run'))
    action=parser.parse_args().action
    if action=='prepare':prepare()
    elif action=='check':check.main()
    else:getattr(v,action)()
