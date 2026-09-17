"""Audit all paired native full curves/fields against the frozen D0 reference."""
import argparse
import json
from pathlib import Path

from analyze_templates_ldmos_contact_sweep import audit_native
from run_templates_ldmos_contact_sweep import read,sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix',type=Path,required=True)
    parser.add_argument('--reference',type=Path,default=Path('reference_staging/templates_ldmos_d0_electrothermal_20260912'))
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();summary=read(args.matrix/'summary.json')
    assert summary['status']=='complete' and summary['repeats']>=2
    native=[r for r in summary['runs'] if r['variant']=='native']
    expected={(r,g) for r in range(summary['repeats']) for g in (4,8)}
    assert len(native)==len(expected) and {(r['repeat'],r['gate_V']) for r in native}==expected
    reports=[]
    for row in native:
        directory=args.matrix/('r%d_vg%d_native'%(row['repeat'],row['gate_V']))
        check=audit_native(directory,args.reference,directory/'normalized')
        reports.append(dict(repeat=row['repeat'],gate_V=row['gate_V'],**check))
    report=dict(status='pass',scope='All paired native curves and potential/QF/temperature fields exactly match frozen qualified reference; geometry identities checked.',
                summary_sha256=sha(args.matrix/'summary.json'),runs=reports)
    args.output.write_text(json.dumps(report,indent=2));print(json.dumps(dict(status='pass',native_runs=len(reports))))


if __name__=='__main__':main()
