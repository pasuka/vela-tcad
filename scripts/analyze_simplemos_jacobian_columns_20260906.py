"""Classify column accuracy separately from FD signal/step qualification."""
from pathlib import Path
import math
import audit_simplemos_jacobian_columns_20260906 as m

a=m.a


def main():
    a.verify(m.OUT/'freeze.json');columns=[];blocks=[];cases=[]
    for c in a.read(m.OUT/'contract.json')['cases']:
        path=Path(c['config']);s=a.read(path.with_suffix('.status.json'));assert s['exit_code']==0
        rows=a.rows(path.parent/'columns.csv')
        assert len(rows)==len(c['columns'])*27==s['rows']
        assert len({(r['column_node'],r['column_block'],r['row_block'],r['step_V']) for r in rows})==len(rows)
        primary=[r for r in rows if abs(float(r['step_V'])/1e-6-1)<1e-10]
        fine={(r['column_node'],r['column_block'],r['row_block']):r for r in rows if abs(float(r['step_V'])/1e-7-1)<1e-10}
        controls=c['variant'].startswith('field_frozen')
        for cb in range(3):
            for rb in range(3):
                group=[r for r in primary if int(r['column_block'])==cb and int(r['row_block'])==rb]
                peak=max(float(r['fd_norm']) for r in group)
                active=stable=passed=0
                for r in group:
                    f=fine[r['column_node'],r['column_block'],r['row_block']]
                    signal=float(r['fd_norm'])
                    floor=max(1e-10*peak,1e-300)
                    is_active=signal>floor
                    is_stable=is_active and float(f['fd_step_difference_norm'])<=1e-4*signal
                    error=float(r['diff_norm'])/max(signal,1e-300)
                    okay=is_stable and error<=1e-3
                    active+=is_active;stable+=is_stable;passed+=okay
                    columns.append(dict(case=c['case'],variant=c['variant'],column_node=r['column_node'],column_block=cb,row_block=rb,
                        active=is_active,fd_stable=is_stable,relative_error=error,passed=okay,
                        max_diff_node=r['max_diff_node'],outside_pattern_relative=float(r['outside_pattern_norm'])/max(signal,1e-300),
                        intentionally_frozen_field=controls))
                norm=math.sqrt(math.fsum(float(r['fd_norm'])**2 for r in group))
                diff=math.sqrt(math.fsum(float(r['diff_norm'])**2 for r in group))
                missing=math.sqrt(math.fsum(float(r['outside_pattern_norm'])**2 for r in group))
                blocks.append(dict(case=c['case'],variant=c['variant'],column_block=cb,row_block=rb,fd_norm=norm,
                    relative_error=diff/max(norm,1e-300),outside_pattern_relative=missing/max(norm,1e-300),
                    active_columns=active,stable_columns=stable,passed_columns=passed,intentionally_frozen_field=controls))
        cases.append(dict(case=c['case'],variant=c['variant'],columns=3*len(c['columns']),exit_code=0))
    invariance=[]
    for case in sorted({c['case'] for c in cases}):
        root=m.LOCAL/case
        same=a.sha(root/'baseline_on/columns.csv.residual.csv')==a.sha(root/'baseline_off/columns.csv.residual.csv')
        assert same;invariance.append(dict(case=case,residual_identical=same))
    root=m.LOCAL/'m65_n23_vd_1p000000_endpoint'
    off=a.rows(root/'field_frozen_off/columns.csv');on=a.rows(root/'field_frozen_on/columns.csv')
    assert off==on
    for name,rows in (('columns',columns),('blocks',blocks),('cases',cases),('invariance',invariance)):a.write_csv(m.OUT/(name+'.csv'),rows)
    baseline=[r for r in columns if r['variant']=='baseline_on' and r['active']]
    a.write(m.OUT/'result.json',dict(input_hashes={a.rel(Path(__file__).resolve()):a.sha(Path(__file__).resolve())},
        active_baseline_blocks=len(baseline),stable_baseline_blocks=sum(r['fd_stable'] for r in baseline),passed_baseline_blocks=sum(r['passed'] for r in baseline),
        frozen_field_switch_identity=True,all_baseline_residuals_identical=True,production_promoted=False,
        boundary='All residual rows, locally complete selected columns only; inactive or unstable cross-block signals are not certified.'))
    print('Column audit analyzed',len(columns),'column/block records',flush=True)


if __name__=='__main__':main()
