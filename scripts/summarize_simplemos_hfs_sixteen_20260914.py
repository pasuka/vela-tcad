"""Combine same-binary high and low HFS controls without hiding failures."""
from pathlib import Path
import simplemos_hfs_curves_20260914 as e
a,d=e.a,e.d

def main():
    for p in (e.O/'control_evidence.json',e.O/'supplemental_source_evidence.json',e.O/'early_low_evidence.json',e.q.O/'comparison_evidence.json'):a.verify(p)
    high=a.rows(e.q.O/'comparison.csv'); low=a.rows(e.O/'low_control_comparison.csv');native=a.rows(e.O/'native_points.csv')
    rows=[]
    for r in high:
        rows.append(dict(device=r['device'],vd=float(r['vd']),vg=float(r['vg']),Id_native=float(r['native_Id_A_per_um']),Id_vela=float(r['vela_Id_A_per_um']),error_percent=float(r['error_percent']),dual_Id_relative=float(r['dual_Id_relative']),**{k:float(r[k]) for k in ('psi_max_V','phin_max_V','phip_max_V','density_max_relative')},qualified=r['dual_qualified']=='True'))
    attempts=a.rows(e.O/'low_control_attempts.csv')
    for r in low:
        ref=next(n for n in native if n['case']==r['case'] and n['index']==r['index'])
        selected=next(x for x in attempts if x['case']==r['case'] and x['index']==r['index'] and x['arm']=='enormal_seed' and x['qualified']=='True')
        rows.append(dict(device=r['device'],vd=float(r['vd']),vg=float(r['vg']),Id_native=float(ref['Id_A_per_um']),Id_vela=float(selected['current_A_per_um']),**{k:float(r[k]) for k in ('error_percent','dual_Id_relative','psi_max_V','phin_max_V','phip_max_V','density_max_relative')},qualified=r['qualified']=='True'))
    rows.sort(key=lambda r:(r['device'],r['vd'],r['vg']))
    attempts+=a.rows(e.q.O/'attempts.csv')
    summary=dict(points=len(rows),qualified_points=sum(r['qualified'] for r in rows),attempts=len(attempts),failed_attempts=sum(r['qualified']!='True' for r in attempts),reloads=sum(int(r['attempt'])>0 for r in attempts),max_row_ratio=max(float(r['max_row_ratio']) for r in attempts),max_port_relative=max(float(r['port_relative']) for r in attempts),max_abs_Id_error_percent=max(abs(r['error_percent']) for r in rows),max_dual_Id_relative=max(r['dual_Id_relative'] for r in rows),max_dual_phi_V=max(r[k] for r in rows for k in ('psi_max_V','phin_max_V','phip_max_V')),max_dual_density_relative=max(r['density_max_relative'] for r in rows),runner=str(e.q.RUNNER),runner_sha256=a.sha(e.q.RUNNER),acceptance_changed=False,all_qualified=len(rows)==16 and all(r['qualified'] for r in rows))
    a.write_csv(e.O/'sixteen_comparison.csv',rows);a.write(e.O/'sixteen_summary.json',summary)
    d.matrix.freeze(e.O/'sixteen_evidence.json',[Path(__file__).resolve(),e.O/'control_evidence.json',e.O/'supplemental_source_evidence.json',e.O/'early_low_evidence.json',e.q.O/'comparison_evidence.json',e.O/'sixteen_comparison.csv',e.O/'sixteen_summary.json'])
    print(summary,flush=True)

if __name__=='__main__':main()
