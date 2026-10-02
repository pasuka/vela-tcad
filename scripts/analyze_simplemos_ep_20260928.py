"""Compare EP PLT components against M60 and frozen Vela; no reference updates."""
import csv
import hashlib
import json
import math
from pathlib import Path
from sentaurus_import import parse_quoted_list, parse_values_block

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'build/simplemos_ep_20260928'
CONTACTS = ('source', 'drain', 'gate', 'substrate')

def readcsv(p):
    with p.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

def table(name, rows):
    if not rows:
        return
    with (OUT/name).open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def currents(path):
    text = path.read_text()
    names = parse_quoted_list(text, 'datasets')
    rows = [dict(zip(names, r)) for r in parse_values_block(text, len(names))]
    result = {}
    for i in range(51):
        hits = [r for r in rows if abs(r['gate OuterVoltage']-i*.05)<1e-10]
        if len(hits)!=1:
            raise ValueError(f'{path}: Vg index {i} has {len(hits)} rows')
        assert abs(hits[0]['drain OuterVoltage']-.05)<1e-10
        result[i] = hits[0]
    return result

def main():
    manifest = json.loads((OUT/'manifest.json').read_text())
    status = dict(line.split() for line in (OUT/'raw/status.tsv').read_text().splitlines())
    data, checks, components = {}, [], []
    for c in manifest['cases']:
        name = c['name']
        folder = OUT/'raw/bundle'/c['device']
        log = folder/f'{name}.console.log'
        text = log.read_text(errors='replace') if log.exists() else ''
        evidence = [line.strip() for line in text.splitlines()
            if any(word in line.lower() for word in ('precision', 'super', 'rhsmin', 'digits', 'linear solver'))]
        check = dict(name=name, exit_code=status.get(name, 'pending'),
            log_sha256=hashlib.sha256(log.read_bytes()).hexdigest() if log.exists() else '',
            math_backend_echo=' | '.join(evidence), points=0)
        if status.get(name)=='0':
            values = currents(folder/f'IdVg_{name}_des.plt')
            check['points'] = len(values)
            data[c['device'], c['algorithm']] = values
            for i, r in values.items():
                for contact in CONTACTS:
                    e,h,t=(r[f'{contact} {k}'] for k in ('eCurrent','hCurrent','TotalCurrent'))
                    components.append(dict(device=c['device'], algorithm=c['algorithm'],
                        index=i, vg=i*.05, contact=contact, electron=e, hole=h, total=t,
                        component_identity_absolute=e+h-t))
        checks.append(check)
    frozen = {(r['device'], round(float(r['vg'])/.05)):r
        for r in readcsv(ROOT/'build/outlier_analysis_20260928/all_points.csv') if float(r['vd'])==.05}
    prior = {(r['device'], int(r['index'])):r for r in readcsv(ROOT/'build/outlier_m60_review_20260928/joined_816.csv') if float(r['vd'])==.05}
    results=[]
    for device in ('n23','n24','n17','n19'):
        if any((device, alg) not in data for alg in ('default','direct')):
            continue
        for i in range(51):
            a,b=data[device,'default'][i],data[device,'direct'][i]
            v=float(frozen[device,i]['current_A_per_um'])
            s=float(prior[device,i]['tight_default_Id_A_per_um'])
            current=a['drain TotalCurrent']
            gap=a['substrate eCurrent']-b['substrate eCurrent']
            kcl=math.fsum(a[f'{c} TotalCurrent'] for c in CONTACTS)
            results.append(dict(device=device,index=i,vg=i*.05,vela_Id=v,m60_default_Id=s,
                ep_default_Id=current,ep_direct_Id=b['drain TotalCurrent'],
                vela_ep_error_percent=100*(v/current-1), ep_m60_change_percent=100*(current/s-1),
                default_KCL=kcl,default_KCL_relative=kcl/abs(current),
                direct_KCL=math.fsum(b[f'{c} TotalCurrent'] for c in CONTACTS),
                substrate_gap=gap,substrate_gap_relative=gap/abs(current),
                direct_default_Id_relative=b['drain TotalCurrent']/current-1))
    table('run_checks.csv',checks)
    table('terminal_components.csv',components)
    table('comparison.csv',results)
    summary=dict(completed_runs=sum(c['points']==51 for c in checks),
        expected_runs=8, paired_points=len(results), acceptance_changed=False,
        actual_precision_and_backend_require_log_review=True,
        max_abs_vela_ep_error_percent=max((abs(r['vela_ep_error_percent']) for r in results),default=None))
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))

if __name__=='__main__':
    main()
