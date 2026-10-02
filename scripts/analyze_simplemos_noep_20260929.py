"""Audit the four-run precision ablation, including failed/partial native runs.

No reference, acceptance criterion, or solver input is modified. Currents are
native A/um, with e+h=total; comparisons require an actual matching PLT point.
"""
import csv
import hashlib
import json
import math
import re
from pathlib import Path

from sentaurus_import import parse_quoted_list, parse_values_block

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'build/simplemos_noep_20260929'
CONTACTS = ('source', 'drain', 'gate', 'substrate')
NUM = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'


def readcsv(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def table(name, rows):
    if rows:
        with (OUT / name).open('w', encoding='utf-8', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def values(path):
    if not path.exists():
        return {}, 0
    text = path.read_text()
    names = parse_quoted_list(text, 'datasets')
    rows = [dict(zip(names, r)) for r in parse_values_block(text, len(names))]
    result = {}
    for row in rows:
        vg, vd = row['gate OuterVoltage'], row['drain OuterVoltage']
        index = round(vg / .05)
        if abs(vg-index*.05) < 1e-10 and abs(vd-.05) < 1e-10 and 0 <= index <= 50:
            if index in result:
                raise ValueError(f'Duplicate exact gate point: {path}, {index}')
            if not all(math.isfinite(v) for v in row.values()):
                raise ValueError(f'Nonfinite PLT row: {path}, {index}')
            result[index] = row
    return result, len(rows)


def log_audit(text, name):
    iterations, blocks = [], []
    context, block = '', None
    for line_no, line in enumerate(text.splitlines(), 1):
        if 'Computing step from' in line or 'Solving ' in line:
            context = line.strip()
        if 'Iteration' in line and '|Rhs|' in line:
            block = dict(name=name, block=len(blocks), context=context,
                         header_line=line_no, last_iteration=None, last_rhs=None,
                         last_error=None, stop_reason='', gate_V=None, drain_V=None)
            blocks.append(block)
        if block is not None and re.fullmatch(r'\s*\d+\s+' + NUM + r'(?:\s+' + NUM + r')+\s*', line):
            tokens = line.split()
            # Initial residual rows have 3 columns; Newton rows have 8.
            if len(tokens) in (3, 8):
                item = dict(name=name, block=block['block'], line=line_no,
                            iteration=int(tokens[0]), rhs=float(tokens[1]),
                            step=float(tokens[3]) if len(tokens)==8 else None,
                            error=float(tokens[4]) if len(tokens)==8 else None)
                iterations.append(item)
                block.update(last_iteration=item['iteration'], last_rhs=item['rhs'], last_error=item['error'])
        contact = re.fullmatch(r'\s*(gate|drain)\s+(' + NUM + r')(?:\s+' + NUM + r'){3}\s*', line)
        if block is not None and contact:
            block[contact[1] + '_V'] = float(contact[2])
        if block is not None and any(s in line.lower() for s in
                ('error smaller than', 'rhs smaller than', 'iterations exceeded',
                 'convergence failed', 'step-size less', 'step size less', 'not converg')):
            block['stop_reason'] += line.strip() + ' | '
    def seconds(label):
        found = re.search(re.escape(label) + r':\s*(' + NUM + r')\s*s', text)
        return float(found[1]) if found else None
    echoes = [line.strip() for line in text.splitlines() if any(s in line.lower() for s in
        ('precision', 'relative error', 'mininum |rhs|', 'minimum |rhs|', 'linear solver :'))]
    checks = dict(math_backend_echo=' | '.join(echoes),
        wallclock_s=seconds('wallclock'), cpu_s=seconds('total cpu'),
        nonlinear_blocks=len(blocks), newton_steps=sum(r['iteration'] > 0 for r in iterations),
        final_rhs=iterations[-1]['rhs'] if iterations else None,
        final_error=iterations[-1]['error'] if iterations else None,
        curve_finished='Curve trace finished.' in text,
        last_stop_reason=blocks[-1]['stop_reason'] if blocks else '')
    return checks, blocks, iterations


def main():
    manifest = json.loads((OUT/'manifest.json').read_text())
    checks, blocks, iterations, components, comparisons, pairs = [], [], [], [], [], []
    data = {}
    frozen = {(r['device'], int(r['index'])):r for r in readcsv(
        ROOT/'build/outlier_analysis_20260928/all_points.csv') if float(r['vd'])==.05}
    m60 = {(r['device'], int(r['index'])):r for r in readcsv(
        ROOT/'build/outlier_m60_review_20260928/joined_816.csv') if float(r['vd'])==.05}
    for case in manifest['cases']:
        name, device, alg = case['name'], case['device'], case['algorithm']
        folder = OUT/'raw/bundle'/device
        log = folder/f'{name}.console.log'
        text = log.read_text(errors='replace')
        checked, b, it = log_audit(text, name)
        blocks.extend(b); iterations.extend(it)
        rows, raw_count = values(folder/f'IdVg_{name}_des.plt')
        data[device, alg] = rows
        exit_code = int((OUT/'raw'/f'{name}.exitcode').read_text())
        checked.update(name=name, exit_code=exit_code, exact_points=len(rows),
            raw_plt_points=raw_count, largest_exact_vg=max(rows)*.05 if rows else None,
            log_sha256=hashlib.sha256(log.read_bytes()).hexdigest())
        checks.append(checked)
        ep, count = values(ROOT/'build/simplemos_ep_20260928/raw/bundle'/device/
                           f"IdVg_{case['source_ep_name']}_des.plt")
        assert count == len(ep) == 51
        for index, row in rows.items():
            current = row['drain TotalCurrent']
            v = float(frozen[device,index]['current_A_per_um'])
            strict = float(m60[device,index]['tight_default_Id_A_per_um'])
            e = ep[index]['drain TotalCurrent']
            kcl = math.fsum(row[f'{c} TotalCurrent'] for c in CONTACTS)
            comparisons.append(dict(device=device,algorithm=alg,index=index,vg=index*.05,
                noep_Id=current,ep_same_observer_Id=e,m60_default_Id=strict,vela_Id=v,
                vela_noep_error_percent=100*(v/current-1),
                vela_m60_default_error_percent=100*(v/strict-1),
                vela_ep_same_observer_error_percent=100*(v/e-1),
                noep_ep_change_percent=100*(current/e-1),
                kcl=kcl,kcl_relative=kcl/abs(current)))
            for contact in CONTACTS:
                ec,hc,total=(row[f'{contact} {s}'] for s in ('eCurrent','hCurrent','TotalCurrent'))
                components.append(dict(device=device,algorithm=alg,index=index,vg=index*.05,
                    contact=contact,electron=ec,hole=hc,total=total,identity_error=ec+hc-total))
    for device in ('n23','n24'):
        a,b=data[device,'default'],data[device,'direct']
        for index in sorted(a.keys() & b.keys()):
            denominator=abs(a[index]['drain TotalCurrent'])
            gap=a[index]['substrate eCurrent']-b[index]['substrate eCurrent']
            pairs.append(dict(device=device,index=index,vg=index*.05,
                substrate_e_gap=gap,substrate_e_gap_over_default_Id=gap/denominator,
                direct_default_Id_relative=b[index]['drain TotalCurrent']/a[index]['drain TotalCurrent']-1))
    table('run_checks.csv', checks)
    table('nonlinear_blocks.csv', blocks)
    table('newton_iterations.csv', iterations)
    table('terminal_components.csv', components)
    table('comparison.csv', comparisons)
    table('observer_pairs.csv', pairs)
    summary=dict(expected_runs=4,finished_processes=len(checks),
        successful_exit_codes=sum(r['exit_code']==0 for r in checks),
        complete_curves=sum(r['exit_code']==0 and r['exact_points']==51 for r in checks),
        compared_points=len(comparisons),paired_points=len(pairs),
        acceptance_changed=False,simulation_retried=False,
        max_abs_default_vela_error_percent=max((abs(r['vela_noep_error_percent'])
            for r in comparisons if r['algorithm']=='default'),default=None),
        max_abs_default_kcl_relative=max((abs(r['kcl_relative'])
            for r in comparisons if r['algorithm']=='default'),default=None),
        accepted_blocks=len(blocks),
        update_stopped_blocks=sum('Error smaller than 1' in r['stop_reason'] for r in blocks),
        blocks_with_rhs_below_setting=sum(r['last_rhs'] is not None and r['last_rhs']<1e-15 for r in blocks),
        dd_terminal_rhs_range=[min(r['last_rhs'] for r in blocks if r['block']>0),
                               max(r['last_rhs'] for r in blocks if r['block']>0)])
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
