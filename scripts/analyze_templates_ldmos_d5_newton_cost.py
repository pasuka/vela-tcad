"""Read-only D5 Newton cost decomposition; never mix with D0 phase counts."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re


def read(p): return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rows(p):
    with Path(p).open(newline='', encoding='utf-8') as f: return list(csv.DictReader(f))


def segment(v):
    # Native log t values are printed at limited precision.
    return 'low_0_9p333' if v <= 28/3+1e-3 else 'mid_9p333_26p667' if v <= 80/3+1e-3 else 'high_26p667_40'


def recovery_counts(curve):
    return {out:sum(int(r.get(field,0)) for r in curve) for out,field in (
        ('internal_recovery_points','carrier_row_recovery_attempted'),
        ('internal_recovery_cycles','carrier_row_recovery_cycles'),
        ('internal_density_passes','carrier_row_recovery_density_passes'))}


def trace_cost(trace, gates, row_gate):
    c = Counter(updates=0); previous = {}; events = Counter()
    for r in trace:
        events[r['event']] += 1
        key = tuple(r[k] for k in ('run_id', 'segment_id', 'attempt_id'))
        before = previous.get(key)
        if r['event'] == 'initial':
            previous[key] = r
            continue
        if r['event'] != 'accepted_iteration': continue
        if before is None: raise ValueError('Accepted iteration has no preceding state')
        c['updates'] += 1
        alpha = float(r['damping'])
        c['damped_updates' if alpha < 1-1e-12 else 'unit_updates'] += 1
        c['tiny_damping_updates'] += alpha < .01
        c['trace_trials'] += int(r['line_search_attempts'])
        c['trace_extra_trials'] += max(0, int(r['line_search_attempts'])-1)
        blocks = all(float(before[field]) <= gates[gate] for field, gate in (
            ('block_psi', 'psi_residual_ceiling'), ('block_phin', 'electron_residual_ceiling'),
            ('block_phip', 'hole_residual_ceiling')))
        c['updates_after_blocks_pass'] += blocks
        c['row_only_tail_updates'] += blocks and float(before['carrier_row_max_ratio']) > row_gate
        c['small_merit_updates'] += float(before['residual_norm']) < 1e-9
        ratio = float(r['residual_norm'])/float(before['residual_norm']) if float(before['residual_norm']) else 0
        c['unit_slow_residual_updates'] += alpha >= 1-1e-12 and ratio > .1
        previous[key] = r
    return dict(c), dict(events)


def analyze_curve(directory):
    ledger = read(directory/'fixed/ledger.json')
    if ledger['status'] != 'completed': raise ValueError('Incomplete curve')
    points = []; totals = Counter(); by_segment = defaultdict(Counter)
    counters = Counter(); by_stage = defaultdict(Counter); guards = Counter()
    hashes = {str(directory/'fixed/ledger.json'): sha(directory/'fixed/ledger.json')}
    for run in ledger['runs']:
        d = directory/run['case']; cfg = read(d/'control.json'); status = read(d/'status.json')
        trace = rows(d/'newton_iterations.csv'); profile = read(d/'performance_profile.json')
        cost, events = trace_cost(trace, cfg['solver']['block_absolute_convergence'],
                                  cfg['solver']['carrier_row_convergence']['eps_row'])
        if cost.get('updates', 0) != run['Newton_updates'] or cost.get('updates', 0) != profile['counters'].get('newton.updates', 0):
            raise ValueError('Update reconciliation failed: '+str(d))
        attempts = status['attempts']
        if sum(int(a['newton_iterations']) for a in attempts) != cost.get('updates', 0):
            raise ValueError('Attempt count mismatch')
        cost.update(services=1, wall_seconds=run['wall_seconds'],
            failed_attempts=sum(a['status'] != 'accepted' for a in attempts),
            updates_in_failed_attempts=sum(int(a['newton_iterations']) for a in attempts if a['status'] != 'accepted'),
            predicted_services=int(bool(run.get('outer_predictor'))),**recovery_counts(status['curve']))
        initial = next((float(r['residual_norm']) for r in trace if r['event'] == 'initial'), None)
        point = dict(case=run['case'], target_V=run['target_V'], parent_V=run['parent_V'],
            step_V=run['target_V']-run['parent_V'], stage=run['stage'],
            predicted=bool(run.get('outer_predictor')), guard=run.get('predictor_guard_reason'),
            initial_residual=initial, events=events, **cost)
        points.append(point); totals.update(cost); by_stage[run['stage']].update(cost)
        by_segment[segment(run['target_V'])].update(cost); counters.update(profile['counters'])
        guards[run.get('predictor_guard_reason') or 'used'] += 1
        for name in ('control.json', 'newton_iterations.csv', 'status.json', 'performance_profile.json'):
            hashes[str(d/name)] = sha(d/name)
    if totals['updates'] != ledger['total_Newton_updates']: raise ValueError('Ledger count mismatch')
    transfers = ledger['transfers']
    return dict(directory=str(directory), totals=dict(totals), segments=dict(by_segment),
        stages=dict(by_stage), profile_counters=dict(counters), prediction_guards=dict(guards),
        advances=len(transfers), rollbacks=len(ledger['rollbacks']),
        max_step=max(t['accepted_step_V'] for t in transfers),
        points=points, hashes=hashes,
        top_updates=sorted(points, key=lambda r:r['updates'], reverse=True)[:12],
        top_row_tail=sorted(points, key=lambda r:r.get('row_only_tail_updates',0), reverse=True)[:12])


def native_tables(text):
    found=[]
    pattern=r'Iteration\s+\|Rhs\|[^\n]*\n-+\n(.*?)Finished, because\.\.\.\n([^\n]+)'
    for m in re.finditer(pattern, text, re.S):
        positive=[]
        for line in m[1].splitlines():
            p=line.split()
            if len(p)==8 and p[0].isdigit() and float(p[2])>0: positive.append(p)
        found.append(dict(updates=len(positive), damped=sum(float(p[2])<1 for p in positive),
            stop=m[2].strip(), success=('smaller than' in m[2].lower())))
    return found


def analyze_native(path):
    text=Path(path).read_text(errors='replace')
    parts=re.split(r'(?m)^Quasistationary \(', text)
    drains=[p for p in parts if 'Contact drain : 40V' in p and 'Computing step from' in p]
    if len(drains)!=2: raise ValueError('Expected two executed drain sweeps')
    result={}
    for gate, part in zip((4,8),drains):
        steps=list(re.finditer(r'Computing step from t=(\S+) to t=(\S+) \(Stepsize: (\S+)\)',part))
        records=[]
        for i,m in enumerate(steps):
            body=part[m.end():steps[i+1].start() if i+1<len(steps) else len(part)]
            tables=native_tables(body)
            if not tables: raise ValueError('Missing native step table')
            # The last block can also contain the next gate's loaded zero solve.
            row=dict(parent_V=40*float(m[1]),target_V=40*float(m[2]),step_V=40*float(m[3]),
                predicted='Extrapolating values' in body[:body.find('Iteration')], **tables[0])
            times=re.search(r'Total time:\s*([\d.]+) s',body)
            row['logged_solve_seconds']=float(times[1]) if times else None
            records.append(row)
        if abs(records[-1]['target_V']-40)>1e-6: raise ValueError('Native sweep incomplete')
        bins=defaultdict(Counter)
        for r in records:
            bins[segment(r['target_V'])].update(advances=1,updates=r['updates'],damped=r['damped'],
                failed_attempts=int(not r['success']),logged_solve_seconds=r['logged_solve_seconds'] or 0)
        result[str(gate)]=dict(advances=len(records),updates=sum(r['updates'] for r in records),
            failed_attempts=sum(not r['success'] for r in records),max_step=max(r['step_V'] for r in records),
            zero_solve_tables=native_tables(part[:steps[0].start()]),segments=dict(bins),points=records)
    return dict(path=str(path),sha256=sha(path),gates=result,
        scope='Executed drain steps only; zero solves separate; prebias excluded; logged solver time is not wall time')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch',type=Path,required=True);p.add_argument('--native-log',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); a.output.mkdir(parents=True,exist_ok=False)
    batch=read(a.batch/'summary.json')
    if batch['status']!='pass': raise ValueError('Timing batch incomplete')
    runs=[]
    for case in batch['cases']:
        if case['config']!='U1': continue
        report=analyze_curve(a.batch/case['name'])
        report.update(gate=case['gate'],round=case['round'],curve_wall_seconds=case['wall_seconds'])
        runs.append(report)
    if len(runs)!=6: raise ValueError('Expected all six UMFPACK curves')
    report=dict(schema='vela.d5_newton_cost.v1',runs=runs,native=analyze_native(a.native_log),
        runner_sha256=batch['runner_sha256'],batch_sha256=sha(a.batch/'summary.json'),
        analyzer_sha256=sha(Path(__file__)),
        limitations=['D5 only, no D0 phase-count transfer',
            'Trace damping has no initial-alpha/component/node limiter attribution',
            'Row-only tail uses original gates, not a proposed acceptance change',
            'Native and Vela stop norms/timing environments differ; do not compare raw residual magnitudes'])
    (a.output/'analysis.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    for r in runs: print(r['gate'],r['round'],r['advances'],r['totals'],r['prediction_guards'],flush=True)


if __name__=='__main__': main()
