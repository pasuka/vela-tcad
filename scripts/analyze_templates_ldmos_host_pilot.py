"""Compare host pilots using verified, relocated exact-point state evidence."""
import argparse
from pathlib import Path

from run_templates_ldmos_linked_d5 import read, write, digest, state_difference


def verified_point(root, point):
    # Exported evidence keeps the relative case layout, not the remote absolute path.
    state = (root/point['case']/'state.csv').resolve()
    if not state.is_relative_to(root.resolve()) or digest(state) != point['sha256']:
        raise ValueError('Missing, unsafe or changed exact state')
    return state


def compare_hosts(local, remote):
    a, b = (read(root/'summary.json') for root in (local, remote))
    if a['status'] != 'pass' or b['status'] != 'pass':
        raise ValueError('Both host pilots must finish successfully')
    for key in ('package_sha256','runner_sha256','schedule','configurations'):
        if a[key] != b[key]: raise ValueError('Host contract mismatch: '+key)
    aa, bb = (read(root/'analysis.json') for root in (local, remote))
    expected = {item['key'] for item in a['schedule']}
    for root,report,analysis in ((local,a,aa),(remote,b,bb)):
        if analysis['source_summary_sha256'] != digest(root/'summary.json'):
            raise ValueError('Stale host analysis')
        if (len(report['cases']) != len(expected)
                or {c['key'] for c in report['cases']} != expected
                or any(c['status'] != 'pass' for c in report['cases'])):
            raise ValueError('Incomplete or unsuccessful host cases')
    stats = [{c['key']:c for c in report['cases']} for report in (aa, bb)]
    rows = []
    for ca in a['cases']:
        cb = next(c for c in b['cases'] if c['key'] == ca['key'])
        da, db = local/ca['name'], remote/cb['name']
        la, lb = (read(d/'fixed/ledger.json') for d in (da,db))
        if len(la['exact_points']) != 8 or len(lb['exact_points']) != 8:
            raise ValueError('Incomplete pilot points')
        points = []
        for x,y in zip(la['exact_points'],lb['exact_points']):
            if x['bias_V'] != y['bias_V']: raise ValueError('Target mismatch')
            points.append(dict(bias_V=x['bias_V'],delta=state_difference(
                verified_point(da,x),verified_point(db,y))))
        maximum = max(p['delta'][k]['max_absolute'] for p in points for k in ('psi','phin','phip'))
        if maximum > 1e-8: raise ValueError('Cross-host state equivalence failed')
        sa,sb = (table[ca['key']] for table in stats)
        row = dict(config=ca['config'],backend=a['configurations'][ca['config']]['backend'],
                   local_wall_s=ca['wall_seconds'],remote_wall_s=cb['wall_seconds'],
                   speed_ratio=ca['wall_seconds']/cb['wall_seconds'],
                   remote_wall_reduction_percent=100*(1-cb['wall_seconds']/ca['wall_seconds']),
                   local_cpu_s=sa['process_tree_cpu_seconds'],remote_cpu_s=sb['process_tree_cpu_seconds'],
                   local_updates=la['total_Newton_updates'],remote_updates=lb['total_Newton_updates'],
                   same_targets=[r['target_V'] for r in la['runs']]==[r['target_V'] for r in lb['runs']],
                   same_updates=[r['Newton_updates'] for r in la['runs']]==[r['Newton_updates'] for r in lb['runs']],
                   max_state_difference_V=maximum,points=points,
                   local_mean_system_busy_percent=sa['mean_system_busy_percent'],
                   remote_mean_system_busy_percent=sb['mean_system_busy_percent'],
                   local_background_cpu_s=sa['background_cpu_seconds_estimate'],
                   remote_background_cpu_s=sb['background_cpu_seconds_estimate'],
                   local_peak_memory_bytes=sa['totals']['peak_solver_process_bytes'],
                   remote_peak_memory_bytes=sb['totals']['peak_solver_process_bytes'],
                   local_stages_s=sa['totals']['seconds'],remote_stages_s=sb['totals']['seconds'])
        rows.append(row)
    return dict(status='pass',scope='Single round, Vg8 first8 D5; not full-curve or repeat qualification',
                package_sha256=a['package_sha256'],runner_sha256=a['runner_sha256'],
                local_host=a['host'],remote_host=b['host'],cases=rows,
                source_hashes={str(root/name):digest(root/name) for root in (local,remote)
                               for name in ('summary.json','analysis.json')})


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('local',type=Path);p.add_argument('remote',type=Path);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();result=compare_hosts(args.local.resolve(),args.remote.resolve())
    write(args.output,result)
    for row in result['cases']:
        print(row['backend'],round(row['local_wall_s'],3),round(row['remote_wall_s'],3),
              'speed',round(row['speed_ratio'],3),'updates',row['local_updates'],row['remote_updates'],
              'max_state_V',row['max_state_difference_V'])
