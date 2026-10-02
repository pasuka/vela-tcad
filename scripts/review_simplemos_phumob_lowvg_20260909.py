"""Review 16-point scope, minority rows, source support and native fields."""
import argparse
import ast
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import validate_simplemos_phumob_lowvg_20260909 as s
import analyze_simplemos_joint_eight_point_20260907 as rowsupport

a, d, OUT, LOCAL = s.a, s.d, s.OUT, s.LOCAL


def post():
    """Qualify derivatives of accepted states independently of A/B completion."""
    a.verify(OUT / 'candidate/validation_freeze.json')
    attempts = []
    for path in (LOCAL / 'candidate/dc/phumob').glob('*/vela/vg_*/attempt_*/result.json'):
        r = a.read(path)
        if r['qualified']:
            a.verify(path.parent / 'input_freeze.json')
            attempts.append(dict(r, qualified='True'))
    jobs, files = [], []
    for case, index in sorted({(r['case'], r['index']) for r in attempts}):
        selected = [r for r in attempts if r['case'] == case and r['index'] == index and
                    r['arm'] == 'vela' and r['qualified'] == 'True']
        assert selected, (case, index)
        r = selected[0]
        src = Path(r['dest'])
        cfg = a.read(src / 'config.json')
        key = f'{case}_vg_{int(index):03d}'
        dest = LOCAL / 'post_jvp' / key
        _, mask = rowsupport.support(r)
        cfg.update(simulation_type='newton_jvp_probe', state_file=str(src / 'state.csv'), output_csv=str(dest / 'jvp.csv'))
        cfg.pop('output_state_file', None)
        cfg['directions'] = [dict(name=f'{mode}_{step:.0e}', mode=mode, amplitude_V=step,
            node_ids=np.where(mask)[0][::3].tolist(), exclude_contacts=True)
            for mode in ('psi', 'phin', 'phip') for step in (1e-4, 3e-5, 1e-5, 3e-6)]
        a.write(dest / 'config.json', cfg)
        jobs.append(dict(key=key, path=str(dest / 'config.json')))
        files += [dest / 'config.json', src / 'state.csv', src / 'result.json', src / 'independent_acceptance.json']
    assert len(jobs) == 8
    a.write(OUT / 'post_contract.json', dict(jobs=jobs, gate=1e-4,
        scope='Eight accepted Vela-initialized candidate states; this derivative gate is independent of native-initialization completion and does not grant dual-initialization acceptance.',
        weak_SRH_cross='Retained in output; excluded from whole-residual finite-difference qualification.'))
    d.matrix.freeze(OUT / 'post_freeze.json', files + [OUT / 'post_contract.json', OUT / 'candidate/validation_freeze.json', s.q.run.RUNNER, Path(__file__).resolve()])
    def one(job):
        path = Path(job['path'])
        status = s.q.run.V.execute(path, s.q.run.RUNNER, s.q.run.V.environment())
        assert status['exit_code'] == 0, status
        rows = []
        for r in a.rows(path.parent / 'jvp.csv'):
            for block in ('psi', 'phin', 'phip'):
                ana, fd = float(r[f'analytic_{block}_norm']), float(r[f'finite_difference_{block}_norm'])
                error = float(r[f'{block}_relative_error']) * max(1., fd) / max(ana, fd, 1e-300)
                weak = (r['mode'], block) in (('phin', 'phip'), ('phip', 'phin'))
                rows.append(dict(key=job['key'], input=r['mode'], output=block, step_V=r['amplitude_V'],
                    relative=error, gated=not weak and float(r['amplitude_V']) < 1e-4, qualified=error <= 1e-4))
        print('Post Jv', job['key'], flush=True)
        return rows
    with ThreadPoolExecutor(max_workers=2) as pool:
        rows = [r for group in pool.map(one, jobs) for r in group]
    a.write_csv(OUT / 'post_jvp.csv', rows)
    gated = [r for r in rows if r['gated']]
    summary = dict(states=len(jobs), gated_checks=len(gated), gated_failures=sum(not r['qualified'] for r in gated), max_relative=max(r['relative'] for r in gated))
    a.write(OUT / 'post_summary.json', summary)
    d.matrix.freeze(OUT / 'post_evidence.json', [OUT / 'post_freeze.json', OUT / 'post_jvp.csv', OUT / 'post_summary.json'] + [p for p in (LOCAL / 'post_jvp').rglob('*') if p.is_file()])
    print(summary, flush=True)


def review():
    for p in (OUT / 'dc_comparison_evidence.json', OUT / 'post_evidence.json',
              s.NOUT / 'export_evidence.json', OUT / 'calibration_evidence.json', OUT / 'export_precision_evidence.json'):
        a.verify(p)
    original_jobs = a.read(s.q.OUT / 'validation_contract.json')['jobs']
    for job in a.read(OUT / 'contract.json')['jobs']:
        original = next(j for j in original_jobs if j['model'] == 'old_slotboom' and
                        (j['case'], j['index'], j['arm']) == (job['case'], job['index'], job['arm']))
        expected = a.read(Path(original['original_config']))
        expected['solver']['mobility'].update(model='phumob', edge_averaging='element_box_phumob')
        assert expected == a.read(Path(job['original_config']))
        assert job['seed'] == original['seed']
    assert a.read(OUT / 'contract.json')['gates'] == a.read(s.previous.OUT / 'contract.json')['gates']
    combined = a.rows(s.previous.OUT / 'dc_comparison.csv') + a.rows(OUT / 'dc_comparison.csv')
    assert len(combined) == len({r['key'] for r in combined}) == 16
    pairs = []
    for vd, vg in sorted({(r['vd'], r['vg']) for r in combined}):
        lo = next(r for r in combined if r['vd'] == vd and r['vg'] == vg and r['device'] == 'n19')
        hi = next(r for r in combined if r['vd'] == vd and r['vg'] == vg and r['device'] == 'n23')
        for arm in ('legacy', 'candidate'):
            nr = float(hi['native_Id_A_per_um']) / float(lo['native_Id_A_per_um'])
            vr = float(hi[arm + '_Id_A_per_um']) / float(lo[arm + '_Id_A_per_um'])
            pairs.append(dict(vd=vd, vg=vg, stage=arm, native_high_over_low=nr, Vela_high_over_low=vr,
                              ratio_error_percent=100 * (vr / nr - 1),
                              qualified=hi[arm + '_dual_qualified'] == lo[arm + '_dual_qualified'] == 'True'))
    ranks, fieldrows, sources, attempts = [], [], [], []
    for origin, root in [('low', OUT), ('high', s.previous.OUT)]:
        for arm in ('legacy', 'candidate'):
            for attempt in a.rows(root / arm / 'attempts.csv'):
                attempts.append(dict(origin=origin, stage=arm, **attempt))
        for r in a.rows(root / 'selected_states.csv'):
            path = Path(r['dest'])
            geo, mask = rowsupport.support(r)
            state = d.ordered(path / 'state.csv', geo.count)
            terms = d.ordered(path / 'all_row.csv', geo.count)
            for carrier, density in [('electron', 'electrons_m3'), ('hole', 'holes_m3')]:
                rr = []
                for i in np.where(mask)[0]:
                    t = terms[i]
                    scale = max(float(t[carrier + '_flux_abs_sum']), abs(float(t[carrier + '_recombination'])), abs(float(t[carrier + '_impact'])))
                    ratio = abs(float(t[carrier + '_residual'])) / scale if scale else math.inf
                    other = 'holes_m3' if carrier == 'electron' else 'electrons_m3'
                    rr.append(dict(node_id=int(i), row_ratio=ratio, carrier_density_m3=float(state[i][density]),
                        is_local_minority=float(state[i][density]) < float(state[i][other]),
                        residual=float(t[carrier + '_residual']), edge_abs_sum=float(t[carrier + '_flux_abs_sum']),
                        srh_row=float(t[carrier + '_recombination'])))
                ranked = sorted(rr, key=lambda x: x['row_ratio'], reverse=True)
                for rank, item in enumerate(ranked[:5]):
                    ranks.append(dict(origin=origin, key=r['key'], stage=r['stage'], initialization=r['arm'],
                                      qualified=r['qualified'], carrier=carrier, rank=rank + 1, **item))
                source = a.read(path / 'independent_acceptance.json')['closure'][carrier]
                sources.append(dict(origin=origin, key=r['key'], stage=r['stage'], initialization=r['arm'], carrier=carrier,
                    absolute_integrated_source=abs(source['integrated_source']), source_floor_active=source['qualified'],
                    unfloored_relative=abs(source['contact_flux'] - source['integrated_source']) / max(abs(source['contact_flux']), abs(source['integrated_source']), 1e-300),
                    frozen_floor_ratio=source['ratio'], accepted_under_frozen_gate=source['satisfied']))
            if origin != 'low' or r['arm'] != 'vela':
                continue
            export = s.NLOCAL / 'exports/phumob' / r['key']
            mapping = s.q.run.w.mapping
            psi, en, ho, spread = mapping.RAW_READER(export, geo)
            assert spread < 1e-12
            fn = mapping.original.m73.scalar(export / 'fields/eQuasiFermiPotential_region0.csv')
            fp = mapping.original.m73.scalar(export / 'fields/hQuasiFermiPotential_region0.csv')
            native = {'psi': psi, 'phin': fn, 'phip': fp, 'electrons_m3': en, 'holes_m3': ho}
            for name, values in native.items():
                ids = [int(i) for i in np.where(mask)[0]]
                actual = [float(rowsupport.physical(state[i], name)) if name in ('psi', 'phin', 'phip') else float(state[i][name]) for i in ids]
                delta = np.array([x - values[i] for x, i in zip(actual, ids)])
                k = int(np.argmax(np.abs(delta)))
                fieldrows.append(dict(key=r['key'], stage=r['stage'], field=name, qualified=r['qualified'],
                    max_absolute=float(np.max(np.abs(delta))), rms=float(np.sqrt(np.mean(delta * delta))), worst_node=ids[k],
                    max_relative=max(abs(x / values[i] - 1) for x, i in zip(actual, ids)) if name.endswith('_m3') else '',
                    units='m^-3' if name.endswith('_m3') else 'V'))
    for name, rr in [('sixteen_points.csv', combined), ('nwell_pairs.csv', pairs), ('minority_row_ranking.csv', ranks),
                     ('native_field_comparison.csv', fieldrows), ('source_support.csv', sources)]:
        a.write_csv(OUT / name, rr)
    s.q.run.v.csv_union(OUT / 'all_attempts.csv', attempts)
    summary = dict(points=16, low_points=8, low_attempts=sum(r['origin'] == 'low' for r in attempts),
        low_failed_attempts=sum(r['origin'] == 'low' and r['qualified'] != 'True' for r in attempts),
        qualified_points={arm: sum(r[arm + '_dual_qualified'] == 'True' for r in combined) for arm in ('legacy', 'candidate')},
        max_absolute_Id_error_percent={arm: max(abs(float(r[arm + '_versus_native_percent'])) for r in combined) for arm in ('legacy', 'candidate')},
        improved_points=sum(abs(float(r['candidate_versus_native_percent'])) < abs(float(r['legacy_versus_native_percent'])) for r in combined),
        source_components=len(sources), source_floor_active=sum(r['source_floor_active'] for r in sources),
        unfloored_source_max_relative=max(r['unfloored_relative'] for r in sources),
        candidate_worst_row=max(r['row_ratio'] for r in ranks if r['stage'] == 'candidate'),
        native_G_algorithm_identified=False, full_curve=False, defaults_changed=False)
    a.write(OUT / 'review_summary.json', summary)
    d.matrix.freeze(OUT / 'review_evidence.json', [Path(__file__).resolve(), OUT / 'dc_comparison_evidence.json',
        OUT / 'post_evidence.json', OUT / 'calibration_evidence.json', OUT / 'export_precision_evidence.json', s.NOUT / 'export_evidence.json',
        s.previous.OUT / 'dc_comparison_evidence.json'] + [OUT / n for n in
        ('sixteen_points.csv', 'nwell_pairs.csv', 'minority_row_ranking.csv', 'native_field_comparison.csv', 'source_support.csv', 'all_attempts.csv', 'review_summary.json')])
    print(summary, flush=True)


def seal():
    import re
    report = s.REPO / 'docs/validation/simplemos_phumob_lowvg_validation_2026-09-09.md'
    scripts = [s.REPO / 'scripts' / n for n in ('validate_simplemos_phumob_lowvg_20260909.py',
        'analyze_simplemos_phumob_lowvg_20260909.py', 'review_simplemos_phumob_lowvg_20260909.py',
        'audit_simplemos_phumob_lowvg_export_precision_20260909.py')]
    for script in scripts:
        ast.parse(script.read_text(encoding='utf8'))
    manifests, identities = [], {}
    for p in OUT.rglob('*.json'):
        obj = a.read(p)
        if not isinstance(obj, dict) or 'input_hashes' not in obj:
            continue
        manifests.append(p)
        for name, value in obj['input_hashes'].items():
            assert name not in identities or identities[name] == value, name
            identities[name] = value
    for name, expected in identities.items():
        assert a.sha(s.REPO / name) == expected, name
    links = re.findall(r'\]\(([^)]+)\)', report.read_text(encoding='utf8'))
    for link in links:
        assert (report.parent / link).resolve().exists(), link
    record = dict(manifests=len(manifests), file_hashes=len(identities), syntax_checks=len(scripts), links=len(links),
                  production_source_changed=False, new_native_simulations=16, full_ctest_repeated=False)
    a.write(OUT / 'completion_review.json', record)
    d.matrix.freeze(OUT / 'completion_evidence.json', manifests + scripts + [report, OUT / 'completion_review.json',
        s.REPO / 'docs/validation/simplemos_branch_status.md'])
    print(record, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('post', 'review', 'seal'))
    globals()[parser.parse_args().action]()
