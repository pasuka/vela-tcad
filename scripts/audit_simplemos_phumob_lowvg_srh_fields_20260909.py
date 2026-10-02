"""Separate low-Vg SRH rate and integration-volume differences; no response fit."""
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import validate_simplemos_phumob_lowvg_20260909 as s
import analyze_simplemos_phumob_native_cells_20260909 as c
import analyze_simplemos_joint_eight_point_20260907 as support

a, d, OUT, LOCAL = s.a, s.d, s.OUT / 'srh_fields', s.LOCAL


def run():
    a.verify(s.NOUT / 'export_evidence.json')
    a.verify(s.OUT / 'candidate/validation_freeze.json')
    jobs = []
    files = [Path(__file__).resolve(), Path(c.__file__), s.NOUT / 'export_evidence.json', s.OUT / 'candidate/validation_freeze.json']
    for path in (LOCAL / 'candidate/dc/phumob').glob('*/vela/vg_*/attempt_0/result.json'):
        job = a.read(path)
        assert job['qualified']
        jobs.append(job)
        files += [path] + [path.parent / n for n in ('state.csv', 'all_row.csv', 'acceptance_edges.csv', 'config.json')]
    assert len(jobs) == 8
    a.write(OUT / 'srh_contract.json', dict(jobs=jobs,
        scope='Accepted candidate Vela-start states versus native PhuMob fields at the same eight biases; free silicon nodes only.',
        units='SRH rates cm^-3 s^-1, volumes m^2, integrated sources particles per metre per second; charge-equivalent source is q*1e-6 A/um.',
        volume='Production SRH retains all-cell mesh volume. Native comparison uses independently mapped runtime element vertex box measures. No volume or source parameter is changed.',
        interpretation='Separate rate change at the same native box volume from volume change at the Vela rate. This bookkeeping identity is not drain-current attribution or an independently calibrated source response.',
        acceptance_changed=False))
    d.matrix.freeze(OUT / 'srh_freeze.json', files + [OUT / 'srh_contract.json'])
    c.s.LOCAL, c.s.OUT = s.NLOCAL, s.NOUT
    native_points = a.rows(s.NOUT / 'native_points.csv')
    rows, summaries = [], []
    for job in jobs:
        ref = next(r for r in native_points if r['model'] == 'phumob' and r['case'] == job['case'] and int(r['index']) == job['index'])
        data = c.load(ref)
        path = Path(job['dest'])
        geo, mask = support.support(job)
        terms = d.ordered(path / 'all_row.csv', geo.count)
        cfg = a.read(path / 'config.json')
        # No external volume perturbation is enabled in the production run.
        assert not any('srh' in k.lower() and 'volume' in k.lower() for k in cfg['solver'])
        edge_rows = a.rows(path / 'acceptance_edges.csv')
        ratios = [float(e['electron_particle_line_flux_per_m_s']) / float(e['electron_flux'])
                  for e in edge_rows if abs(float(e['electron_flux'])) > 1e-100]
        conversion = ratios[0]
        assert max(abs(r / conversion - 1) for r in ratios) < 1e-12
        native_volume = defaultdict(float)
        for vs in data['cells'].values():
            for v in vs:
                native_volume[data['mapping'][int(v['vertex'])]] += float(v['measure_um2']) * 1e-12
        export = s.NLOCAL / 'exports/phumob' / ref['key']
        manifest = a.read(export / 'field_manifest.json')['fields']
        field = next(r for r in manifest if r['name'] == 'srhRecombination' and r['region'] == 0)
        assert field['unit'] == 'cm^-3*s^-1' and field['support_kind'] == 'node' and field['mapping_status'] == 'complete', field
        rates = c.scalar(export / 'fields/srhRecombination_region0.csv')
        native_total, vela_total, rate_delta, volume_delta = [], [], [], []
        group = []
        for i in np.where(mask)[0]:
            i = int(i)
            vol = float(geo.volumes['all_cell'][i])
            assert vol > 0 and native_volume[i] > 0
            actual = float(terms[i]['electron_recombination']) * conversion / (vol * 1e6)
            native_source = rates[i] * native_volume[i] * 1e6
            vela_source = actual * vol * 1e6
            rate_part = (actual - rates[i]) * native_volume[i] * 1e6
            volume_part = actual * (vol - native_volume[i]) * 1e6
            item = dict(key=ref['key'], node_id=i, native_SRH_cm3_s=rates[i], Vela_SRH_cm3_s=actual,
                rate_absolute_difference_cm3_s=actual - rates[i],
                native_volume_m2=native_volume[i], Vela_volume_m2=vol,
                Vela_over_native_volume=vol / native_volume[i],
                native_source_per_m_s=native_source, Vela_source_per_m_s=vela_source,
                rate_contribution_per_m_s=rate_part, volume_contribution_per_m_s=volume_part)
            rows.append(item)
            group.append(item)
            native_total.append(native_source)
            vela_total.append(vela_source)
            rate_delta.append(rate_part)
            volume_delta.append(volume_part)
        nr, vr, rd, vd = map(math.fsum, (native_total, vela_total, rate_delta, volume_delta))
        closure = abs((vr - nr) - (rd + vd)) / max(abs(nr), abs(vr), 1e-300)
        assert closure < 1e-12
        max_rate = max(group, key=lambda r: abs(r['rate_absolute_difference_cm3_s']))
        summaries.append(dict(key=ref['key'], device=job['device'], vg=job['vg'], vd=job['vd'],
            native_integral_per_m_s=nr, Vela_integral_per_m_s=vr, integrated_source_relative=vr / nr - 1,
            rate_contribution_per_m_s=rd, volume_contribution_per_m_s=vd, bookkeeping_closure=closure,
            max_absolute_rate_difference_cm3_s=abs(max_rate['rate_absolute_difference_cm3_s']), max_rate_difference_node=max_rate['node_id'],
            rate_contribution_charge_equivalent_A_per_um=rd * 1.602176634e-19 * 1e-6,
            volume_contribution_charge_equivalent_A_per_um=vd * 1.602176634e-19 * 1e-6,
            drain_response_calibrated=False))
    a.write_csv(OUT / 'srh_nodes.csv', rows)
    a.write_csv(OUT / 'srh_summary.csv', summaries)
    d.matrix.freeze(OUT / 'srh_evidence.json', [OUT / 'srh_freeze.json', OUT / 'srh_nodes.csv', OUT / 'srh_summary.csv'])
    print(summaries, flush=True)


if __name__ == '__main__':
    run()
