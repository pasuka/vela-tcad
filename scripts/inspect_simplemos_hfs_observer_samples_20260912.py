"""Match actual HFS observer samples to final states and adjacent cell hypotheses.

No nearest-cell choice is promoted to a unique identity. Every adjacent cell
and both candidate fields are retained; numerical errors are diagnostic.
"""
import argparse
import math
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import check_simplemos_hfs_native_identity_20260912 as identity

a, d, L, O, calibration = identity.a, identity.d, identity.L, identity.O, identity.calibration


def main(stage):
    a.verify(O / f'{stage}_identity_evidence.json')
    assert a.read(O / f'{stage}_identity_summary.json')['all_identity_qualified']
    samples, candidates, summaries = [], [], []
    gates = a.read(O / 'contract.json')['gates']
    for job in identity.jobs(stage):
        if job['arm'] != 'observed':
            continue
        source = L / f'{stage}_exports' / job['name']
        geometry = calibration.geometry(job['device'])
        raw_geometry = calibration.g.geometry(job['device'])[0]
        contacts = set(raw_geometry.contact_nodes)
        xy = {int(r['id']): (float(r['x_um']), float(r['y_um']))
              for r in a.rows(source / 'nodes.csv')}
        assert geometry['xy'] == xy
        nodes = sorted(xy)
        tree = cKDTree([xy[k] for k in nodes])
        fields = {key: calibration.scalar(source, field + '_region0.csv')
                  for key, field in [('psi', 'ElectrostaticPotential'),
                                     ('phin', 'eQuasiFermiPotential'),
                                     ('phip', 'hQuasiFermiPotential'),
                                     ('n', 'eDensity'), ('p', 'hDensity'),
                                     ('nd', 'DonorConcentration'),
                                     ('na', 'AcceptorConcentration')]}
        cell_fields = {}
        for cid, cell in geometry['cells'].items():
            ns = cell['nodes']
            gradients = {}
            for field in ('psi', 'phin', 'phip'):
                value = fields[field]
                gradients[field] = np.array([value[ns[1]] - value[ns[0]],
                                             value[ns[2]] - value[ns[0]]]) @ geometry['gradient'][cid] * 1e4
            cell_fields[cid] = dict(E=float(np.linalg.norm(gradients['psi'])),
                                    Fn=float(np.linalg.norm(gradients['phin'])),
                                    Fp=float(np.linalg.norm(gradients['phip'])),
                                    Enormal=abs(float(gradients['psi'] @ geometry['gd'][cid])),
                                    contact_vertices=sum(k in contacts for k in ns))
        local = []
        for b, carrier in enumerate(('e', 'h')):
            rows = a.rows(L / f'{stage}_raw' / 'bundle' / job['name'] /
                          f'hfs_observer_{carrier}_Silicon_1.csv')
            seen = set()
            for row in rows:
                distance, index = tree.query([float(row['x_um']), float(row['y_um'])])
                node = nodes[int(index)]
                assert distance <= 1e-12 and node in fields['n'] and node not in seen
                seen.add(node)
                temperature = float(row['temperature_K'])
                assert temperature == 300 and float(row['carrier_temperature_K']) == 300
                potential = abs(float(row['potential_V']) - fields['psi'][node])
                nr = abs(float(row['n_cm3']) / fields['n'][node] - 1)
                pr = abs(float(row['p_cm3']) / fields['p'][node] - 1)
                F, low, mu = (float(row[k]) for k in ('F_V_cm', 'mulow_cm2_V_s', 'mu_cm2_V_s'))
                assert F >= 0 and low >= 0 and all(math.isfinite(x) for x in (F, low, mu))
                vsat, beta = ((1.07e7, 1.109), (8.37e6, 1.213))[b]
                expected = low / (1 + (low * F / vsat) ** beta) ** (1 / beta)
                qualified_state = (potential <= gates['observer_psi_V'] and
                                   max(nr, pr) <= gates['observer_density_relative'])
                sample = dict(key=job['key'], carrier=carrier, node=node,
                              potential_error_V=potential, n_relative=nr, p_relative=pr,
                              final_state_qualified=qualified_state, F_V_cm=F,
                              mulow_cm2_V_s=low, mu_cm2_V_s=mu,
                              canali_mu_relative=abs(expected / mu - 1) if mu else abs(expected),
                              driving_force_enum=int(row['driving_force_enum']),
                              interface_distance_error_um=abs(float(row['interface_distance_um']) - geometry['distance'][node]),
                              adjacent_cells=len(geometry['nodecells'][node]))
                samples.append(sample)
                local.append(sample)
                bulk = calibration.bulk({k: fields[k][node] for k in ('nd', 'na', 'n', 'p')}, b)
                for cid in geometry['nodecells'][node]:
                    cell = cell_fields[cid]
                    fq = cell['Fn' if carrier == 'e' else 'Fp']
                    low_prediction = 1 / (1 / bulk + calibration.inverse_surface(
                        cell['Enormal'], geometry['distance'][node],
                        fields['nd'][node] + fields['na'][node], b))
                    candidates.append(dict(key=job['key'], carrier=carrier, node=node, cell=cid,
                                           **cell, observed_F=F,
                                           E_error_V_cm=abs(cell['E'] - F),
                                           QF_error_V_cm=abs(fq - F),
                                           any_contact_switch_error_V_cm=abs((cell['E'] if cell['contact_vertices'] else fq) - F),
                                           contact_edge_switch_error_V_cm=abs((cell['E'] if cell['contact_vertices'] >= 2 else fq) - F),
                                           mathematical_floor_mulow=low_prediction,
                                           mulow_relative=low_prediction / low - 1 if low else None))
            assert seen
            summaries.append(dict(key=job['key'], carrier=carrier, samples=len(seen),
                                  silicon_nodes=len(fields['n']), missing_nodes=len(set(fields['n']) - seen)))
        assert local
    assert samples
    for name, rows in [('samples', samples), ('cell_hypotheses', candidates), ('coverage', summaries)]:
        a.write_csv(O / f'{stage}_observer_{name}.csv', rows)
    summary = dict(samples=len(samples), final_state_qualified=sum(r['final_state_qualified'] for r in samples),
                   all_observed_final_states_qualified=all(r['final_state_qualified'] for r in samples),
                   canali_sample_max_relative=max(r['canali_mu_relative'] for r in samples),
                   cell_formation_qualified=False, production_changed=False,
                   scope='All adjacent hypotheses retained. Coverage is last-call only. Mathematical PhuMob floor is unfitted and retains its existing hole qualification gap. No field or cell identity is inferred from a minimum error alone.')
    a.write(O / f'{stage}_observer_summary.json', summary)
    d.matrix.freeze(O / f'{stage}_observer_evidence.json',
                    [Path(__file__).resolve(), O / f'{stage}_identity_evidence.json'] +
                    [O / f'{stage}_observer_{n}' for n in
                     ('samples.csv', 'cell_hypotheses.csv', 'coverage.csv', 'summary.json')])
    print(summary, flush=True)
    assert summary['all_observed_final_states_qualified']


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', choices=('pilot', 'rest'), required=True)
    main(parser.parse_args().stage)
