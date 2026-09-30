"""HFS native formation hypotheses; preserve all failed comparisons.

Boundary projection is specified by the local T-2022.03 manual pp. 452-453.
The 1 V/cm cutoff is an empirical hypothesis, not a documented parameter.
Last-call samples do not identify a unique cell. Cell replay independently
tests every Si cell and retains the existing unfitted PhuMob floor gap.
"""
import argparse
import math
from collections import defaultdict
from pathlib import Path
import numpy as np
import check_simplemos_hfs_native_identity_20260912 as native

a, d, L, O, c = native.a, native.d, native.L, native.O, native.calibration
OUT = O / 'formation_20260914'


def geometry(device):
    g = c.geometry(device)
    raw, xy, cells, _, _, _, _, _ = c.g.geometry(device)
    adj = defaultdict(list)
    for cid, cell in cells.items():
        ns = cell['nodes']
        for k in range(3):
            adj[tuple(sorted((ns[k], ns[(k+1) % 3])))].append(cid)
    boundary = defaultdict(list)
    for edge, ids in adj.items():
        si = [i for i in ids if cells[i]['material'] == 'Si']
        if len(si) == 1 and (len(ids) == 1 or any(cells[i]['material'] != 'Si' for i in ids)):
            u = np.array(xy[edge[1]]) - np.array(xy[edge[0]])
            t = u / np.linalg.norm(u)
            boundary[si[0]].append((edge, np.outer(t, t)))
    return g, set(raw.contact_nodes), boundary


def fields(source):
    return {k: c.scalar(source, v + '_region0.csv') for k, v in
            [('nd', 'DonorConcentration'), ('na', 'AcceptorConcentration'),
             ('n', 'eDensity'), ('p', 'hDensity'), ('psi', 'ElectrostaticPotential'),
             ('phin', 'eQuasiFermiPotential'), ('phip', 'hQuasiFermiPotential')]}


def canali(m, F, b):
    v, beta = ((1.07e7, 1.109), (8.37e6, 1.213))[b]
    return m * math.exp(-math.log1p((m*F/v)**beta)/beta)


def main(stage):
    assert not (OUT / f'{stage}_evidence.json').exists()
    a.verify(O / f'{stage}_observer_evidence.json')
    assert a.read(O / f'{stage}_observer_summary.json')['all_observed_final_states_qualified']
    a.write(OUT / f'{stage}_contract.json', dict(
        source='T-2022.03 manual pp. 448,452-453; original default contact and PartialLayer settings.',
        hypotheses=['raw_qf_any_contact', 'boundary_qf_any_contact', 'boundary_qf_any_contact_cut1',
                    'boundary_qf_edge_contact_cut1'],
        cutoff_status='1 V/cm hard cutoff inferred from pilot samples; not a documented or production-qualified parameter.',
        field_gate_absolute_V_cm=1e-4, field_gate_relative=1e-9,
        cell_gate_relative=1e-7, existing_PhUMob_hole_failure_preserved=True,
        production_changed=False, acceptance_changed=False))
    sample_rows, cell_rows, summaries, coverage = [], [], [], []
    samples = a.rows(O / f'{stage}_observer_samples.csv')
    for job in native.jobs(stage):
        if job['arm'] != 'observed':
            continue
        source = L / f'{stage}_exports' / job['name']
        g, contacts, boundary = geometry(job['device'])
        f = fields(source)
        for b, carrier in enumerate(('e', 'h')):
            bulk = {k: c.bulk({q: f[q][k] for q in ('nd', 'na', 'n', 'p')}, b) for k in f['n']}
            native_mu = c.scalar(source, carrier + 'Mobility_region0_cells.csv', 'cell_id')
            models = {}
            for cid, cell in g['cells'].items():
                ns = cell['nodes']
                gradients = {key: np.array([f[key][ns[1]]-f[key][ns[0]], f[key][ns[2]]-f[key][ns[0]]]) @ g['gradient'][cid]*1e4
                             for key in ('psi', 'phin', 'phip')}
                ef, qf = gradients['psi'], gradients['phin' if b == 0 else 'phip']
                E, Q = float(np.linalg.norm(ef)), float(np.linalg.norm(qf))
                ncontact = sum(k in contacts for k in ns)
                # Retain corner cells explicitly; no arbitrary choice of one edge.
                bs = boundary[cid]
                project = np.eye(2)
                if bs:
                    project = bs[0][1]
                    if any(np.linalg.norm(p-project) > 1e-10 for _, p in bs[1:]):
                        project = np.zeros((2, 2))
                BQ = float(np.linalg.norm(project @ qf))
                fv = {'raw_qf_any_contact': E if ncontact else Q,
                      'boundary_qf_any_contact': E if ncontact else BQ,
                      'boundary_qf_any_contact_cut1': E if ncontact else (BQ if BQ >= 1 else 0),
                      'boundary_qf_edge_contact_cut1': E if ncontact >= 2 else (BQ if BQ >= 1 else 0)}
                en = abs(float(ef @ g['gd'][cid]))
                low = [1/(1/bulk[k] + c.inverse_surface(en, g['distance'][k], f['nd'][k]+f['na'][k], b)) for k in ns]
                w = g['weights'][cid]
                models[cid] = dict(fv=fv, low=low, ncontact=ncontact, boundary_edges=len(bs), raw_qf=Q, boundary_qf=BQ)
                for mode, F in fv.items():
                    mu = math.fsum(wi*canali(mi, F, b) for wi, mi in zip(w, low))
                    meanfirst = canali(math.fsum(wi*mi for wi, mi in zip(w, low)), F, b)
                    cell_rows.append(dict(key=job['key'], carrier=carrier, cell=cid, mode=mode,
                        contact_vertices=ncontact, boundary_edges=len(bs), F_V_cm=F,
                        native_mu=native_mu[cid], replay_mu=mu, relative_error=mu/native_mu[cid]-1,
                        mean_before_saturation_relative=meanfirst/native_mu[cid]-1))
            local = [s for s in samples if s['key'] == job['key'] and s['carrier'] == carrier]
            for s in local:
                k, F = int(s['node']), float(s['F_V_cm'])
                for mode in next(iter(models.values()))['fv']:
                    matches=[]
                    for cid in g['nodecells'][k]:
                        z=models[cid]; predicted=z['fv'][mode]
                        er=abs(predicted-F)
                        if er <= 1e-4+1e-9*abs(F): matches.append(cid)
                    best=min(g['nodecells'][k],key=lambda cid:abs(models[cid]['fv'][mode]-F))
                    z=models[best]
                    sample_rows.append(dict(key=job['key'],carrier=carrier,node=k,mode=mode,
                        observed_F=F,best_field_error_V_cm=abs(z['fv'][mode]-F),
                        matches=len(matches),candidate_cells=';'.join(map(str,matches)),
                        best_cell_diagnostic=best,contact_vertices=z['ncontact'],boundary_edges=z['boundary_edges'],
                        raw_qf=z['raw_qf'],boundary_qf=z['boundary_qf'],
                        uniquely_identified=False))
            for mode in next(iter(models.values()))['fv']:
                rr=[r for r in cell_rows if r['key']==job['key'] and r['carrier']==carrier and r['mode']==mode]
                ss=[r for r in sample_rows if r['key']==job['key'] and r['carrier']==carrier and r['mode']==mode]
                summaries.append(dict(key=job['key'],carrier=carrier,mode=mode,cells=len(rr),samples=len(ss),
                    max_cell_relative=max(abs(r['relative_error']) for r in rr),
                    max_mean_first_relative=max(abs(r['mean_before_saturation_relative']) for r in rr),
                    max_sample_field_V_cm=max(r['best_field_error_V_cm'] for r in ss),
                    samples_without_field_match=sum(r['matches']==0 for r in ss),
                    original_cell_gate_passed=max(abs(r['relative_error']) for r in rr)<=1e-7))
            coverage.append(dict(key=job['key'],carrier=carrier,cells=len(models),
                corners_with_multiple_boundary_edges=sum(x['boundary_edges']>1 for x in models.values()),
                single_contact_vertex_cells=sum(x['ncontact']==1 for x in models.values())))
        print(job['key'], 'formation replay complete', flush=True)
    for name, rows in [('samples',sample_rows),('cells',cell_rows),('summary',summaries),('coverage',coverage)]:
        a.write_csv(OUT / f'{stage}_{name}.csv', rows)
    a.write(OUT/f'{stage}_summary.json',dict(points=len(coverage)//2,cell_comparisons=len(cell_rows),
        sample_comparisons=len(sample_rows),production_changed=False,
        inference_not_native_algorithm_proof=True,full_device_derivatives_qualified=False))
    d.matrix.freeze(OUT/f'{stage}_evidence.json',[Path(__file__).resolve(),Path(c.__file__),
        O/f'{stage}_observer_evidence.json',L.parents[0]/'m79_research/sdevice_ug_local_2022.txt']+
        list(OUT.glob(f'{stage}_*.csv'))+[OUT/f'{stage}_contract.json',OUT/f'{stage}_summary.json'])
    print(summaries,flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=('pilot','rest'),required=True)
    main(p.parse_args().stage)
