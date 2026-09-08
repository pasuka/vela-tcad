"""Read-only comparison of sealed curves and native/strict endpoint fields."""
from decimal import Decimal
from pathlib import Path
import math
import numpy as np
import decompose_simplemos_calibrated_transport as d

a = d.a
OUT = d.ROOT/'physics_jacobian_review_20260906_v2'


def main():
    a.verify(d.FREEZE)
    a.verify(d.OUT/'evidence.json')
    paths = [Path(__file__).resolve(), d.FREEZE, d.OUT/'evidence.json']
    summaries, sources, ranges, nodes = [], [], [], []
    for c in a.read(d.CONTRACT)['cases']:
        key, count = c['case'], c['nodes']
        statepath = d.matrix.spatial.precision.LOCAL/key/'state.csv'
        paths.append(statepath)
        state = d.ordered(statepath, count)
        geo = d.matrix.spatial.m73.Geometry(c['device'])
        masks, _, xy = d.matrix.spatial.old.m78.supports(c['device'], geo, .05)
        fields = {}
        for name in ('ElectrostaticPotential', 'eQuasiFermiPotential', 'hQuasiFermiPotential', 'eDensity', 'hDensity', 'srhRecombination', 'EffectiveIntrinsicDensity'):
            path = Path(c['export'])/'fields'/(name+'_region0.csv')
            paths.append(path)
            fields[name] = {int(r['node_id']):float(r['component0']) for r in a.rows(path)}
        ids = np.array(sorted(fields['eDensity']))
        assert set(ids) == set(np.flatnonzero(masks['all_si']))
        assert all(set(f) == set(ids) for f in fields.values())
        native = {name:np.array([values[i] for i in ids]) for name, values in fields.items()}
        vela = {name:np.array([float(state[i][col]) for i in ids]) for name, col in (
            ('ElectrostaticPotential','psi'), ('eDensity','electrons_m3'), ('hDensity','holes_m3'))}
        vela['eDensity'] *= 1e-6
        vela['hDensity'] *= 1e-6
        for name, carrier in (('eQuasiFermiPotential','electron'), ('hQuasiFermiPotential','hole')):
            vela[name] = np.array([float(Decimal(state[i][carrier+'_qf_reference_V'])+Decimal(state[i][carrier+'_qf_increment_V'])) for i in ids])
        # Re-evaluate the actual configured generalized Boltzmann SRH rate,
        # then check it against the production carrier-term probe.
        terms_path = d.LOCAL/key/'strict/carrier.csv'
        config_path = d.LOCAL/key/'strict/carrier.json'
        edge_path = d.LOCAL/key/'strict/edges.csv'
        paths += [terms_path, config_path, edge_path]
        terms = d.ordered(terms_path, count)
        cfg = a.read(config_path)['solver']['srh_doping_dependence']
        # UnitScalingConfig::concentrationToInternal is identity: despite the
        # legacy *_m3 keys, this deck and carrier probe use internal cm^-3.
        density = np.array([float(terms[i]['donors_m3'])+float(terms[i]['acceptors_m3']) for i in ids])
        lifetimes = {}
        for carrier in ('electron', 'hole'):
            q = cfg[carrier]
            lifetimes[carrier] = q['tau_min_s']+(q['tau_max_s']-q['tau_min_s'])/(1+(density/q['reference_doping_m3'])**q['gamma'])
        ni = np.array([float(terms[i]['ni_eff_m3']) for i in ids])
        split = np.array([float(Decimal(state[i]['hole_qf_reference_V'])+Decimal(state[i]['hole_qf_increment_V'])-Decimal(state[i]['electron_qf_reference_V'])-Decimal(state[i]['electron_qf_increment_V'])) for i in ids])
        rate = ni**2*np.expm1(split/d.fixed.upstream.VT)/(lifetimes['hole']*(vela['eDensity']+ni)+lifetimes['electron']*(vela['hDensity']+ni))
        vela['srhRecombination'] = rate
        # terms are normalized continuity residuals; calibration from production
        # SG edge line flux yields physical charge-current conversion.
        edges = a.rows(edge_path)
        edge = max(edges, key=lambda r:abs(float(r['electron_flux'])))
        factor = float(edge['electron_particle_line_flux_per_m_s'])*d.fixed.Q*1e-6/float(edge['electron_flux'])
        actual = np.array([float(terms[i]['electron_recombination']) for i in ids])*factor
        expected = rate*geo.volumes['all_cell'][ids]*d.fixed.Q
        rate_check = float(np.sum(abs(actual-expected))/np.sum(abs(actual)))
        assert rate_check < 1e-7, rate_check
        area = geo.volumes['barycentric_si'][ids]
        for region in ('all_si', 'channel', 'gate_interface', 'source', 'drain', 'substrate'):
            selected = masks[region][ids]
            weights = area[selected]
            for name in ('ElectrostaticPotential','eQuasiFermiPotential','hQuasiFermiPotential','eDensity','hDensity'):
                nv, vv = native[name][selected], vela[name][selected]
                diff = vv-nv if 'Density' not in name else np.log10(vv/nv)
                peak = int(np.argmax(abs(diff)))
                summaries.append(dict(case=key, region=region, field=name, count=int(sum(selected)),
                    units='V' if 'Density' not in name else 'dex', mean_signed=float(np.dot(weights,diff)/sum(weights)),
                    weighted_rms=float(np.sqrt(np.dot(weights,diff**2)/sum(weights))),
                    max_abs=float(abs(diff[peak])), max_node=int(ids[selected][peak])))
            nr, vr = native['srhRecombination'][selected], rate[selected]
            sources.append(dict(case=key, region=region, units='A/um',
                native_integral=float(d.fixed.Q*np.dot(weights,nr)), vela_integral=float(d.fixed.Q*np.dot(weights,vr)),
                delta_integral=float(d.fixed.Q*np.dot(weights,vr-nr)),
                normalized_L1=float(np.dot(weights,abs(vr-nr))/np.dot(weights,abs(nr))),
                native_max_abs_rate_cm3_s=float(max(abs(nr))), vela_max_abs_rate_cm3_s=float(max(abs(vr))),
                production_rate_reconstruction_L1=rate_check))
        for carrier in ('electron','hole'):
            mu = np.array([float(r[carrier+'_mobility_m2_V_s'])*1e4 for r in edges])
            mu = mu[mu > 0]
            ranges.append(dict(case=key, carrier=carrier, units='cm2/V/s', support='conducting_edges',
                minimum=float(min(mu)), median=float(np.median(mu)), maximum=float(max(mu))))
        for index, node in enumerate(ids):
            nodes.append(dict(case=key,node_id=int(node),x_um=xy[node,0],y_um=xy[node,1],
                **{name+'_native':float(native[name][index]) for name in native},
                **{name+'_vela':float(value[index]) for name,value in vela.items()}))
    curvepath = d.ROOT/'matched_ni_full_curve/m66_matched_ni_full_curve_point_ledger.csv'
    curve_evidence = d.ROOT/'simplemos_m66_matched_ni_full_curve_evidence.json'
    assert a.sha(curvepath) == a.read(curve_evidence)['artifacts'][a.rel(curvepath)]
    paths += [curvepath, curve_evidence]
    curves = a.rows(curvepath)
    bins = []
    for device in ('high_all', 'n23'):
        for vd in (.05, 1.):
            for lo,hi in ((0.,.5),(.55,.9),(.55,1.),(1.,1.),(1.05,2.5)):
                group = [r for r in curves if (r['device'] in ('n21','n22','n23','n24') if device == 'high_all' else r['device'] == device) and float(r['drain_voltage_V']) == vd and lo <= float(r['gate_voltage_V']) <= hi]
                errors = [float(r['absolute_relative_error_percent']) for r in group]
                worst = group[int(np.argmax(errors))]
                bins.append(dict(device=device,vd=vd,vg_low=lo,vg_high=hi,points=len(group),
                    median_absolute_error_percent=float(np.median(errors)), max_absolute_error_percent=max(errors),
                    worst_device=worst['device'],worst_vg=float(worst['gate_voltage_V']),
                    worst_native_A_per_um=float(worst['sentaurus_no_bgn_A_per_um']),
                    worst_vela_A_per_um=float(worst['vela_matched_ni_A_per_um'])))
    for name, rows in (('fields',summaries),('srh',sources),('mobility_ranges',ranges),('nodes',nodes),('curve_bins',bins)):
        a.write_csv(OUT/(name+'.csv'),rows)
    for name in ('src/equation/CoupledDDAssembler.cpp','src/physics/MobilityModel.cpp','src/physics/RecombinationModel.cpp','include/vela/physics/MobilityModel.h','include/vela/solver/NewtonSolver.h','include/vela/core/UnitScaling.h'):
        paths.append(d.REPO/name)
    a.write(OUT/'review.json', dict(input_hashes={a.rel(p):a.sha(p) for p in sorted(set(paths))},
        new_solver_calls=0, sign='Vela minus Sentaurus; density differences are log10(Vela/Sentaurus)',
        rate_definition='Current configured generalized Boltzmann SRH, verified against production carrier terms; common Si barycentric volume for native/strict comparison',
        depth_um=.05, channel_half_width_um=.125, strict_field_points=4,
        full_curve_source='M66 frozen historical matched-ni curve, not a new strict 0..1 V rerun'))
    print('Read-only field, source and curve review complete')


if __name__ == '__main__':
    main()
