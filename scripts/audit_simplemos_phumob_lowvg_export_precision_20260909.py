"""Bound sensitivity of native SG port replay to binary64 QF rounding.

This is an export conditioning audit, not a change to solver acceptance or
proof of Sentaurus internal precision. Native DC terminal currents stay intact.
"""
import math
from collections import defaultdict
from pathlib import Path

import validate_simplemos_phumob_lowvg_20260909 as s
import analyze_simplemos_phumob_native_cells_20260909 as c

a, d, OUT = s.a, s.d, s.OUT


def run():
    a.verify(OUT / 'calibration_evidence.json')
    c.s.LOCAL, c.s.OUT = s.NLOCAL, s.NOUT
    a.write(OUT / 'export_precision_contract.json', dict(
        scope='Frozen native densities, electrostatic potential, ni, geometry and native cell mobility; vary only exported electron/hole quasi-Fermi potentials.',
        diagnostic='Aggregate signed drain sensitivity per node before taking absolute values. Estimate half-ULP round-to-nearest uncertainty; verify with one-ULP nextafter perturbations in the maximizing/minimizing gradient directions.',
        limitations='A first-order conditioning estimate under an explicit rounding assumption, not reconstruction of hidden native values or proof of the export cause. No new acceptance threshold and no adjustment to reference current.',
        native_terminal_relative_gate=1e-6))
    d.matrix.freeze(OUT / 'export_precision_freeze.json', [Path(__file__).resolve(), Path(c.__file__),
        OUT / 'calibration_evidence.json', s.NOUT / 'export_evidence.json', OUT / 'export_precision_contract.json'])
    rows = []
    for job in a.rows(s.NOUT / 'native_points.csv'):
        data = c.load(job)
        f = data['fields']
        drain = data['contacts']['drain']
        gradient = defaultdict(float)
        terms = []
        for (i, j), parts in data['parts'].items():
            sign = int(i in drain) - int(j in drain)
            if not sign:
                continue
            coef = {car: math.fsum(g * data['native'][car][cid] for cid, g in parts) for car in ('e', 'h')}
            eta = (f['psi'][j] - f['psi'][i]) / c.VT
            common = -c.Q * 1e-4 * sign
            e = common * -c.VT * f['n'][i] * c.B(-eta - math.log(f['ni'][j] / f['ni'][i])) * coef['e']
            h = common * c.VT * f['p'][i] * c.B(eta + math.log(f['ni'][i] / f['ni'][j])) * coef['h']
            terms.append((i, j, e, h))
            dn = e / c.VT * math.exp((f['fn'][i] - f['fn'][j]) / c.VT)
            dp = h / c.VT * math.exp((f['fp'][j] - f['fp'][i]) / c.VT)
            gradient['fn', i] += dn
            gradient['fn', j] -= dn
            gradient['fp', i] -= dp
            gradient['fp', j] += dp
        def current(values):
            return math.fsum(e * math.expm1((values['fn'][i] - values['fn'][j]) / c.VT) +
                             h * math.expm1((values['fp'][j] - values['fp'][i]) / c.VT)
                             for i, j, e, h in terms)
        base = current(f)
        bound = math.fsum(abs(g) * math.ulp(f[key][node]) / 2 for (key, node), g in gradient.items())
        checks = []
        for sign in (-1, 1):
            perturbed = {key: dict(f[key]) for key in ('fn', 'fp')}
            predicted = []
            for (key, node), g in gradient.items():
                if not g:
                    continue
                value = math.nextafter(f[key][node], math.copysign(math.inf, sign * g))
                perturbed[key][node] = value
                predicted.append(g * (value - f[key][node]))
            predicted = math.fsum(predicted)
            observed = current(perturbed) - base
            checks.append((observed, predicted, abs(observed - predicted) / max(abs(predicted), 1e-300)))
        reference = float(job['Id_A_per_um'])
        mismatch = abs(base - reference)
        rows.append(dict(key=job['key'], model=job['model'], native_Id_A_per_um=reference,
            replay_Id_A_per_um=base, replay_relative=base / reference - 1,
            replay_absolute_mismatch_A_per_um=mismatch,
            half_ulp_QF_bound_A_per_um=bound, half_ulp_bound_over_Id=bound / abs(reference),
            mismatch_over_half_ulp_bound=mismatch / bound if bound else math.inf,
            minus_one_ulp_response_A_per_um=checks[0][0], plus_one_ulp_response_A_per_um=checks[1][0],
            one_ulp_linearization_max_relative=max(x[2] for x in checks),
            mismatch_within_diagnostic_half_ulp_bound=mismatch <= bound,
            replay_qualified=abs(base / reference - 1) <= 1e-6))
    a.write_csv(OUT / 'export_precision.csv', rows)
    summary = dict(states=len(rows), replay_failures=sum(not r['replay_qualified'] for r in rows),
        mismatches_within_half_ulp_bound=sum(r['mismatch_within_diagnostic_half_ulp_bound'] for r in rows),
        max_one_ulp_linearization_relative=max(r['one_ulp_linearization_max_relative'] for r in rows),
        min_half_ulp_bound_over_Id=min(r['half_ulp_bound_over_Id'] for r in rows),
        max_half_ulp_bound_over_Id=max(r['half_ulp_bound_over_Id'] for r in rows),
        export_cause_proven=False, native_algorithm_identified=False, acceptance_changed=False)
    a.write(OUT / 'export_precision_summary.json', summary)
    d.matrix.freeze(OUT / 'export_precision_evidence.json', [OUT / 'export_precision_freeze.json', OUT / 'export_precision.csv', OUT / 'export_precision_summary.json'])
    print(summary, flush=True)


if __name__ == '__main__':
    run()
