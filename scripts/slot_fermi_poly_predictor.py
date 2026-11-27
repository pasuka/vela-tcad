"""Unqualified rigid-bias predictor for isolated, source-free Fermi PolySilicon.

Shifting psi and both physical QFs by the same bias increment preserves their
differences and thus Fermi occupations at fixed material/doping/temperature.
This is not a Poisson solution or a device acceptance operation.
"""
from pathlib import Path
import hashlib
import math
import logging
import numpy as np
import slot_state_archive as states
from slot_poly_equilibrium import _components

# Advisory monitor only. The continuation policy owns its hard step limit.
GATE_BIAS_WARNING_THRESHOLD_V = .1


def monitor_bias_change(before, target):
    if not math.isfinite(before) or not math.isfinite(target):
        raise ValueError('Nonfinite predictor gate bias')
    span=abs(target-before)
    exceeded=span>GATE_BIAS_WARNING_THRESHOLD_V+1e-12
    if exceeded:
        logging.getLogger(__name__).warning(
            'WARNING gate_predictor_bias_change: |delta_Vg|=%.12g V exceeds advisory %.12g V; '
            'continuing unqualified prediction; Poisson and full Newton acceptance required',
            span,GATE_BIAS_WARNING_THRESHOLD_V)
    return dict(bias_change_V=target-before,warning_threshold_V=GATE_BIAS_WARNING_THRESHOLD_V,
                warning_threshold_exceeded=exceeded,warning_only=True)


def shift_fields(fields, node_deltas):
    if any(name in fields for name in states.archive.PACKED + states.archive.PACKED_LOW):
        raise ValueError('Packed solver coordinates require a separate predictor')
    result = {name: np.array(value, copy=True) for name, value in fields.items()}
    for node, delta in node_deltas.items():
        if not math.isfinite(delta):
            raise ValueError('Nonfinite bias shift')
        for name in ('psi', 'phin', 'phip'):
            result[name][node] += delta
        if states.archive.SPLIT <= fields.keys():
            for carrier, qf in (('electron', 'phin'), ('hole', 'phip')):
                result[carrier+'_qf_reference_V'][node] += delta
                result[qf][node] = (result[carrier+'_qf_reference_V'][node]
                                    + result[carrier+'_qf_increment_V'][node])
    if any(not np.isfinite(value).all() for value in result.values()):
        raise ValueError('Nonfinite predicted field')
    return result


def predict(document, source, output, root, target):
    solver = document['solver']
    if (solver.get('carrier_statistics') != 'fermi_dirac'
            or solver.get('temperature_K', 300) != 300
            or solver.get('mobility', {}).get('constant_materials') != ['PolySilicon']
            or solver.get('recombination_materials') != ['Si']
            or solver.get('impact_ionization', {}).get('current_approximation')
                != 'phumob_element_density_grad_qf'
            or solver.get('band_to_band', 'none') != 'none'
            or document['sweep']['contact'] != 'gate'):
        raise ValueError('Requires the explicit source-free Poly Fermi IdVd profile')
    if any(key in owner for owner in (document, solver) for key in
           ('optical_generation', 'generation', 'thermal', 'electrothermal',
            'electron_quantum_potential', 'btbt')):
        raise ValueError('Unsupported source/quantum/thermal predictor profile')
    fields, meta = states.read(source, document, root)
    if meta.get('state_role') != 'accepted' or meta.get('qualified') is False:
        raise ValueError('Predictor requires a qualified-parent archive')
    for key in ('mesh_file', 'materials_file', 'node_doping_file'):
        actual = hashlib.sha256((Path(root)/document[key]).read_bytes()).hexdigest()
        if meta.get('input_file_sha256', {}).get(key) != actual:
            raise ValueError('Parent input binding mismatch: '+key)
    before = float(meta['contact_biases_V']['gate'])
    bias_monitor=monitor_bias_change(before,target)
    old, components, _, _ = _components(document, Path(root), before)
    new, _, _, _ = _components(document, Path(root), target)
    # Reject an imported nonequilibrium Poly parent rather than resetting QFs.
    if any(abs(float(fields[q][n])-v) > 1e-10
           for n,v in old.items() for q in ('phin', 'phip')):
        raise ValueError('Parent Poly QFs are not at component contact equilibrium')
    deltas = {n: new[n]-v for n,v in old.items()}
    predicted = shift_fields(fields, deltas)
    states.write_candidate(output, predicted, meta, document, [source], target,
                           'fermi_poly_rigid_bias_predictor')
    invariant_error = max(float(np.max(np.abs(
        (predicted['psi']-predicted[q])-(fields['psi']-fields[q]))))
        for q in ('phin', 'phip'))
    return dict(qualified=False, gate_from_V=before, gate_to_V=target,bias_monitor=bias_monitor,
                components=components, shifted_nodes=sum(d != 0 for d in deltas.values()),
                max_psi_minus_qf_roundoff_V=invariant_error,
                densities_bitwise_unchanged=all(np.array_equal(fields[k], predicted[k])
                    for k in ('electrons_m3', 'holes_m3')))
