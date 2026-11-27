"""Explicit electrical-frame candidates for the fixed SLOT Si/Poly profile.

Coordinates and metadata are translated together. This grants no convergence
qualification, even if the parent archive was accepted in another frame.
"""
import copy
import json
import math
from pathlib import Path

import slot_state_archive as states


def finite(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError("reference-frame voltages must be finite numbers")
    return float(value)


def translate_document(document, target_origin):
    target_origin = finite(target_origin)
    old = finite(document.get('potential_origin_V', 0.))
    result = copy.deepcopy(document)
    circuit = result.get('sweep', {}).get('external_circuit', {})
    if circuit.get('enabled', bool(circuit)) and circuit.get('mode', 'none') != 'none':
        raise ValueError('reference-frame continuation requires direct device bias')
    if result.get('boundaries'):
        raise ValueError('explicit boundaries require a separate frame contract')
    contacts = result.get('contacts', [])
    names = [c['name'] for c in contacts]
    if len(set(names)) != len(names) or not {'source', 'drain'} <= set(names):
        raise ValueError('requires unique SLOT source/drain contacts')
    if any(c.get('type') != 'ohmic' for c in contacts):
        raise ValueError('only the validated SLOT ohmic contact profile is supported')
    if result['sweep'].get('contact') != 'drain':
        raise ValueError('requires drain-driven SLOT sweep')
    for c in contacts:
        c['bias'] = math.fsum((finite(c.get('bias', 0.)), old, -target_origin))
    sweep = result['sweep']
    for key in ('start', 'stop'):
        if key in sweep:
            sweep[key] = math.fsum((finite(sweep[key]), old, -target_origin))
    if 'bias_points' in sweep:
        sweep['bias_points'] = [math.fsum((finite(v), old, -target_origin)) for v in sweep['bias_points']]
    # Resource relocation is done separately. Do not carry a stale seed/output.
    result.pop('state_file', None)
    for key in ('initial_state_file', 'write_state_file', 'write_state_every_point_prefix'):
        sweep.pop(key, None)
    result['potential_origin_V'] = target_origin
    result['state_format'] = 'hdf5'
    return result


def translate_candidate(source, destination, source_document, directory,
                        target_document, physical_drain_voltage):
    source_origin = finite(source_document.get('potential_origin_V', 0.))
    target_origin = finite(target_document.get('potential_origin_V', 0.))
    physical_drain_voltage = finite(physical_drain_voltage)
    fields, metadata = states.read(source, source_document, directory)
    count, identity, _ = states.identity(source_document, directory)
    target_count, target_identity, _ = states.identity(target_document, directory)
    if (count, identity) != (target_count, target_identity):
        raise ValueError('frame conversion cannot change mesh identity')
    mesh = json.loads((Path(directory) / source_document['mesh_file']).read_text(encoding='utf-8-sig'))
    allowed = {'Si', 'PolySilicon', 'SiO2', 'Oxide', 'Nitride'}
    if any(r['material'] not in allowed for r in mesh['regions']):
        raise ValueError('unsupported SLOT material activity classification')
    region_ids = {r['id'] for r in mesh['regions'] if r['material'] in {'Si', 'PolySilicon'}}
    active = {i for t in mesh['triangles'] if t['region_id'] in region_ids for i in t['node_ids']}
    target_biases = {c['name']: finite(c.get('bias', 0.)) for c in target_document['contacts']}
    target_biases['drain'] = physical_drain_voltage - target_origin
    old_biases = metadata.get('contact_biases_V', {})
    if set(old_biases) != set(target_biases) or metadata.get('bias_contact') != 'drain':
        raise ValueError('seed requires complete drain/contact bias provenance')
    for name, value in old_biases.items():
        physical = math.fsum((finite(value), source_origin))
        expected = math.fsum((target_biases[name], target_origin))
        if not math.isclose(physical, expected, rel_tol=0., abs_tol=1e-12):
            raise ValueError('seed terminal voltage does not match declared physical bias')
    if not math.isclose(finite(metadata.get('bias_V')), finite(old_biases['drain']), rel_tol=0., abs_tol=1e-12):
        raise ValueError('inconsistent seed drain bias metadata')
    delta = target_origin - source_origin
    inactive_contacts = {}
    for c in mesh['contacts']:
        if c['name'] not in target_biases:
            raise ValueError('mesh contact missing from frame configuration')
        for i in c['node_ids']:
            if i not in active:
                value = target_biases[c['name']]
                if i in inactive_contacts and inactive_contacts[i] != value:
                    raise ValueError('conflicting inactive contact constraints')
                inactive_contacts[i] = value
    if set(c['name'] for c in mesh['contacts']) != set(target_biases):
        raise ValueError('configuration/mesh contacts do not match')
    for i in range(count):
        fields['psi'][i] = math.fsum((float(fields['psi'][i]), -delta))
        live = i in active
        if not live and (fields['electrons_m3'][i] != 0 or fields['holes_m3'][i] != 0):
            raise ValueError('nontransport node has nonzero carrier density')
        for carrier, key in (('electron', 'phin'), ('hole', 'phip')):
            ref, inc = carrier+'_qf_reference_V', carrier+'_qf_increment_V'
            if ref not in fields or inc not in fields:
                raise ValueError('SLOT frame conversion requires split QF fields')
            a, b = float(fields[ref][i]), float(fields[inc][i])
            if live:
                new = math.fsum((a, -delta))
                fields[ref][i] = new
                fields[inc][i] = math.fsum((a, -delta, -new, b))
                fields[key][i] = math.fsum((new, float(fields[inc][i])))
            else:
                # Insulator contact dummy QF equations follow terminal bias;
                # interior inactive equations use zero. Neither is a density.
                fields[key][i] = fields[ref][i] = inactive_contacts.get(i, 0.)
                fields[inc][i] = 0.
    metadata['potential_origin_V'] = target_origin
    metadata['frame_translation'] = dict(source_origin_V=source_origin,
        target_origin_V=target_origin, inactive_contact_nodes=len(inactive_contacts))
    states.write_candidate(destination, fields, metadata, target_document, [source],
                           target_biases['drain'], 'reference_frame_candidate')
    return dict(source_origin_V=source_origin, target_origin_V=target_origin,
                physical_drain_voltage_V=physical_drain_voltage,
                solver_drain_voltage_V=target_biases['drain'],
                Vds_V=target_biases['drain']-target_biases['source'],
                inactive_contact_nodes=len(inactive_contacts), qualified=False)
