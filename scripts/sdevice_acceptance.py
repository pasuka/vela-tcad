#!/usr/bin/env python3
"""WP10: fail-closed, same-mesh SDevice frontend acceptance (stdlib only).

References must be staged and registered first with reference_artifacts.py.
Missing VELA_REFERENCE_ARTIFACTS_DIR returns 77 for an explicit CTest skip.
No physics defaults, per-case sign corrections, or post-result thresholds.
"""
from __future__ import annotations
import argparse
import base64
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import shutil
import struct
import tempfile
import subprocess
import sys
import xml.etree.ElementTree as ET

import reference_artifacts as ra
from compare_reference_curves import compare_series, aligned_values, interpolate_at
from vela_sdevice.plt import read_plt_document
from vela_sdevice.plt.document import CURRENT_SIGN_FACTOR
from vela_sdevice.orchestrator.outputs import vtu_metadata
from analyze_pn2d_avalanche_on_bv_parity import CurvePoint, continuous_breakpoint, adjacent_slopes, slope_knee

REPO = Path(__file__).resolve().parents[1]
CONTRACT = REPO / 'configs/sdevice_acceptance_contract.json'


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite(value):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError('nonfinite measurement')
    return result


def csv_rows(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        rows = list(csv.DictReader(stream))
    if not rows or any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError('empty or incomplete CSV: ' + str(path))
    return rows


def qualified(report):
    """Partial, offline or skipped verification can never admit a reference."""
    checks = report.get('checks', [])
    return (report.get('status') == 'pass' and not report.get('partial_references')
            and report.get('mode') == 'online' and report.get('compared_with_source') and len(checks) == 6
            and {c.get('number') for c in checks} == set(range(1, 7))
            and all(c.get('status') == 'pass' and not c.get('skipped')
                    and all(i.get('status') in ('pass','info') for i in c.get('items', [])) for c in checks))


def augmented_source(original, archived):
    """Allow only Parameter insertion and necessary vector Plot output."""
    if ra.deck_tokens(original.encode()) == ra.deck_tokens(archived.encode()):
        return 'identical'
    declaration = re.compile(r'\bParameter\s*=\s*"models\.par"')
    if declaration.search(ra.strip_comments(original)) or len(declaration.findall(ra.strip_comments(archived))) != 1:
        raise ValueError('source differs beyond authorized Parameter insertion')
    stripped = declaration.sub('',archived)
    augmentation = 'Parameter="models.par" inserted; all remaining tokens identical'
    if 'ElectricField/Vector' not in ra.strip_comments(original) and 'ElectricField/Vector' in ra.strip_comments(stripped):
        # Restrict output augmentation to the existing Plot block.
        if ra.strip_comments(stripped).count('ElectricField/Vector') != 1 or not re.search(r'\bPlot\s*\{[^}]*\bElectricField/Vector',ra.strip_comments(stripped)):
            raise ValueError('vector output augmentation outside Plot')
        stripped = stripped.replace('ElectricField/Vector','')
        augmentation = 'Parameter inserted and ElectricField/Vector Plot output added; all other tokens identical'
    if ra.deck_tokens(original.encode()) != ra.deck_tokens(stripped.encode()):
        raise ValueError('physical or control source tokens differ')
    return augmentation


def native_rows(path, contact):
    # Explicit factor=1 preserves native Sentaurus signs, even after C6 freezes.
    _, rows = read_plt_document(Path(path).read_text(), sign_factor=1.0)
    return [{'bias_V': finite(r[contact]['outer_voltage_V']),
             'I': finite(r[contact]['total_current_A_per_um']), 'time': finite(r['time'])} for r in rows]


def raw_rows(paths, contact):
    output = []
    for path in paths:
        for row in csv_rows(path):
            if row.get('converged') not in ('1', 'true'):
                continue
            currents = [finite(v) for k, v in row.items() if k.startswith('contact_') and k.endswith('_total_A_per_um')]
            if len(currents) < 2:
                raise ValueError('wide terminal currents required')
            output.append({'bias_V': finite(row['bias_V']),
                           'I': finite(row[f'contact_{contact}_total_A_per_um']), 'currents': currents})
    if not output:
        raise ValueError('no accepted raw current rows')
    return output


def sign_evidence(reference, candidate):
    positive = sorted({r['bias_V']: r for r in reference if r['bias_V'] > 0}.values(), key=lambda r:r['bias_V'])[:3]
    if len(positive) != 3:
        raise ValueError('three positive reference points required')
    pairs = []
    for ref in positive:
        raw = interpolate_at(sorted((r['bias_V'],r['I']) for r in candidate),ref['bias_V'],'log_current')
        if raw is None:
            raise ValueError('first three reference biases outside raw candidate coverage')
        if raw == 0 or ref['I'] == 0:
            raise ValueError('zero current cannot establish a sign')
        pairs.append({'bias_V': ref['bias_V'], 'reference_A_per_um': ref['I'], 'raw_A_per_um': raw,
                      'sampling':'log-current interpolation between accepted raw rows',
                      'factor': 1.0 if raw * ref['I'] > 0 else -1.0})
    factors = {p['factor'] for p in pairs}
    if len(factors) != 1:
        raise ValueError('inconsistent sign at the first three points')
    return {'factor': pairs[0]['factor'], 'pairs': pairs}


def curve_metrics(reference, candidate, case, contract):
    lo, hi = case['bias_window_V']
    for series in (reference, candidate):
        if not series or min(r['bias_V'] for r in series) > lo + 1e-9 or max(r['bias_V'] for r in series) < hi - 1e-9:
            raise ValueError('curve does not cover the frozen comparison window')
    limits = contract['curve']
    kcl = max(abs(sum(r['currents'])) / max(max(map(abs, r['currents'])), limits['current_floor_A_per_um']) for r in candidate)
    if lo == hi:
        absolute = max(abs(r['I']) for r in reference + candidate)
        return {'status': 'pass' if absolute <= limits['current_floor_A_per_um'] and kcl <= limits['max_terminal_current_imbalance'] else 'fail',
                'max_absolute_current_A_per_um': absolute, 'max_kcl_relative': kcl}
    metrics = compare_series('iv', reference, candidate, CURRENT_SIGN_FACTOR, 1.0, 0.0, lo, hi,
                             reference_column='I', candidate_column='I', interpolation='log_current')
    # OOM alone uses magnitudes; opposite signs must fail independently.
    _, ref_values, cand_values = aligned_values(reference,candidate,'I','I',CURRENT_SIGN_FACTOR,1.0,0.0,lo,hi,'log_current')
    sign_ok = bool(ref_values) and all(r*c > 0 for r,c in zip(ref_values,cand_values))
    passed = metrics['available'] and metrics['points_compared'] >= limits['min_points'] and sign_ok \
        and metrics['orders_of_magnitude'] is not None and metrics['orders_of_magnitude'] <= limits['max_orders_of_magnitude'] \
        and (metrics['trend_match'] or not limits['require_trend_match']) and kcl <= limits['max_terminal_current_imbalance']
    metrics.update(status='pass' if passed else 'fail', max_kcl_relative=kcl, same_sign=sign_ok)
    metrics['derived_max_relative_error'] = 10**limits['max_orders_of_magnitude']-1
    if metrics['max_relative_error'] is None or metrics['max_relative_error'] > metrics['derived_max_relative_error']:
        metrics['status']='fail'
    return metrics


def breakdown_knee_metrics(reference, candidate, settings):
    biases = settings['biases_V']
    curves = []
    for label, rows in [('reference', reference), ('candidate', candidate)]:
        curve = {}
        for bias in biases:
            matches = [r for r in rows if abs(r['bias_V']-bias) <= 1e-9]
            if len(matches) != 1:
                raise ValueError('missing or duplicate exact breakdown-knee bias: '+label)
            curve[bias] = CurvePoint(bias, finite(matches[0]['I']))
        curves.append(curve)
    reference_knee = continuous_breakpoint(curves[0], biases, 'reference')
    candidate_knee = continuous_breakpoint(curves[1], biases, 'candidate')
    error = abs(candidate_knee-reference_knee)
    slope_knots = [slope_knee(adjacent_slopes(curve, biases, label))
                   for curve, label in zip(curves, ('reference', 'candidate'))]
    same_crossing = (slope_knots[0] is None) == (slope_knots[1] is None)
    return {'status':'pass' if error <= settings['maximum_V_break_absolute_error_V'] and same_crossing else 'fail',
            'definition':settings['definition'], 'threshold_source':settings['source'], 'biases_V':biases,
            'reference_fitted_knee_V':reference_knee, 'candidate_fitted_knee_V':candidate_knee,
            'absolute_error_V':error, 'maximum_absolute_error_V':settings['maximum_V_break_absolute_error_V'],
            'slope_crossing_V':{'reference':slope_knots[0], 'candidate':slope_knots[1]},
            'slope_outcome':'shared_no_crossing' if all(s is None for s in slope_knots) else 'both_cross' if same_crossing else 'mismatch'}


def array_values(element):
    if element is None:
        raise ValueError('missing VTU array')
    fmt = element.get('format', 'ascii')
    if fmt == 'ascii':
        values = [finite(v) for v in (element.text or '').split()]
    elif fmt == 'binary':
        raw = base64.b64decode(''.join((element.text or '').split()), validate=True)
        size = struct.unpack('<I', raw[:4])[0]
        if size != len(raw)-4 or size % 8 or element.get('type') != 'Float64':
            raise ValueError('malformed inline Float64 VTU')
        values = [finite(v) for v in struct.unpack('<' + 'd' * (size//8), raw[4:])]
    else:
        raise ValueError('unsupported VTU encoding')
    components = int(element.get('NumberOfComponents', '1'))
    if components < 1 or len(values) % components:
        raise ValueError('invalid VTU array shape')
    return [values[i:i+components] for i in range(0, len(values), components)]


def percentile(values, q):
    data = sorted(values)
    if not data:
        raise ValueError('empty active field support')
    at = (len(data)-1)*q; lo = int(at); hi = min(lo+1, len(data)-1)
    return data[lo] + (data[hi]-data[lo])*(at-lo)


def field_metrics(export, vtu, mesh, limits, expected_bias=0.0, expected_sha=None):
    tree = ET.parse(vtu); meta = vtu_metadata(vtu)
    if tree.getroot().get('byte_order','LittleEndian')!='LittleEndian' or tree.getroot().get('header_type','UInt32')!='UInt32' or tree.getroot().get('compressor'):
        raise ValueError('unsupported VTU byte order, header or compression')
    if not meta or meta.get('Converged') != 1 or abs(meta.get('Bias_V', math.inf)-expected_bias) > 1e-9:
        raise ValueError('VTU is not the requested accepted bias')
    if expected_sha and meta.get('CmdSHA256') != expected_sha:
        raise ValueError('VTU belongs to another command')
    nodes = csv_rows(export/'nodes.csv'); ids = [int(n['id']) for n in nodes]
    points = array_values(tree.find('.//Points/DataArray'))
    mesh_nodes = mesh['nodes']
    if len(set(ids)) != len(ids) or len(points) != len(nodes) or [n['id'] for n in mesh_nodes] != ids:
        raise ValueError('node identity is not bijective')
    for r, p, m in zip(nodes, points, mesh_nodes):
        if len(p)!=3 or abs(p[2])>limits['coordinate_tolerance_um']:raise ValueError('candidate mesh is not 2D')
        if max(abs(finite(r['x_um'])-p[0]),abs(finite(r['y_um'])-p[1]),abs(m['x']-p[0]),abs(m['y']-p[1])) > limits['coordinate_tolerance_um']:
            raise ValueError('node coordinates differ')
    manifest = read_json(export/'field_manifest.json')['fields']
    result = {}; errors = []
    units = {'ElectrostaticPotential': 'V', 'eDensity': 'cm^-3', 'hDensity': 'cm^-3', 'ElectricField': 'V/cm'}
    for field, unit in units.items():
        definitions = [f for f in manifest if f['name'] == field and f['region'] == 0]
        if field == 'ElectricField' and definitions:
            highest = max(f['components'] for f in definitions)
            definitions = [f for f in definitions if f['components']==highest]
        if len(definitions) != 1 or definitions[0]['mapping_status'] != 'complete' or definitions[0]['support_kind'] != 'node':
            raise ValueError('incomplete or ambiguous native field: '+field)
        definition = definitions[0]
        if definition['unit'].replace('V*cm^-1','V/cm') != unit or meta.get('Unit:'+field) != unit:
            raise ValueError('field units differ: '+field)
        records = csv_rows(export/'fields'/definition['csv_file'])
        component_keys={k for k in records[0] if k.startswith('component')}
        if component_keys!={f'component{i}' for i in range(definition['components'])}:raise ValueError('exported field component shape differs from manifest')
        if [int(r['node_id']) for r in records] != ids:
            raise ValueError('native field node mapping differs')
        reference_vectors = [[finite(r[f'component{i}']) for i in range(definition['components'])] for r in records]
        reference = [math.sqrt(sum(x*x for x in p)) if field == 'ElectricField' else p[0] for p in reference_vectors]
        candidates = array_values(tree.find(f'.//PointData/DataArray[@Name="{field}"]'))
        if len(candidates) != len(ids):
            raise ValueError('candidate field is incomplete')
        if field!='ElectricField' and any(len(p)!=1 for p in candidates):raise ValueError('scalar field shape differs')
        candidate = [math.sqrt(sum(x*x for x in p)) if field == 'ElectricField' else p[0] for p in candidates]
        if field.endswith('Density'):
            if any(v <= 0 for v in reference + candidate):
                raise ValueError('nonpositive carrier density')
            difference = [math.log10(c/r) for r,c in zip(reference,candidate) if r >= limits['carrier_active_reference_min_cm3']]
            rms_limit, p95_limit = limits['carrier_log_rmse_dex'], limits['carrier_log_p95_dex']
            if len(difference) < limits.get('carrier_active_min_nodes', 1):
                raise ValueError('insufficient active carrier support: '+field)
        else:
            difference = [c-r for r,c in zip(reference,candidate)]
            rms_limit, p95_limit = limits['potential_rmse_V'], limits['potential_p95_V']
        if not difference:raise ValueError('empty active field support')
        rms = math.sqrt(sum(v*v for v in difference)/len(difference))
        metric = {'nodes': len(difference), 'rms': rms, 'absolute_p95': percentile(list(map(abs,difference)),.95), 'maximum_diagnostic': max(map(abs,difference))}
        if field == 'ElectricField':
            # Same barycentric nodal area weighting as the qualified variants.
            positions = {n['id']:(n['x'],n['y']) for n in mesh_nodes}; weights = dict.fromkeys(ids,0.0)
            for tri in mesh['triangles']:
                a,b,c = [positions[i] for i in tri['node_ids']]
                area = abs((b[0]-a[0])*(c[1]-a[1])-(c[0]-a[0])*(b[1]-a[1]))/6
                for i in tri['node_ids']: weights[i] += area
            selected = [j for j,p in enumerate(points) if abs(p[0]-limits['junction_x_um']) <= limits['junction_half_width_um']]
            vector = definition['components'] > 1
            if vector and any(len(c) < len(r) or any(x != 0 for x in c[len(r):]) for r,c in zip(reference_vectors,candidates)):
                raise ValueError('electric-field component dimensions differ')
            squared_errors = [sum((c[i]-r[i])**2 for i in range(len(r))) if vector else (cv-rv)**2
                              for r,c,rv,cv in zip(reference_vectors,candidates,reference,candidate)]
            for label, support in [('global',list(range(len(ids)))),('junction',selected)]:
                denominator = sum(weights[ids[j]]*reference[j]**2 for j in support)
                if not support or denominator <= 0: raise ValueError('empty electric-field support')
                metric[label+'_relative_l2'] = math.sqrt(sum(weights[ids[j]]*squared_errors[j] for j in support)/denominator)
            metric['comparison'] = 'vector' if vector else 'magnitude_only'
            metric['pass'] = vector and metric['global_relative_l2'] <= limits['field_relative_l2'] \
                and metric['junction_relative_l2'] <= limits.get('junction_field_relative_l2', math.inf)
        else:
            metric['pass'] = rms <= rms_limit and metric['absolute_p95'] <= p95_limit
        result[field] = metric
        if not metric['pass']: errors.append(field)
    return {'status':'fail' if errors else 'pass', 'fields':result, 'failed_fields':errors, 'nodes':len(ids)}


def prepare_reference(base, case, destination):
    source = base/'source'/case['cmd']
    original = REPO/'reference_tcad'/case['archive']/'source'/case['cmd']
    relation = augmented_source(original.read_text(), source.read_text())
    destination.mkdir(parents=True, exist_ok=False)
    for path in (source, base/'source/models.par', base/'mesh'/case['grid']):
        shutil.copy2(path, destination/path.name)
    provenance = {'source_relation':relation, 'original_sha256':sha(original), 'cmd_sha256':sha(source),
                  'parameter_sha256':sha(base/'source/models.par'), 'grid_sha256':sha(base/'mesh'/case['grid'])}
    sidecar = original.with_suffix('.vela.json')
    if sidecar.is_file():
        shutil.copy2(sidecar,destination/sidecar.name)
        provenance['sidecar_sha256'] = sha(sidecar)
        provenance['sidecar'] = read_json(sidecar)
    return destination/source.name, provenance


def check_candidate_inputs(directory, case, provenance):
    cmd = directory/case['cmd']
    for path, expected in ((cmd,provenance['cmd_sha256']), (directory/'models.par',provenance['parameter_sha256']),
                           (directory/case['grid'],provenance['grid_sha256'])):
        if sha(path) != expected:
            raise ValueError('reused candidate input hash differs')
    sidecar = cmd.with_suffix('.vela.json')
    if ('sidecar_sha256' in provenance and (not sidecar.is_file() or sha(sidecar)!=provenance['sidecar_sha256'])) \
            or ('sidecar_sha256' not in provenance and sidecar.exists()):
        raise ValueError('reused candidate numerical sidecar differs')


def stage_policy_sha256(deck):
    policy = {key: deck[key] for key in ('solver', 'mesh_geometry', 'contacts')}
    return hashlib.sha256(json.dumps(policy, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def bv_spatial_qualification(contract, name, provenance, version, manifest_sha256, decks, mesh):
    """Admit only the independently qualified native inputs and numerical policy."""
    qualification = contract.get('biased_spatial_qualification', {}).get(name)
    if not qualification:
        return None, 'No independent avalanche spatial qualification.'
    if set(qualification.get('input_sha256', {})) != {'cmd_sha256', 'parameter_sha256', 'grid_sha256', 'sidecar_sha256'} or any(provenance.get(key) != value
            for key, value in qualification['input_sha256'].items()):
        return None, 'Avalanche spatial qualification input fingerprints differ.'
    if version != qualification['reference_version'] or manifest_sha256 != qualification['reference_manifest_sha256']:
        return None, 'Avalanche spatial qualification reference identity differs.'
    if len(decks) != 2 or stage_policy_sha256(decks[1]) != qualification['stage_policy_sha256']:
        return None, 'Avalanche spatial qualification solver/geometry/contact policy differs.'
    expected = qualification['mesh']
    if (mesh['nodes']['count'] != expected['nodes'] or mesh['elements']['count'] != expected['cells']
            or mesh['qualification']['non_obtuse'] != expected['non_obtuse']):
        return None, 'Avalanche spatial qualification mesh scope differs.'
    anchors = qualification['anchors']
    if not anchors or len({a['bias_V'] for a in anchors}) != len(anchors) or any(a['bias_V'] >= 0 for a in anchors):
        raise ValueError('invalid avalanche qualification anchors')
    required = qualification['required_biases_V']
    if not required or len(set(required)) != len(required) or not set(required).issubset({a['bias_V'] for a in anchors}):
        raise ValueError('invalid required avalanche qualification biases')
    required_limits = {'potential_rmse_V', 'potential_p95_V', 'carrier_active_reference_min_cm3',
        'carrier_active_min_nodes', 'carrier_log_rmse_dex', 'carrier_log_p95_dex', 'field_relative_l2',
        'junction_field_relative_l2', 'coordinate_tolerance_um', 'junction_x_um', 'junction_half_width_um'}
    if set(qualification['spatial']) != required_limits:
        raise ValueError('incomplete independent avalanche qualification limits')
    if any(finite(v) <= 0 for key, v in qualification['spatial'].items() if key != 'junction_x_um'):
        raise ValueError('invalid avalanche qualification limits')
    return qualification, None


def biased_field_metrics(name, case, base, output, importer, work, mesh, provenance, limits, qualification=None, reason=None):
    anchors = qualification['anchors'] if qualification else [{'bias_V': case['field_bias_V']}]
    primary = case['field_bias_V']
    if primary not in [a['bias_V'] for a in anchors]:
        raise ValueError('qualified anchors omit the required field bias')
    required = qualification['required_biases_V'] if qualification else [primary]
    if primary not in required:
        raise ValueError('qualification does not require the original field bias')
    # Inspect metadata once, then require exactly one accepted VTU per anchor.
    states = [(p, vtu_metadata(p)) for p in work.glob('stage*.vtu')]
    measurements = []
    for anchor in anchors:
        target = anchor['bias_V']
        matches = [p for p, m in states if m and m.get('Converged') == 1 and abs(m.get('Bias_V', math.inf)-target) <= 1e-9]
        native = base/'tdr'/anchor.get('native_tdr', f'{name}_multibias_{round(abs(target)/.05):04d}_des.tdr')
        if len(matches) != 1 or not native.is_file():
            measurements.append({'status':'fail', 'error':'bias-matched native TDR or accepted VTU missing', 'bias_V':target})
            continue
        if qualification and sha(native) != anchor['native_tdr_sha256']:
            raise ValueError('qualified native TDR fingerprint differs')
        export = output/'biased_reference_fields'/f'bias_{target:g}V'
        imported = subprocess.run([str(importer),'--tdr',str(native),'--export-dir',str(export)],capture_output=True,text=True)
        if imported.returncode: raise ValueError(imported.stderr)
        metric = field_metrics(export, matches[0], mesh, limits, target, provenance['cmd_sha256'])
        metric.update(bias_V=target, reference_tdr_sha256=sha(native), candidate_vtu_sha256=sha(matches[0]))
        measurements.append(metric)
    for measurement in measurements:
        measurement['acceptance_required'] = measurement['bias_V'] in required
    result = dict(next(m for m in measurements if m['bias_V'] == primary))
    result['anchors'] = measurements
    result['status'] = 'pass' if all(m['status'] == 'pass' for m in measurements if m['acceptance_required']) else 'fail'
    diagnostics = [m for m in measurements if not m['acceptance_required']]
    result['diagnostic_status'] = ('pass' if all(m['status'] == 'pass' for m in diagnostics) else 'fail') if diagnostics else 'not_required'
    if name == 'pn2d_bv':
        if qualification:
            result['qualification'] = {'id':qualification['id'], 'scope':qualification['scope'],
                'evidence':qualification['evidence'], 'limits':qualification['spatial'], 'required_biases_V':required}
        else:
            result.update(status='unqualified_threshold', reason=reason)
    return result


def run_case(name, staging, output, launcher, importer, contract, candidate_dir=None, timeout=7200):
    case = contract['cases'][name]; base = staging/case['archive']/contract['reference_version']
    manifest_path = base/'artifact_manifest.json'
    if not manifest_path.is_file():
        manifest_path = REPO/'reference_tcad'/case['archive']/f"artifact_manifest_{contract['reference_version']}.json"
    if not manifest_path.is_file():raise ValueError('registered manifest missing: '+str(manifest_path))
    manifest = read_json(manifest_path)
    if manifest['case'] != case['archive'] or manifest['version'] != contract['reference_version']:
        raise ValueError('wrong archive identity')
    # Pin ALL archived SDevice sources to the authorized repository augmentation.
    sources = output/'verified_source'; sources.mkdir(parents=True)
    for item in manifest['artifacts']:
        path = base/item['path']
        if item['path'].startswith('source/') and path.suffix == '.cmd':
            original = REPO/'reference_tcad'/case['archive']/'source'/path.name
            if ra.Deck(path.name,path.read_bytes()).kind == 'sdevice':
                augmented_source(original.read_text(), path.read_text())
            elif original.is_file() and not ra.decks_equivalent(original.read_bytes(),path.read_bytes()):
                raise ValueError('mesh source differs from repository: '+path.name)
            shutil.copy2(path, sources/path.name)
    shutil.copy2(base/'source/models.par',sources/'models.par')
    verification = ra.verify_case(manifest, manifest_path=manifest_path, backend=ra.LocalDirBackend(staging),
        staging_root=staging, report_dir=output/'reference_verification', repo_source_dir=sources, tdr_importer=[str(importer)])
    if not qualified(verification):
        return {'status':'unqualified_reference','verification':verification}
    cmd, provenance = prepare_reference(base, case, output/'candidate')
    if candidate_dir:
        cmd = candidate_dir/case['cmd']
        check_candidate_inputs(candidate_dir,case,provenance)
        code = read_json(candidate_dir/(cmd.stem+'_vela')/'run_manifest.json')['exit_code']
        if code is None:raise ValueError('candidate is still running')
    else:
        env = dict(os.environ, VELA_PYTHON=sys.executable)
        with (output/'candidate_console.log').open('w',encoding='utf-8') as log:
            proc = subprocess.Popen([str(launcher),'device',str(cmd),'--allow-lossy','--exit-on-failure','--threads','1'],
                env=env,stdout=log,stderr=subprocess.STDOUT,
                creationflags=subprocess.CREATE_NEW_PROCESS_GROUP if os.name=='nt' else 0)
            try:
                code = proc.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                proc.send_signal(signal.CTRL_BREAK_EVENT if os.name=='nt' else signal.SIGINT)
                proc.wait(timeout=60)
                code = 124
    work = cmd.parent/(cmd.stem+'_vela'); run = read_json(work/'run_manifest.json')
    plan = read_json(work/'stage_plan.json')
    if run['sha256'] != provenance['cmd_sha256']:raise ValueError('candidate manifest command fingerprint differs')
    goals = code == 0 and run['exit_code'] == 0 and len(run['stages']) == len(plan['stages']) \
        and all(s['goal_reached'] and s['runner_exit_code'] == 0 for s in run['stages'])
    reference = native_rows(base/'plt'/case['plt'],case['contact'])
    candidate = raw_rows([work/s['csv'] for s in plan['stages']],case['contact'])
    try:
        curves = curve_metrics(reference,candidate,case,contract)
    except ValueError as error:
        curves = {'status':'fail','error':str(error),'candidate_bias_range':[min(r['bias_V'] for r in candidate),max(r['bias_V'] for r in candidate)]}
    # Equilibrium VTU is retained as the first stage's accepted output.
    vtu = next(iter(sorted(work.glob('stage00_*.vtu'))),None)
    tdr = base/'tdr'/case['equilibrium_tdr']; export = output/'reference_fields'
    imported = subprocess.run([str(importer),'--tdr',str(tdr),'--export-dir',str(export)],capture_output=True,text=True)
    if imported.returncode or not vtu: raise ValueError('equilibrium field export missing: '+imported.stderr)
    fields = field_metrics(export,vtu,read_json(work/'mesh.json'),contract['spatial'],expected_sha=provenance['cmd_sha256'])
    biased = {'status':'not_required'}
    decks = [read_json(work/s['deck']) for s in plan['stages']]
    if case['require_biased_fields']:
        qualification, reason = (bv_spatial_qualification(contract, name, provenance, manifest['sentaurus_version'],
            sha(manifest_path), decks, read_json(work/'mesh_qualification.json')) if name == 'pn2d_bv' else (None, None))
        biased = biased_field_metrics(name, case, base, output, importer, work, read_json(work/'mesh.json'), provenance,
            qualification['spatial'] if qualification else contract['spatial'], qualification, reason)
    result = {'status':'pass' if goals and curves['status']=='pass' and fields['status']=='pass' and biased['status'] in ('pass','not_required') else 'fail',
              'goals_reached':goals,'curve':curves,'equilibrium_fields':fields,'provenance':provenance,
              'reference_version':manifest['sentaurus_version'],'reference_manifest_sha256':sha(manifest_path),
              'reference_plt_sha256':sha(base/'plt'/case['plt']), 'candidate_run':run,
              'biased_fields':biased}
    if name=='pn2d_bv':
        fields['threshold_scope']='equilibrium control only; does not qualify avalanche biased states'
        try:
            result['breakdown_knee'] = breakdown_knee_metrics(reference, candidate, case['breakdown_check'])
        except ValueError as error:
            result['breakdown_knee'] = {'status':'fail', 'error':str(error)}
        if result['breakdown_knee']['status'] != 'pass': result['status'] = 'fail'
    result['cmd_exit_code']=code
    result['candidate_csv_sha256']={s['csv']:sha(work/s['csv']) for s in plan['stages'] if (work/s['csv']).is_file()}
    result['mesh']={'nodes':len(read_json(work/'mesh.json')['nodes']),'cells':len(read_json(work/'mesh.json')['triangles'])}
    result['known_differences'] = [item for item in contract['known_difference_ledger'] if item['case'] == name]
    if name == 'pn2d_iv': result['sign_evidence'] = sign_evidence(reference,candidate)
    result['backend'] = {'runner_sha256':sha(launcher.parent/('vela_example_runner.exe' if os.name=='nt' else 'vela_example_runner')),
                         'materials_report':read_json(work/'materials_report.json'),
                         'stage_decks': decks}
    return result


def time_semantics(base):
    files = (base/'plt/bjt_m0_base_bjt_m0_des.plt',base/'plt/bjt_m0_collector_bjt_m0_des.plt')
    evidence = []
    for path,contact in zip(files,('base','collector')):
        rows = native_rows(path,contact);times = [r['time'] for r in rows]
        if abs(times[0]) > 1e-9 or abs(times[-1]-1.) > 1e-9 or any(b<a for a,b in zip(times,times[1:])):
            raise ValueError('native Time does not reset to 0..1 at each stage')
        evidence.append({'path':str(path),'sha256':sha(path),'rows':len(rows),'first_time':times[0],'last_time':times[-1]})
    return {'semantics':'reset_per_quasistationary_stage','evidence':evidence}


def run_gates(staging,output,launcher):
    results=[];env=dict(os.environ,VELA_PYTHON=sys.executable)
    for archive,commands,grid in (
        ('schottky_charon_sentaurus2018',['schottky_n_des.cmd'],'schottky_n_msh.tdr'),
        ('bvmethods_sentaurus2018',['aba_poisson_sdevice.cmd','aba_coupled_sdevice.cmd','external_resistor_sdevice.cmd','voltage_to_current_sdevice.cmd'],'n1_msh.tdr')):
        base=staging/archive/'v1'
        meshes=list(base.rglob(grid))
        if len(meshes)!=1:raise ValueError('check-only mesh missing or ambiguous: '+grid)
        manifest=read_json(REPO/'reference_tcad'/archive/'artifact_manifest.json')
        listed=[a for a in manifest['artifacts'] if a['path']==meshes[0].relative_to(base).as_posix()]
        if len(listed)!=1 or sha(meshes[0])!=listed[0]['sha256']:raise ValueError('check-only grid hash mismatch')
        for name in commands:
            directory=output/name;directory.mkdir()
            source=REPO/'reference_tcad'/archive/'source'/name
            text=source.read_text()
            original_grid=re.search(r'\bGrid\s*=\s*"([^"]+)"',text).group(1)
            shutil.copy2(meshes[0],directory/original_grid)
            parameter=REPO/'reference_tcad'/archive/'source/models.par'
            if not parameter.is_file():
                parameter=base/'source/models.par'
            if not parameter.is_file():
                from vela_sdevice.materials import default_models_par_path
                parameter=default_models_par_path()
            if 'Parameter' not in ra.strip_comments(text):text=text.replace('File {','File {\n Parameter="models.par"',1)
            shutil.copy2(parameter,directory/'models.par')
            cmd=directory/name;cmd.write_text(text,encoding='utf-8')
            completed=subprocess.run([str(launcher),'device',str(cmd),'--check','--allow-lossy'],capture_output=True,text=True,encoding='utf-8',env=env)
            results.append({'cmd':name,'expected_exit_code':0,'exit_code':completed.returncode,'diagnostics':completed.stdout+completed.stderr,'pass':completed.returncode==0})
    fixture=REPO/'tests/fixtures/sdevice/transportmodels_sentaurus2022/pp6_des.cmd'
    rejected=output/'transport_rejection';rejected.mkdir()
    cmd=rejected/fixture.name;shutil.copy2(fixture,cmd)
    completed=subprocess.run([str(launcher),'device',str(cmd),'--check'],capture_output=True,text=True,encoding='utf-8',env=env)
    diagnostic=completed.stdout+completed.stderr
    from vela_sdevice.cmd_parser import parse_cmd_file
    ids=parse_cmd_file(cmd)['verdict']['fatal_ids']
    results.append({'cmd':str(cmd),'expected_exit_code':1,'exit_code':completed.returncode,'diagnostics':diagnostic,
                    'fatal_ids':ids,'pass':completed.returncode==1 and 'E105' in diagnostic and set(ids)=={'E105','E109'}})
    return {'status':'pass' if all(r['pass'] for r in results) else 'fail','gates':results}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',choices=('pn2d_0v','pn2d_iv','pn2d_bv','gates'),required=True)
    parser.add_argument('--reference-dir',type=Path,default=os.environ.get('VELA_REFERENCE_ARTIFACTS_DIR'))
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--launcher',type=Path,required=True)
    parser.add_argument('--tdr-importer',type=Path,required=True)
    parser.add_argument('--candidate-dir',type=Path,help='Audit an already finished run, after exact cmd/par/grid/sidecar hash checks')
    parser.add_argument('--timeout-seconds',type=float,default=7200.,help='Bound each fresh candidate run; timeout never passes (full BV retains 401 Plot biases)')
    args = parser.parse_args(argv)
    if args.timeout_seconds<=0:parser.error('--timeout-seconds must be positive')
    if args.reference_dir is None:
        print('SKIP: VELA_REFERENCE_ARTIFACTS_DIR is not set; native references unavailable.')
        return 77
    args.output_dir.mkdir(parents=True,exist_ok=True)
    args.output_dir = Path(tempfile.mkdtemp(prefix=args.case+'_',dir=args.output_dir))
    try:
        result = run_gates(args.reference_dir.resolve(),args.output_dir.resolve(),args.launcher.resolve()) if args.case=='gates' else \
            run_case(args.case,args.reference_dir.resolve(),args.output_dir.resolve(),args.launcher.resolve(),args.tdr_importer.resolve(),read_json(CONTRACT),
                     args.candidate_dir.resolve() if args.candidate_dir else None,args.timeout_seconds)
    except (ValueError,OSError,KeyError,StopIteration,ET.ParseError,subprocess.TimeoutExpired,struct.error) as error:
        result = {'status':'fail','error':str(error)}
    result.update(case=args.case,contract_sha256=sha(CONTRACT),units=read_json(CONTRACT)['units'])
    (args.output_dir/'acceptance_report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'case':args.case,'status':result['status'],'report':str(args.output_dir/'acceptance_report.json')}))
    return 0 if result['status']=='pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())
