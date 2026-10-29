"""Frozen candidate audit. Conservation acceptance never qualifies a device."""
from pathlib import Path
import math
import run_slot_ldmos_loadline_continuation as ll


def validate_result(result, policy):
    if (not isinstance(result, dict)
            or result.get('schema') != 'frozen_candidate_conservation_v1'
            or result.get('role') != 'initial_guess_only'
            or result.get('device_qualified') is not False
            or result.get('packing') != 'preserve_referenced_contact_low_parts'
            or result.get('equations') != 'before_dirichlet_replacement'
            or result.get('policy') != policy):
        raise ll.EvidenceError('candidate conservation protocol/policy mismatch')
    accepted = True
    for name in ('electron', 'hole'):
        row = result.get(name, {})
        if not isinstance(row, dict):
            raise ll.EvidenceError('invalid candidate conservation carrier record')
        values = [row.get(k) for k in ('flux', 'source', 'mismatch', 'ratio')]
        if any(type(v) not in (int, float) or not math.isfinite(v) for v in values):
            raise ll.EvidenceError('invalid candidate conservation values')
        flux, source, mismatch, ratio = values
        scale = max(abs(flux), abs(source), policy['source_floor'])
        expected = abs(flux-source)/scale if scale else 0.
        relative = max(abs(flux), abs(source)) >= policy['source_floor']
        passed = (expected <= policy['tolerance'] if relative else
                  abs(flux-source) <= policy['absolute_tolerance'])
        if (not math.isclose(mismatch, flux-source, rel_tol=1e-12, abs_tol=0.)
                or not math.isclose(ratio, expected, rel_tol=1e-12, abs_tol=0.)
                or row.get('relative_criterion') is not relative
                or row.get('satisfied') is not passed):
            raise ll.EvidenceError('candidate conservation arithmetic/decision mismatch')
        accepted = accepted and passed
    if result.get('satisfied') is not accepted:
        raise ll.EvidenceError('candidate conservation aggregate mismatch')
    return accepted


def audit(backend, state, inner, folder, runner):
    folder = Path(folder).resolve()
    folder.mkdir(parents=True, exist_ok=False)
    document, identity = ll.freeze_inputs(backend.base, backend.repo, Path(runner))
    gate = document['solver']['global_continuity_closure']
    if gate.get('mode') != 'enforce':
        raise ll.EvidenceError('candidate audit requires unchanged final closure enforcement')
    policy = dict(tolerance=gate.get('tolerance', .01),
                  source_floor=gate.get('source_floor', 1e-12),
                  absolute_tolerance=gate.get('absolute_tolerance', 1e-12))
    driven = document.pop('sweep')['contact']
    for contact in document['contacts']:
        if contact['name'] == driven:
            contact['bias'] = inner
    state = Path(state).resolve()
    result_file = folder/'result.json'
    document.update(simulation_type='newton_candidate_conservation_probe',
                    state_format='hdf5', state_file=str(state), output_csv=str(result_file))
    config = folder/'config.json'
    ll.atomic_json(config, document)
    inputs = dict(identity['files'])
    inputs.update(ll.hashes([state, config, Path(__file__)]))
    code, timed_out, elapsed = backend._execute(Path(runner), backend.repo, config, folder, backend.timeout)
    ll.verify_hashes(inputs)
    if code != 0 or timed_out:
        raise ll.EvidenceError('candidate conservation probe did not complete')
    import json
    result = json.loads(result_file.read_text(encoding='utf-8'))
    passed = validate_result(result, policy)
    summary = dict(role='initial_guess_only', device_qualified=False,
                   conservation_satisfied=passed, strict_reclosure_required=True,
                   disposition='strict_reclosure_required' if passed else 'conservation_failed_repair_required',
                   inner_V=inner, elapsed_s=elapsed, result=result,
                   evidence={**inputs, **ll.hashes([result_file])})
    ll.atomic_json(folder/'summary.json', summary)
    return summary
