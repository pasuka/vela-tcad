#!/usr/bin/env python3
"""Serial SLOT resistor campaign from local equilibrium, with independent proofs.

Outer coordinates are adaptive: a coupled proposal is independently reclosed as
a device, then its actual load-line coordinate is confirmed. This is not proof
of solving each originally proposed external voltage. Native geometry resources
remain explicit inputs; no native solution is used for initialization.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import run_slot_ldmos_loadline_continuation as ll
import slot_reference_frame as frame
import slot_loadline_step_control as steps
from dataclasses import asdict


class Stagnation(ll.NumericalFailure):
    """Qualified device state without resolvable forward continuation progress."""


def progress_policy(args):
    if not getattr(args, 'progress_guard', False):
        return None
    policy = dict(outer_floor_V=10*ll.TOLERANCE, inner_floor_V=1e-9,
                  current_relative=1e-5)
    return policy


def assess_progress(previous, current, previous_outer, outer, policy):
    """Explicit engineering guard, not a proven floating-point error bound."""
    di = current['current_A_per_um'] - previous['current_A_per_um']
    dv = current['inner_V'] - previous['inner_V']
    du = outer - previous_outer
    if not all(math.isfinite(v) for v in (di, dv, du)):
        raise ll.EvidenceError('nonfinite continuation progress')
    current_floor = policy['current_relative']*abs(previous['current_A_per_um'])
    outer_floor = max(policy['outer_floor_V'], ll.RESISTANCE*current_floor)
    # A nearly vertical, increasing-current branch remains eligible.
    resolved_state = dv >= policy['inner_floor_V'] or di > current_floor
    return dict(resolved=du >= outer_floor and resolved_state,
                delta_inner_V=dv, delta_outer_V=du, delta_current_A_per_um=di,
                outer_floor_V=outer_floor, current_floor_A_per_um=current_floor)


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def origin_for(voltage):
    # Previously used origins, scheduled before the historical cancellation
    # barriers. Their qualification must be repeated for this configuration.
    return (36. if voltage >= 35 else 35. if voltage >= 33.5 else
            25. if voltage >= 30 else 15. if voltage >= 28 else
            5. if voltage >= 3 else 0.)


def proposed_span(inner, current):
    if inner < .01:
        return .025
    if inner < 3:
        return .1
    if inner < 35:
        return .25
    return max(.025, .4 * ll.RESISTANCE * abs(current))


def check_step(previous, current, outer, previous_outer):
    if not all(math.isfinite(v) for v in (previous, current, outer, previous_outer)):
        raise ll.EvidenceError("nonfinite branch coordinate")
    if outer <= previous_outer:
        raise ll.NumericalFailure("external voltage did not increase")
    # Near equilibrium a relative current ratio is meaningless.
    if previous > 1e-13 and not 1 < current / previous <= 2:
        raise ll.NumericalFailure("current growth outside (1, 2]")


def backend_for(base, args):
    backend = ll.SolverBackend(base, args.repo_root, args.runner,
        timeout_seconds=args.timeout_seconds, max_newton_iterations=40,
        carrier_corrector_passes=8, carrier_corrector_mode="density_then_qf",
        restore_poly_equilibrium=True)
    backend.set_deadline(args.deadline)
    backend.corrector_stagnation_passes = getattr(args, 'corrector_stagnation_passes', 0)
    return backend


def stop_reason(current, physical_outer, criterion, outer_goal):
    if not all(math.isfinite(v) for v in (current, physical_outer, criterion, outer_goal)):
        raise ll.EvidenceError("nonfinite stopping coordinate")
    if abs(current) >= criterion:
        return "current_criterion"
    if physical_outer >= outer_goal - ll.TOLERANCE:
        return "outer_voltage_goal"
    return None


def predictor_parent(records, base, point, span, ratio, *, outer=None):
    if ratio is None:
        return None
    if outer is not None:
        eligible = [r for r in records if r['base'] == base and r['outer_V'] < outer-1e-6
                    and r['point']['current_A_per_um'] < point['current_A_per_um']]
        if not eligible:
            return None
        stable = [r for r in eligible if ratio*(outer-r['outer_V']) >= span]
        return (max(stable,key=lambda r:r['outer_V']) if stable else
                min(eligible,key=lambda r:r['outer_V']))
    eligible = [r for r in records if r['base'] == base
                and r['point']['inner_V'] < point['inner_V'] - 1e-10]
    if not eligible:
        return None
    stable = [r for r in eligible if ratio * (point['inner_V']-r['point']['inner_V']) >= span]
    return (max(stable, key=lambda r: r['point']['inner_V']) if stable
            else min(eligible, key=lambda r: r['point']['inner_V']))


def seed_scalar(control, backend, point, previous, outer):
    """Import only revalidated parents, then obtain an explicit circuit proof."""
    parents = ([previous] if previous is not None else []) + [point]
    for original in parents:
        backend.validate(original)
        if original['active_ids'] != point['active_ids']:
            raise ll.EvidenceError("predictor parent active rows differ")
    for original in parents:
        p = ll.clone(original)
        p['id'] = len(control.snapshot['inner_ledger'])
        control.snapshot['inner_ledger'].append(p)
    control.snapshot['active_ids'] = point['active_ids']
    control.save()
    control.accept_outer(control.snapshot['inner_ledger'][-1], outer)


def actual_outer_candidate(backend, point, outer, target, folder, *, previous,
                           predictor_ratio, max_inner_step, predicted_inner,
                           failed_trials):
    """Propose a device state; only the caller's independent circuit proof accepts it.

    No exact proposed-supply closure is asserted. One predictor attempt and one
    qualified-parent fallback are allowed; the outer controller owns shrinking.
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    bias = min(point['inner_V'] + max_inner_step,
               predicted_inner if predicted_inner is not None else
               point['inner_V'] + target - outer)
    if not math.isfinite(bias) or bias <= point['inner_V'] + 1e-9:
        raise ll.NumericalFailure('actual_outer: no resolvable forward device step')
    modes = ['predicted', 'parent'] if previous is not None and predictor_ratio is not None else ['parent']
    attempts = []
    for mode in modes:
        parent_hash = previous['point']['state_sha256'] if mode == 'predicted' else None
        key = dict(seed_sha256=point['state_sha256'], previous_sha256=parent_hash, mode=mode)
        duplicate = any(all(f.get(k) == v for k, v in key.items()) and
                        abs(f['inner_V']-bias) <= 1e-9 for f in failed_trials)
        if duplicate:
            attempts.append(dict(**key, inner_V=bias, failure='previously_failed_same_seed_voltage'))
            continue
        try:
            if mode == 'predicted':
                p = backend.solve_predicted(bias, point, previous['point'], folder / mode, predictor_ratio)
            else:
                p = backend.solve(bias, point, folder / mode)
            backend.validate(p)
        except ll.NumericalFailure as exc:
            entry = dict(**key, inner_V=bias, failure=str(exc))
            failed_trials.append(entry)
            attempts.append(entry)
            ll.atomic_json(folder / 'device_trials.json', attempts)
            continue
        attempts.append(dict(**key, inner_V=bias, qualified=True))
        p['id'] = point.get('id', 0) + 1
        qualified_file = folder / 'qualified_device.json'
        ll.atomic_json(qualified_file, p)
        log = folder / 'device_trials.json'
        ll.atomic_json(log, attempts)
        return p, dict(role='qualified_device_at_actual_outer', proposed_outer_V=target,
            actual_outer_V=p['inner_V'] + ll.RESISTANCE*p['current_A_per_um'],
            recovered=len(attempts) > 1, evidence=ll.hashes([log, qualified_file]))
    ll.atomic_json(folder / 'device_trials.json', attempts)
    raise ll.NumericalFailure('actual_outer: bounded device attempts exhausted')


def bounded_candidate(backend, point, outer, target, folder, *, candidate_seconds=None,
                      actual_outer=False, predicted_inner=None, failed_trials=None, **kwargs):
    global_deadline = backend.deadline
    local_deadline = (min(global_deadline, time.monotonic()+candidate_seconds)
                      if candidate_seconds is not None else global_deadline)
    backend.set_deadline(local_deadline)
    try:
        if actual_outer and point['inner_V'] + backend.base.get('potential_origin_V', 0.) < 38.4:
            return actual_outer_candidate(backend, point, outer, target, folder,
                previous=kwargs.get('previous'), predictor_ratio=kwargs.get('predictor_ratio'),
                max_inner_step=kwargs.get('max_inner_step', .1), predicted_inner=predicted_inner,
                failed_trials=failed_trials if failed_trials is not None else [])
        return candidate(backend, point, outer, target, folder, **kwargs)
    except ll.BudgetExhausted as exc:
        if time.monotonic() >= global_deadline:
            raise
        raise ll.NumericalFailure('candidate_budget_exhausted: ' + str(exc)) from exc
    finally:
        backend.set_deadline(global_deadline)


def candidate(backend, point, outer, target, folder, *, previous=None, predictor_ratio=None,
              max_inner_step=.25, coupled_inner_step=None, zero_converged_device_rhs=False, conservation_runner=None,
              state_predictor='secant', predictor_outer_coordinate=False):
    if state_predictor not in ('secant', 'tangent'):
        raise ValueError('candidate state predictor must be secant or tangent')
    high = point["inner_V"] + backend.base.get("potential_origin_V", 0.) >= 38.4
    if not high:
        # The audited scalar controller routes device failures through the same
        # bounded density/QF recovery as high-voltage independent qualification.
        # The executable's nested controller does not have that Python recovery.
        base, identity = ll.freeze_inputs(backend.base, backend.repo, backend.runner)
        remaining = backend.deadline - time.monotonic()
        control = ll.Campaign(folder / 'scalar', identity, [target], point['state_file'],
            point['inner_V'], 1, backend, initial_outer=outer,
            max_attempts=40, max_shots=8, max_inner_retries=4,
            max_inner_step=min(max_inner_step, target-outer), max_total_seconds=remaining,
            bracket_margin=1e-6, predictor_max_ratio=predictor_ratio)
        if predictor_ratio is not None:
            seed_scalar(control, backend, point, previous['point'] if previous else None, outer)
        if not control.run():
            failure = control.snapshot['failure']
            if failure.startswith(('EvidenceError', 'runtime_environment')):
                raise ll.EvidenceError(failure)
            raise ll.NumericalFailure('audited scalar loadline: ' + failure)
        entry = control.snapshot['outer_ledger'][-1]
        qualified = control.snapshot['inner_ledger'][entry['device_id']]
        backend.validate(qualified)
        backend.validate_confirmation(qualified, target, entry['confirmation'])
        return qualified, dict(config=entry['confirmation']['config_file'],
            updates=qualified['coupled_newton_updates'],
            recovered=bool(control.snapshot['events']),
            evidence=ll.hashes([control.output / 'checkpoint.json']))
    if state_predictor == 'secant' and predictor_ratio is not None and previous is not None:
        backend.validate(point)
        backend.validate(previous['point'])
        gap = outer - previous['outer_V']
        if gap <= 0:
            raise ll.EvidenceError("invalid high-voltage predictor parent spacing")
        ratio = min(predictor_ratio, (target-outer)/gap)
        predicted_inner = point['inner_V'] + ratio*(point['inner_V']-previous['point']['inner_V'])
        predicted = folder / 'predicted.h5'
        backend.helper.write_secant_restart(Path(previous['point']['state_file']), Path(point['state_file']),
            predicted, previous['point']['inner_V'], point['inner_V'], predicted_inner, predictor_ratio,
            document=backend.base, config_directory=backend.repo,
            **({'ratio_override': ratio} if predictor_outer_coordinate else {}))
        metadata = folder / 'predictor.json'
        ll.atomic_json(metadata, dict(qualified=False, role='initial_guess_only', ratio=ratio,
            target_inner_V=predicted_inner, parent_hashes=ll.hashes([
                previous['point']['state_file'], point['state_file']]), predicted_sha256=ll.sha256(predicted)))
        if predictor_outer_coordinate:
            details=load(metadata)
            details.update(ratio_coordinate='outer', previous_outer_V=previous['outer_V'],
                           parent_outer_V=outer, target_outer_V=target)
            ll.atomic_json(metadata,details)
        try:
            p, proof = candidate(backend, {**point, 'state_file': str(predicted), 'inner_V': predicted_inner}, outer, target,
                                 folder / 'predicted_trial', max_inner_step=max_inner_step, coupled_inner_step=coupled_inner_step,
                                 zero_converged_device_rhs=zero_converged_device_rhs, conservation_runner=conservation_runner)
        except ll.NumericalFailure as exc:
            failure = folder / 'prediction_failure.json'
            ll.atomic_json(failure, dict(failure=str(exc), qualified=False))
            p, proof = candidate(backend, point, outer, target, folder / 'parent_retry', max_inner_step=max_inner_step, coupled_inner_step=coupled_inner_step,
                                 zero_converged_device_rhs=zero_converged_device_rhs, conservation_runner=conservation_runner)
            proof['recovered'] = True
            proof['evidence'].update(ll.hashes([failure]))
        proof['evidence'].update(ll.hashes([metadata, predicted, previous['point']['state_file'], point['state_file']]))
        return p, proof
    doc, config, directory = backend.helper.render_attempt(backend.base, folder / "c",
        target, Path(point["state_file"]), 1, 40)
    span = target - outer
    doc["sweep"]["external_circuit"].update(enabled=True, solver="coupled_newton",
        initial_inner_voltage_V=point["inner_V"], max_iterations=24,
        coupled_positive_density_updates=True, coupled_equation_tolerance=1e-8,
        coupled_line_search_mode="merit", coupled_linear_solver="direct_bordered_qr",
        coupled_damping_factor=1., coupled_max_line_search_steps=12,
        coupled_apply_device_update_limit=True, max_inner_voltage_step_V=.02 if coupled_inner_step is None else coupled_inner_step,
        coupled_initial_outer_step_V=span, coupled_max_outer_step_V=span,
        coupled_min_outer_step_V=min(1e-6, span), coupled_max_step_retries=0)
    if state_predictor == 'tangent':
        # The executable constructs the full-state tangent at the qualified
        # parent, never at the Python secant guess. Strict reclosure is unchanged.
        doc['sweep']['external_circuit']['coupled_state_predictor'] = 'tangent'
    if zero_converged_device_rhs:
        # Candidate-only tangent policy. Keep the original residual tolerance;
        # final device reclosure and circuit confirmation use the untouched base.
        doc["sweep"]["external_circuit"]["coupled_inexact_device_forcing"] = dict(
            enabled=True, max_equation_tolerance=1e-8, load_activation_ratio=100.,
            zero_converged_residual=True)
    # The restart is paired with its own inner coordinate. Solve the new
    # supply directly; a predicted state must not be reclosed at the old supply.
    doc["sweep"]["bias_points"] = [target]
    doc["sweep"]["continuation"] = {"arclength": {"enabled": False}}
    _, config, directory, passed = backend._write_run((doc, config, directory))
    if not passed:
        raise ll.NumericalFailure("coupled resistor candidate did not converge")
    row = ll.read_rows(directory / "iv.csv")[-1]
    inner = ll.finite(row["inner_voltage_V"])
    conservation = None
    if conservation_runner is not None:
        import slot_candidate_conservation
        conservation = slot_candidate_conservation.audit(backend, directory / "final_state.h5",
            inner, folder / "conservation", conservation_runner)
    projected = folder / "poly.h5"
    backend.helper.restore_poly_candidate(backend.base, directory / "final_state.h5",
        projected, backend.repo, inner)
    seed = dict(state_file=str(projected), inner_V=inner, active_ids=point["active_ids"])
    qualified = backend.solve(inner, seed, folder / "q", same_bias=True)
    backend.validate(qualified)
    if int(row["newton_iterations"]) + qualified["coupled_newton_updates"] <= 0:
        raise ll.NumericalFailure("new point has no full Newton updates")
    evidence = ll.hashes([config, directory / "iv.csv", directory / "final_state.h5", projected])
    if conservation is not None:
        ll.verify_hashes(conservation['evidence'])
        evidence.update(conservation['evidence'])
        evidence.update(ll.hashes([folder / "conservation/summary.json"]))
    return qualified, dict(config=str(config), updates=int(row["newton_iterations"]),
                          conservation=conservation, evidence=evidence)


def allow_driver_update(previous, current):
    """Only this orchestrator may change; physical inputs and helpers stay frozen."""
    old = ll.clone(previous)
    path = str(Path(__file__).resolve())
    if path not in old['files'] or path not in current['files']:
        return False
    old['files'][path] = current['files'][path]
    return old == current


def confirm(backend, point, folder):
    outer = point["inner_V"] + ll.RESISTANCE * point["current_A_per_um"]
    proof = backend.confirm(point, outer, folder)
    backend.validate_confirmation(point, outer, proof)
    return outer, proof


def step_policy(args):
    if not getattr(args, 'port_step_control', False):
        return None
    return steps.Policy(outer_max_V=getattr(args, 'outer_max_step', .2),
        accepted_inner_max_V=getattr(args, 'accepted_inner_max_step', .2),
        trial_inner_max_V=getattr(args, 'trial_inner_max_step', .1),
        current_fraction=getattr(args, 'step_current_fraction', .25),
        current_floor_A_per_um=getattr(args, 'step_current_floor', 1e-16),
        growth=getattr(args, 'adaptive_growth', 1.4),
        target_iterations=getattr(args, 'step_target_iterations', 20.),
        target_seconds=getattr(args, 'step_target_seconds', 120.),
        predictor=getattr(args, 'port_predictor', 'legacy'),
        accepted_inner_decrease_max_V=getattr(args, 'accepted_inner_decrease_max_step', 0.))


def growth_policy(args, policy):
    if not getattr(args, 'current_growth_control', False):
        return None
    if (policy is None or not getattr(args, 'progress_guard', False) or
            getattr(args, 'candidate_conservation_runner', None) is None):
        raise ValueError('current growth control requires port control, progress guard and candidate audit')
    mode=getattr(args,'current_proposal_mode','requested')
    control = steps.GrowthPolicy(proposal_mode=mode,
        version='current_target_fold_v1' if mode=='current_target' else 'qualified_actual_growth_v1')
    control.validate(policy)
    return control


def feedback_mode(args, policy):
    mode = getattr(args, 'port_feedback_mode', 'proposed_outer')
    if mode not in ('proposed_outer', 'matched_actual_outer'):
        raise ValueError('unknown port feedback mode')
    if mode == 'matched_actual_outer' and (policy is None or policy.predictor != 'limited_rational'):
        raise ValueError('matched outer feedback requires port control and limited_rational predictor')
    return mode


def growth_metrics(point, proof):
    """Read bound candidate and strict evidence, never placeholder closure fields."""
    ll.verify_hashes(proof['evidence'])
    raw_file = Path(load(proof['config'])['output_csv'])
    if str(raw_file) not in proof['evidence']:
        raise ll.EvidenceError('growth observation lacks bound candidate current')
    raw = ll.finite(ll.read_rows(raw_file)[-1]['current_total_A_per_um'])
    strict = ll.finite(point['current_A_per_um'])
    if strict <= 0:
        raise ll.EvidenceError('current growth mode requires positive qualified current')
    audit = proof.get('conservation')
    if not audit:
        raise ll.EvidenceError('current growth observation lacks frozen conservation audit')
    from slot_candidate_conservation import validate_result
    validate_result(audit['result'], audit['result']['policy'])
    row = ll.read_rows(point['iv_file'])[-1]
    return dict(candidate_error=abs(raw/strict-1),
        candidate_closure=max(audit['result'][k]['ratio'] for k in ('electron','hole')),
        row_ratio=ll.finite(row['carrier_row_max_ratio']))


def quality_repair_policy(args, control):
    if not getattr(args, 'row_quality_repair', False):
        return None
    if control is None:
        raise ValueError('row quality repair requires current growth control')
    return dict(version='one_pass_row_quality_reuse_v2', max_passes=1,
                mode='carrier_only', final_gates='unchanged')


def repair_row_quality(backend, point, proof, folder, control, *, eligible):
    """Optional one-pass initial-guess repair; an original qualified state is fallback.

    Evidence/environment/budget errors propagate. A numerical rejection retains
    the original state, suppresses promotion, and still contributes full cost.
    """
    before = growth_metrics(point, proof)
    if (not eligible or before['row_ratio'] <= control.promotion_row_ratio or
            before['candidate_error'] > control.candidate_target or
            before['candidate_closure'] > control.candidate_target):
        return point, dict(status='not_needed_or_ineligible', attempted=False, before=before)
    backend.validate(point)
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=False)
    report = dict(status='running', attempted=True, before=before,
                  initial_state_sha256=point['state_sha256'])
    summary = folder/'summary.json'
    ll.atomic_json(summary, report)
    started = time.monotonic()
    old_passes = backend.corrector_passes
    try:
        backend.lineage = []
        # The backend's helper routes every nested process through its deadline.
        correction = backend.helper.correct_local_candidate(
            backend.base, Path(point['state_file']), point['inner_V'], backend.runner,
            backend.repo, folder/'block', backend.timeout, backend.repo,
            max_passes=1, mode='carrier_only', restore_poly_equilibrium=backend.restore_poly)
        ll.atomic_json(folder/'block_result.json', correction)
        report['block_passed'] = bool(correction['passed'])
        if not correction['passed']:
            raise ll.NumericalFailure('one-pass quality block did not qualify')
        # The helper has already performed the full strict solve and row probe.
        # Bind those exact files; do not solve/audit the same result again.
        repaired = backend.adopt_validation(correction['state_file'], point['inner_V'],
            point, same_bias=True)
        repaired['evidence'].update(ll.hashes([p for p in (folder/'block').rglob('*') if p.is_file()]))
        repaired['evidence'].update(point['evidence'])
        backend.validate(repaired)
        if repaired['inner_V'] != point['inner_V']:
            raise ll.EvidenceError('quality repair changed fixed internal bias')
        after = growth_metrics(repaired, proof)
        report['after'] = after
        if (after['row_ratio'] >= before['row_ratio'] or
                after['candidate_error'] > control.candidate_target):
            raise ll.NumericalFailure('quality repair did not improve row quality within candidate target')
        report.update(status='improved', final_state_sha256=repaired['state_sha256'])
    except ll.NumericalFailure as exc:
        report.update(status='numerical_rejection', failure=str(exc))
        repaired = point
    except Exception as exc:
        report.update(status='aborted', failure=type(exc).__name__+': '+str(exc))
        raise
    finally:
        backend.corrector_passes = old_passes
        report['elapsed_s'] = time.monotonic()-started
        ll.atomic_json(summary, report)
    if report['status'] == 'improved':
        repaired['evidence'].update(ll.hashes([summary, folder/'block_result.json']))
        backend.validate(repaired)
    return repaired, report


def run(args):
    output = args.output.resolve()
    args.deadline = time.monotonic() + args.max_seconds
    supplied, identity = ll.freeze_inputs(load(args.base_config), args.base_config.parent, args.runner)
    ll.validate_profile(supplied)
    if supplied["sweep"]["external_circuit"]["current_direction"] != 1:
        raise ValueError("this positive drain branch requires current_direction=+1")
    identity["files"].update(ll.hashes([args.base_config, Path(__file__),
        *[Path(__file__).with_name(n + ".py") for n in (
            "run_slot_ldmos_loadline_continuation", "run_slot_ldmos_noavalanche_continuation",
            "slot_reference_frame", "slot_state_archive", "state_archive", "slot_poly_equilibrium",
            "audit_slot_ldmos_local_continuity")]]))
    identity["criterion"] = args.criterion
    goal = getattr(args, 'outer_goal', 200000.)
    predictor = getattr(args, 'state_predictor_max_ratio', None)
    adaptive = getattr(args, 'adaptive_initial_step', None)
    growth = getattr(args, 'adaptive_growth', 1.4)
    low_cap = getattr(args, 'low_voltage_max_step', .1)
    policy = step_policy(args)
    growth_control = growth_policy(args, policy)
    if getattr(args,'current_proposal_mode','requested') != 'requested' and growth_control is None:
        raise ValueError('current-target proposals require current growth control')
    if policy is not None and policy.accepted_inner_decrease_max_V and (
            growth_control is None or growth_control.proposal_mode != 'current_target'):
        raise ValueError('internal decrease policy requires current-target control')
    row_repair = quality_repair_policy(args, growth_control)
    feedback_coordinate = feedback_mode(args, policy)
    proposal_cap = getattr(args, 'proposal_max_step', None)
    if proposal_cap is not None and (not math.isfinite(proposal_cap) or proposal_cap <= 0):
        raise ValueError('proposal step cap must be finite and positive')
    initial_checkpoint = getattr(args, 'initial_loadline_checkpoint', None)
    identity['continuation_policy'] = dict(outer_goal_V=goal, predictor_max_ratio=predictor,
        adaptive_initial_step_V=adaptive, growth=growth, low_voltage_max_step_V=low_cap)
    identity['continuation_policy']['coupled_candidate_protocol'] = 'matched_inner_direct_target_v1'
    if proposal_cap is not None:
        identity['continuation_policy']['proposal_max_step_V'] = proposal_cap
    conservation_runner = getattr(args, 'candidate_conservation_runner', None)
    if conservation_runner is not None:
        conservation_runner = Path(conservation_runner).resolve()
        identity['files'].update(ll.hashes([conservation_runner,
            Path(__file__).with_name('slot_candidate_conservation.py')]))
        identity['continuation_policy']['candidate_conservation'] = 'frozen_candidate_conservation_v1'
    progress = progress_policy(args)
    if progress is not None:
        identity['continuation_policy']['progress_guard'] = progress
    if getattr(args, 'candidate_zero_converged_device_rhs', False):
        identity['continuation_policy']['candidate_zero_converged_device_rhs'] = True
    if getattr(args, 'candidate_state_predictor', 'secant') != 'secant':
        identity['continuation_policy']['candidate_state_predictor'] = 'bordered_tangent_v1'
    if policy is not None:
        identity['continuation_policy']['port_step_control'] = asdict(policy)
        if not policy.accepted_inner_decrease_max_V:
            identity['continuation_policy']['port_step_control'].pop('accepted_inner_decrease_max_V')
        identity['files'].update(ll.hashes([Path(steps.__file__)]))
    if feedback_coordinate != 'proposed_outer':
        identity['continuation_policy']['port_feedback'] = dict(
            mode=feedback_coordinate, version='matched_actual_outer_v1',
            fallback='original_error_no_growth_no_promotion')
    if growth_control is not None:
        identity['continuation_policy']['current_growth_control'] = ll.clone(asdict(growth_control))
        if growth_control.proposal_mode == 'requested':
            identity['continuation_policy']['current_growth_control'].pop('proposal_mode')
        else:
            identity['continuation_policy']['predictor_history_coordinate'] = 'outer_current_v1'
    if row_repair is not None:
        identity['continuation_policy']['row_quality_repair'] = row_repair
    actual_outer = getattr(args, 'actual_outer_continuation', False)
    candidate_seconds = getattr(args, 'candidate_max_seconds', None)
    stagnation_passes = getattr(args, 'corrector_stagnation_passes', 0)
    reference_floor = getattr(args, 'reference_origin_floor', None)
    if reference_floor is not None:
        if not math.isfinite(reference_floor) or initial_checkpoint is None:
            raise ValueError('reference origin floor requires a finite value and audited seed window')
        identity['continuation_policy']['reference_origin_floor_V'] = reference_floor
    if actual_outer or candidate_seconds is not None or stagnation_passes:
        if actual_outer and policy is None:
            raise ValueError('actual outer continuation requires port step control')
        identity['continuation_policy']['bounded_recovery'] = dict(
            actual_outer=actual_outer, candidate_max_seconds=candidate_seconds,
            corrector_stagnation_passes=stagnation_passes, failed_trial_voltage_tolerance_V=1e-9)
    if initial_checkpoint is not None:
        identity['files'].update(ll.hashes([initial_checkpoint]))
    checkpoint = output / "checkpoint.json"
    if args.resume:
        result = load(checkpoint)
        if result["identity"] != identity:
            if not getattr(args, 'accept_driver_update', False) or not allow_driver_update(result['identity'], identity):
                raise ll.EvidenceError("resume input/tool identity mismatch")
            result.setdefault('driver_updates', []).append(dict(previous_identity=result['identity'],
                new_driver_sha256=identity['files'][str(Path(__file__).resolve())]))
            result['identity'] = identity
        ll.verify_hashes(identity["files"])
        if not result["records"]:
            raise ll.EvidenceError("no qualified restart; use a fresh output for cold initialization")
        # Validate every accepted record in its original frame, not only the tail.
        for record in result["records"]:
            check = backend_for(record["base"], args)
            check.validate(record["point"])
            check.validate_confirmation(record["point"], record["outer_V"], record["confirmation"])
        for attempt in result["attempts"]:
            if "candidate" in attempt:
                ll.verify_hashes(attempt["candidate"]["evidence"])
            if attempt["status"] == "running":
                attempt["status"] = "interrupted"
        result.update(completed=False, failure="")
    else:
        output.mkdir(parents=True, exist_ok=False)
        result = dict(schema="vela.slot.bv_campaign.v1", identity=identity,
            completed=False, failure="", records=[], attempts=[])
    save = lambda: ll.atomic_json(checkpoint, result)
    if growth_control is not None:
        result.setdefault('growth_state', steps.initial_growth_state())
    result.setdefault('invocations', []).append(dict(started_unix=time.time(),
        max_seconds=args.max_seconds, timeout_seconds=args.timeout_seconds,
        max_attempts=args.max_attempts, environment=ll.solver_environment()))
    save()
    try:
        if not result['records'] and initial_checkpoint is not None:
            imported = load(initial_checkpoint)
            if not imported.get('completed') or imported.get('failure'):
                raise ll.EvidenceError('initial loadline campaign is not complete')
            if imported['identity']['config_sha256'] != identity['config_sha256']:
                raise ll.EvidenceError('initial loadline physical configuration mismatch')
            backend = backend_for(supplied, args)
            for e in imported['outer_ledger']:
                p = imported['inner_ledger'][e['device_id']]
                backend.validate(p)
                backend.validate_confirmation(p, e['outer_V'], e['confirmation'])
                result['records'].append(dict(base=supplied, point=p, outer_V=e['outer_V'],
                    confirmation=e['confirmation'], imported=True))
            if not result['records']:
                raise ll.EvidenceError('empty initial loadline ledger')
            result['initial_loadline_checkpoint'] = str(initial_checkpoint.resolve())
            save()
        if not result["records"]:
            base = frame.translate_document(supplied, 0.)
            if any(c["bias"] != 0 for c in base["contacts"]):
                raise ValueError("cold start requires all physical contact biases zero")
            backend = backend_for(base, args)
            _, _, cold, passed = backend._write_run(
                backend.helper.render_cold_equilibrium(base, output / "cold"))
            if not passed:
                raise ll.NumericalFailure("cold local Poisson initialization failed")
            point = backend.solve(0., dict(state_file=str(cold / "final_state.h5"), inner_V=0.),
                                  output / "zero/q", same_bias=True)
            backend.validate(point)
            outer, proof = confirm(backend, point, output / "zero/k")
            if abs(outer) > ll.TOLERANCE:
                raise ll.NumericalFailure("zero-bias state does not close the zero external supply")
            result["records"].append(dict(base=base, point=point, outer_V=outer, confirmation=proof))
            save()
        retries = 0
        span = result.get('next_span_V', adaptive)
        for _ in range(args.max_attempts):
            previous = result["records"][-1]
            base, point, outer = previous["base"], previous["point"], previous["outer_V"]
            origin = base["potential_origin_V"]
            physical = point["inner_V"] + origin
            reason = stop_reason(point['current_A_per_um'], outer+origin, args.criterion, goal)
            if reason:
                result["completed"] = True
                result['stop_reason'] = reason
                break
            if time.monotonic() >= args.deadline:
                raise ll.BudgetExhausted("campaign wall time budget exhausted")
            ll.verify_hashes(identity["files"])
            backend = backend_for(base, args)
            backend.validate(point)
            number = len(result["attempts"]) + 1
            folder = output / f"a{number:04d}"
            folder.mkdir(exist_ok=False)
            attempt = dict(status="running", seed_sha256=point["state_sha256"], folder=str(folder))
            result["attempts"].append(attempt)
            save()
            started_attempt = time.monotonic()
            try:
                desired_origin = origin_for(physical)
                if reference_floor is not None:
                    desired_origin = max(desired_origin, reference_floor)
                if desired_origin != origin:
                    attempt["kind"] = "frame_requalification"
                    save()
                    shifted = frame.translate_document(base, desired_origin)
                    state = folder / "frame.h5"
                    frame.translate_candidate(point["state_file"], state, base, args.repo_root,
                                              shifted, physical)
                    backend = backend_for(shifted, args)
                    p = backend.solve(physical - desired_origin,
                        dict(state_file=str(state), inner_V=physical - desired_origin, active_ids=point["active_ids"]),
                        folder / "q", same_bias=True)
                    backend.validate(p)
                    attempt["kind"] = "frame_requalification"
                    attempt['frame_current_change_A_per_um'] = p['current_A_per_um'] - point['current_A_per_um']
                    attempt['previous_origin_V'] = origin
                    attempt['new_origin_V'] = desired_origin
                    base = shifted
                else:
                    # Enabled port control owns proposals; legacy high-bias
                    # heuristics must not silently override this policy.
                    span = ((policy.outer_max_V if policy is not None else
                             proposed_span(physical, point["current_A_per_um"])) if span is None else span)
                    if proposal_cap is not None:
                        span = min(span, proposal_cap)
                    if adaptive is not None and policy is None:
                        span = min(span, low_cap if physical < 35 else proposed_span(physical, point['current_A_per_um']))
                    if policy is not None:
                        # Only qualified same-frame history enters the port fit.
                        history = ([] if growth_control is not None and growth_control.proposal_mode=='current_target' else
                                   steps.port_samples(result['records']) if policy.predictor == 'limited_rational'
                                   else steps.samples(result['records']))
                        for sample in history:
                            record = result['records'][sample['record_index']]
                            backend.validate(record['point'])
                            backend.validate_confirmation(record['point'], record['outer_V'], record['confirmation'])
                        decision = (steps.growth_propose(result['records'], policy, growth_control,
                            result['growth_state'], span, goal-origin-outer,
                            retry_span=span if retries else None, proposal_cap=proposal_cap) if growth_control is not None else
                            steps.propose(result['records'], policy, span, goal-origin-outer))
                        span = decision['span_V']
                        attempt['step_proposal'] = decision
                    target = min(outer + span, goal-origin)
                    attempt.update(kind="coupled_resistor" if physical >= 38.4 else "scalar_loadline",
                                   proposed_outer_V=target, span_V=span)
                    save()
                    parent = predictor_parent(result['records'], base, point, target-outer, predictor,
                        **({'outer': outer} if growth_control is not None and
                           growth_control.proposal_mode=='current_target' else {}))
                    attempt['state_predictor_initialization'] = dict(
                        mode=getattr(args,'candidate_state_predictor','secant'),
                        history_coordinate='outer' if growth_control is not None and
                            growth_control.proposal_mode=='current_target' else 'inner',
                        previous_state_sha256=parent['point'].get('state_sha256') if parent else None,
                        fallback_reason='no_eligible_history' if predictor is not None and parent is None else None)
                    candidate_options = dict(
                        previous=parent, predictor_ratio=predictor,
                        max_inner_step=policy.trial_inner_max_V if policy else (low_cap if adaptive is not None else .25),
                        coupled_inner_step=policy.trial_inner_max_V if policy else None)
                    if getattr(args, 'candidate_zero_converged_device_rhs', False):
                        candidate_options['zero_converged_device_rhs'] = True
                    if conservation_runner is not None:
                        candidate_options['conservation_runner'] = conservation_runner
                    if getattr(args, 'candidate_state_predictor', 'secant') != 'secant':
                        candidate_options['state_predictor'] = args.candidate_state_predictor
                    if growth_control is not None and growth_control.proposal_mode=='current_target':
                        candidate_options['predictor_outer_coordinate'] = True
                    if actual_outer or candidate_seconds is not None:
                        p, candidate_proof = bounded_candidate(backend, point, outer, target, folder,
                            actual_outer=actual_outer, candidate_seconds=candidate_seconds,
                            predicted_inner=decision['predicted_inner_V']-origin if policy else None,
                            failed_trials=result.setdefault('failed_device_trials', []), **candidate_options)
                    else:
                        p, candidate_proof = candidate(backend, point, outer, target, folder, **candidate_options)
                    attempt["candidate"] = candidate_proof
                    if row_repair is not None:
                        provisional_outer = p['inner_V'] + ll.RESISTANCE*p['current_A_per_um']
                        resolved = assess_progress(point,p,outer,provisional_outer,progress)['resolved']
                        provisional_feedback = steps.feedback(decision, p['inner_V']+origin,
                            p['current_A_per_um'], provisional_outer-outer, policy,
                            steps.measure_cost(folder, time.monotonic()-started_attempt),
                            candidate_proof.get('recovered',False) or retries>0,
                            mode=feedback_coordinate, actual_reference_origin=base.get('potential_origin_V',0.)) if resolved else {}
                        eligible = (not candidate_proof.get('recovered',False) and retries==0 and
                            provisional_feedback.get('promotion_eligible',True) and
                            provisional_feedback.get('normalized_port_prediction_error',math.inf) <= 1. and
                            resolved and
                            steps.check_actual(decision,p['inner_V']+origin,provisional_outer+origin,
                                               p['current_A_per_um'],policy) is None and
                            provisional_outer+origin <= goal+ll.TOLERANCE)
                        p, repair_report = repair_row_quality(backend,p,candidate_proof,
                            folder/'quality_repair',growth_control,eligible=eligible)
                        attempt['row_quality_repair'] = repair_report
                        if repair_report['status']=='numerical_rejection':
                            candidate_proof['recovered'] = True
                    actual = p["inner_V"] + ll.RESISTANCE * p["current_A_per_um"]
                    if progress is not None:
                        observation = assess_progress(point, p, outer, actual, progress)
                        attempt['progress'] = observation
                        if not observation['resolved']:
                            # Retain the strict proof, but do not pollute the
                            # accepted ledger, predictor history or feedback.
                            attempt['stalled_point'] = p
                            raise Stagnation('qualified candidate has no resolvable forward progress')
                    if policy is None:
                        check_step(point["current_A_per_um"], p["current_A_per_um"], actual, outer)
                    elif growth_control is not None and point['current_A_per_um'] > policy.current_floor_A_per_um and p['current_A_per_um'] <= point['current_A_per_um']:
                        raise ll.NumericalFailure('port step control: positive branch current did not increase')
                    if actual+origin > goal+ll.TOLERANCE:
                        raise ll.NumericalFailure('actual external voltage exceeds requested goal')
                    if policy is not None:
                        rejected = steps.check_actual(decision, p['inner_V']+origin, actual+origin,
                                                      p['current_A_per_um'], policy)
                        if rejected:
                            raise ll.NumericalFailure('port step control: '+rejected)
                    attempt["proposal_residual_V"] = actual - target
                if growth_control is not None and attempt['kind'] != 'frame_requalification':
                    attempt['growth_metrics'] = growth_metrics(p, candidate_proof)
                new_outer, proof = confirm(backend, p, folder / "k")
                result["records"].append(dict(base=base, point=p, outer_V=new_outer, confirmation=proof))
                attempt["status"] = "accepted"
                if policy is not None and attempt['kind'] != 'frame_requalification':
                    cost = steps.measure_cost(folder, time.monotonic()-started_attempt)
                    attempt['attempt_cost'] = dict(cost)
                    # Include preceding rejected outer attempts, not only the
                    # final successful solve. Nested retries are in the files.
                    for rejected in reversed(result['attempts'][:-1]):
                        if rejected['status'] == 'accepted':
                            break
                        for key, value in rejected.get('attempt_cost', {}).items():
                            cost[key] += value
                    attempt['step_feedback'] = steps.feedback(decision, p['inner_V']+origin,
                        p['current_A_per_um'], new_outer-outer, policy, cost,
                        candidate_proof.get('recovered', False) or retries > 0,
                        minimum_actual_span=attempt.get('progress', {}).get('outer_floor_V', 0.),
                        mode=feedback_coordinate, actual_reference_origin=base.get('potential_origin_V',0.))
                    span = attempt['step_feedback']['next_span_V']
                    if growth_control is not None:
                        if growth_control.proposal_mode=='current_target':
                            result['growth_state']['proposal_feedback_factor']=attempt['step_feedback']['proposal_feedback_factor']
                        attempt['growth_observation'] = steps.growth_observe(
                            result['growth_state'], growth_control, decision, **attempt['growth_metrics'],
                            actual_current=p['current_A_per_um'], recovered=candidate_proof.get('recovered',False) or retries>0,
                            prediction_error=attempt['step_feedback']['normalized_port_prediction_error'],
                            prediction_valid=attempt['step_feedback'].get('promotion_eligible',True))
                elif adaptive is not None and attempt['kind'] != 'frame_requalification':
                    accepted_span = new_outer-outer
                    span = accepted_span * (1. if candidate_proof.get('recovered', False) else growth)
                elif adaptive is None:
                    span = None
                retries = 0
                result['next_span_V'] = span
                print(json.dumps(dict(accepted=True, inner_V=p["inner_V"]+base["potential_origin_V"],
                    outer_V=new_outer+base["potential_origin_V"], current_A_per_um=p["current_A_per_um"])), flush=True)
            except Stagnation as exc:
                if growth_control is not None:
                    steps.growth_failed(result['growth_state'])
                attempt.update(status='stagnation', failure=str(exc))
                attempt['attempt_cost'] = steps.measure_cost(folder, time.monotonic()-started_attempt)
                retries += 1
                # Shrink the proposal, never a numerically tiny achieved span.
                next_span = attempt['span_V']/2
                floor = attempt['progress']['outer_floor_V']
                if retries >= 6 or next_span < floor:
                    result['stop_reason'] = 'unresolved_progress'
                    raise
                span = next_span
                result['next_span_V'] = span
                print(json.dumps(dict(accepted=False, status='stagnation', next_span_V=span)), flush=True)
            except ll.NumericalFailure as exc:
                if growth_control is not None:
                    steps.growth_failed(result['growth_state'])
                attempt.update(status="timeout" if str(exc).startswith("timeout:") else "numerical_failure", failure=str(exc))
                if policy is not None:
                    attempt['attempt_cost'] = steps.measure_cost(folder, time.monotonic()-started_attempt)
                retries += 1
                span = (span if span is not None else proposed_span(physical, point["current_A_per_um"])) / 2
                result['next_span_V'] = span
                print(json.dumps(dict(accepted=False, failure=str(exc), next_span_V=span)), flush=True)
                if retries >= 6 or span < 1e-6 or attempt.get("kind") == "frame_requalification":
                    raise
            save()
            last = result['records'][-1]
            reason = stop_reason(last['point']['current_A_per_um'],
                last['outer_V']+last['base']['potential_origin_V'], args.criterion, goal)
            if reason:
                result.update(completed=True, stop_reason=reason)
                break
        else:
            raise ll.BudgetExhausted("campaign attempt budget exhausted")
    except (ll.NumericalFailure, ll.BudgetExhausted, ll.EvidenceError, ValueError, RuntimeError, OSError) as exc:
        result["failure"] = type(exc).__name__ + ": " + str(exc)
        if result["attempts"] and result["attempts"][-1]["status"] == "running":
            result["attempts"][-1].update(status="environment_or_evidence_failure", failure=str(exc))
    except KeyboardInterrupt:
        result["failure"] = "interrupted"
        if result["attempts"] and result["attempts"][-1]["status"] == "running":
            result["attempts"][-1]["status"] = "interrupted"
    save()
    print(json.dumps(dict(completed=result["completed"], failure=result["failure"],
                          accepted=len(result["records"]))), flush=True)
    return 0 if result["completed"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("base-config", "runner", "repo-root", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--accept-driver-update", action="store_true",
                        help="explicitly accept this driver update only; revalidate every previous state")
    parser.add_argument("--criterion", type=float, default=1e-7)
    parser.add_argument('--outer-goal', type=float, default=200000.)
    parser.add_argument('--initial-loadline-checkpoint', type=Path)
    parser.add_argument('--state-predictor-max-ratio', type=float)
    parser.add_argument('--adaptive-initial-step', type=float)
    parser.add_argument('--proposal-max-step', type=float,
                        help='optional proposal-only external step cap; actual acceptance caps remain unchanged')
    parser.add_argument('--adaptive-growth', type=float, default=1.4)
    parser.add_argument('--low-voltage-max-step', type=float, default=.1)
    parser.add_argument('--port-step-control', action='store_true', help='opt-in three-point port controller')
    parser.add_argument('--current-growth-control', action='store_true',
        help='opt-in staged current proposals, current-scaled outer cap and local candidate-error feedback')
    parser.add_argument('--current-proposal-mode', choices=('requested','current_target'), default='requested',
        help='opt-in current-derived supply target, independent of monotone internal voltage history')
    parser.add_argument('--accepted-inner-decrease-max-step', type=float, default=0.,
        help='explicit negative inner increment cap for current-target fold continuation; final gates unchanged')
    parser.add_argument('--row-quality-repair', action='store_true',
        help='opt-in one carrier block for qualified states missing growth row quality; final gates unchanged')
    parser.add_argument('--port-predictor', choices=('legacy','limited_rational'), default='legacy',
                        help='limited endpoint baseline with causally selected rational candidate')
    parser.add_argument('--candidate-state-predictor', choices=('secant','tangent'), default='secant',
                        help='opt-in full bordered tangent at the qualified high-voltage parent')
    parser.add_argument('--port-feedback-mode', choices=('proposed_outer','matched_actual_outer'),
                        default='proposed_outer',
                        help='opt-in same-supply feedback using the frozen port model; final gates unchanged')
    parser.add_argument('--actual-outer-continuation', action='store_true',
                        help='accept independently confirmed actual supply coordinates on adaptive low branch')
    parser.add_argument('--candidate-max-seconds', type=float)
    parser.add_argument('--candidate-conservation-runner', type=Path,
        help='Opt-in frozen raw candidate audit runner; strict reclosure is always required.')
    parser.add_argument('--candidate-zero-converged-device-rhs', action='store_true',
        help='opt-in high-voltage coupled candidate tangent policy at the original 1e-8 tolerance; final qualification unchanged')
    parser.add_argument('--progress-guard', action='store_true',
        help='opt-in engineering resolution guard; stalled strict states do not enter continuation history')
    parser.add_argument('--corrector-stagnation-passes', type=int, default=0)
    parser.add_argument('--reference-origin-floor', type=float,
                        help='minimum numerical reference origin after an explicitly audited seed import')
    parser.add_argument('--outer-max-step', type=float, default=.2)
    parser.add_argument('--accepted-inner-max-step', type=float, default=.2)
    parser.add_argument('--trial-inner-max-step', type=float, default=.1)
    parser.add_argument('--step-current-fraction', type=float, default=.25)
    parser.add_argument('--step-current-floor', type=float, default=1e-16)
    parser.add_argument('--step-target-iterations', type=float, default=20.)
    parser.add_argument('--step-target-seconds', type=float, default=120.)
    parser.add_argument("--max-attempts", type=int, default=500)
    parser.add_argument("--max-seconds", type=float, default=28800)
    parser.add_argument("--timeout-seconds", type=float, default=300)
    args = parser.parse_args()
    if args.proposal_max_step is not None and (not math.isfinite(args.proposal_max_step) or args.proposal_max_step <= 0):
        parser.error('proposal step cap must be finite and positive')
    if (args.candidate_max_seconds is not None and
            (not math.isfinite(args.candidate_max_seconds) or args.candidate_max_seconds <= 0)):
        parser.error('candidate budget must be finite and positive')
    if args.corrector_stagnation_passes < 0 or args.corrector_stagnation_passes == 1:
        parser.error('corrector stagnation passes must be zero or at least two')
    if any(not math.isfinite(v) or v <= 0 for v in
           (args.criterion, args.outer_goal, args.low_voltage_max_step,
            args.max_attempts, args.max_seconds, args.timeout_seconds)):
        parser.error("budgets and current criterion must be finite and positive")
    if (not math.isfinite(args.adaptive_growth) or args.adaptive_growth <= 1
            or any(v is not None and (not math.isfinite(v) or v <= 0)
                   for v in (args.state_predictor_max_ratio, args.adaptive_initial_step))
            or (args.adaptive_initial_step is not None and args.adaptive_initial_step > (args.outer_max_step if args.port_step_control else args.low_voltage_max_step))):
        parser.error('invalid predictor/adaptive stepping parameters')
    try:
        growth_policy(args, step_policy(args))
        quality_repair_policy(args, growth_policy(args, step_policy(args)))
        feedback_mode(args, step_policy(args))
    except ValueError as exc:
        parser.error(str(exc))
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
