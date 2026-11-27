"""Opt-in port-based step proposals, not convergence or physical error bounds.

Only qualified, same-frame campaign records may be supplied. Negative or
inconsistent slopes disable resistor amplification; they are never clipped
into an asserted monotone device law. Full state tangents are not used.
"""
from dataclasses import asdict, dataclass, replace
import csv
import math
import statistics
from pathlib import Path


@dataclass(frozen=True)
class Policy:
    outer_max_V: float = .2
    accepted_inner_max_V: float = .2
    trial_inner_max_V: float = .1
    current_fraction: float = .25
    current_floor_A_per_um: float = 1e-16
    growth: float = 1.4
    target_iterations: float = 20.
    target_seconds: float = 120.
    predictor: str = 'legacy'
    accepted_inner_decrease_max_V: float = 0.

    def __post_init__(self):
        if self.predictor not in ('legacy', 'limited_rational'):
            raise ValueError('unknown port predictor')
        if any(not math.isfinite(v) or v <= 0 for k,v in asdict(self).items()
               if k not in ('predictor', 'accepted_inner_decrease_max_V')):
            raise ValueError('step-control parameters must be finite and positive')
        if (not math.isfinite(self.accepted_inner_decrease_max_V) or
                not 0 <= self.accepted_inner_decrease_max_V <= self.accepted_inner_max_V):
            raise ValueError('internal decrease cap must be finite, nonnegative and within internal cap')
        if not 1 < self.growth <= 2 or self.current_fraction > 1:
            raise ValueError('growth must be in (1,2], current fraction in (0,1]')


@dataclass(frozen=True)
class GrowthPolicy:
    """Opt-in engineering proposal strategy; never relaxes device qualification."""
    fractions: tuple = (.04, .08, .20)
    promotion_window: int = 5
    candidate_target: float = .005
    promotion_row_ratio: float = .0009
    outer_current_scale: float = .25
    history_window: int = 5
    version: str = 'qualified_actual_growth_v1'
    proposal_mode: str = 'requested'

    def validate(self, policy):
        if self.proposal_mode not in ('requested', 'current_target'):
            raise ValueError('unknown growth proposal mode')
        if (policy.predictor != 'limited_rational' or
                not all(0 < x < policy.current_fraction for x in self.fractions) or
                tuple(sorted(set(self.fractions))) != self.fractions):
            raise ValueError('growth stages require hybrid predictor and headroom below final current guard')


def initial_growth_state():
    return dict(stage=0, clean_streak=0, multiplier=1., observations=[], failures=0)


def growth_propose(records, policy, control, state, requested, remaining, resistance=1e12,
                   *, retry_span=None, proposal_cap=None):
    """One policy governs proposed increments and strict actual acceptance.

    Candidate error/span is only a local empirical forecast. The actual current
    guard remains policy.current_fraction, independent of promotion or backoff.
    """
    control.validate(policy)
    origin = records[-1]['base'].get('potential_origin_V', 0.)
    if 'reference_origin_V' in state and state['reference_origin_V'] != origin:
        # A newly qualified reference frame must not inherit empirical error
        # coefficients or a promotion streak from the previous representation.
        state.clear()
        state.update(initial_growth_state())
    state['reference_origin_V'] = origin
    current = abs(records[-1]['point']['current_A_per_um'])
    outer_cap = max(policy.outer_max_V, control.outer_current_scale*resistance*current)
    fraction = control.fractions[state['stage']]*state['multiplier']
    forecast_policy = replace(policy, outer_max_V=outer_cap, current_fraction=fraction)
    observations = state['observations'][-control.history_window:]
    coefficient = max((o['candidate_error']/o['span_V'] for o in observations), default=0.)
    local_cap = control.candidate_target/coefficient if coefficient > 0 else None
    if control.proposal_mode == 'current_target':
        # A new target is based on current, not on the last achieved/proposed
        # ratio. Explicit retry/proposal caps are still authoritative.
        if any(not math.isfinite(v) or v <= 0 for v in (remaining, resistance, current)):
            raise ValueError('current target requires finite positive current, resistance and remaining supply')
        correction = state.get('proposal_feedback_factor', 1.)
        if not math.isfinite(correction) or not .5 <= correction <= 1.:
            raise ValueError('invalid bounded current-target feedback')
        nominal = fraction*resistance*current
        limits = dict(current_target_V=nominal*correction, outer_V=outer_cap, remaining_V=remaining)
        for key, value in (('candidate_error_V', local_cap), ('retry_V', retry_span),
                           ('proposal_cap_V', proposal_cap)):
            if value is not None:
                if not math.isfinite(value) or value <= 0:
                    raise ValueError('invalid current-target cap')
                limits[key] = value
        h = min(limits.values())
        r = records[-1]
        p = dict(inner_V=r['point']['inner_V']+origin, outer_V=r['outer_V']+origin,
                 current_A_per_um=r['point']['current_A_per_um'])
        # Constant V is an explicit port forecast, not a fitted monotonicity
        # claim. Full-state initialization and final checks remain independent.
        d = dict(span_V=h, limiting=min(limits,key=limits.get), limits=limits,
                 parent=p, predicted_inner_V=p['inner_V'],
                 predicted_current_A_per_um=p['current_A_per_um']+h/resistance,
                 current_allowance_A_per_um=fraction*current,
                 proposal_protocol='current_target_v1', reference_origin_V=origin,
                 resistance_ohm_um=resistance, nominal_target_span_V=nominal,
                 fit=dict(trusted=False,reason='explicit_constant_inner_forecast'))
    else:
        limited = min(requested, local_cap) if local_cap is not None else requested
        d = propose(records, forecast_policy, limited, remaining, resistance)
    d['prediction_current_allowance_A_per_um'] = d['current_allowance_A_per_um']
    d['current_allowance_A_per_um'] = policy.current_fraction*max(current,policy.current_floor_A_per_um)
    d['accepted_outer_max_V'] = outer_cap
    d['growth_policy'] = dict(stage=state['stage'], proposed_current_fraction=fraction,
        requested_span_V=requested, candidate_error_cap_V=local_cap,
        local_error_per_V=coefficient, actual_current_fraction_guard=policy.current_fraction)
    return d


def growth_observe(state, control, proposal, *, candidate_error, candidate_closure,
                   row_ratio, actual_current, recovered, prediction_error, prediction_valid=True):
    """Update only after strict device and actual load-line confirmation."""
    values=(candidate_error,candidate_closure,row_ratio,actual_current,prediction_error)
    if any(not math.isfinite(v) for v in values) or min(candidate_error,candidate_closure,row_ratio,prediction_error)<0:
        raise ValueError('invalid growth observation')
    clean = (prediction_valid and not recovered and candidate_error <= control.candidate_target and
             candidate_closure <= control.candidate_target and row_ratio <= control.promotion_row_ratio
             and prediction_error <= 1.)
    observation=dict(span_V=proposal['span_V'],candidate_error=candidate_error,
        candidate_closure=candidate_closure,row_ratio=row_ratio,clean=clean,
        actual_current_fraction=(actual_current-proposal['parent']['current_A_per_um'])/
                                max(abs(proposal['parent']['current_A_per_um']),1e-300),
        proposed_current_fraction=proposal['growth_policy']['proposed_current_fraction'])
    state['observations']=(state['observations']+[observation])[-control.history_window:]
    state['clean_streak']=state['clean_streak']+1 if clean else 0
    state['multiplier']=min(1.,state['multiplier']*1.4) if clean else state['multiplier']*.5
    if state['clean_streak']>=control.promotion_window:
        state['stage']=min(state['stage']+1,len(control.fractions)-1)
        state['clean_streak']=0
    return observation


def growth_failed(state):
    state.update(clean_streak=0,multiplier=state['multiplier']*.5,failures=state['failures']+1)


def samples(records):
    """Keep the newest three separated points in the current reference frame."""
    last = records[-1]
    selected = []
    for index in range(len(records)-1, -1, -1):
        r = records[index]
        if r['base'] != last['base']:
            break
        origin = r['base'].get('potential_origin_V', 0.)
        p = r['point']
        value = dict(inner_V=p['inner_V']+origin, outer_V=r['outer_V']+origin,
                     current_A_per_um=p['current_A_per_um'])
        if not all(math.isfinite(x) for x in value.values()):
            raise ValueError('nonfinite qualified port coordinate')
        if selected:
            dv = selected[-1]['inner_V']-value['inner_V']
            du = selected[-1]['outer_V']-value['outer_V']
            if dv < -1e-6 or du < -1e-6:
                break  # a turn is not a positive differential-resistance branch
            if dv <= 1e-6 or du <= 1e-6:
                continue  # duplicate/frame requalification/noise-sized spacing
        value['record_index'] = index
        selected.append(value)
        if len(selected) == 3:
            break
    return selected[::-1]


def estimate(points):
    if len(points) < 3:
        return dict(trusted=False, reason='need_three_separated_same_frame_points')
    # Shift coordinates before weighted least squares to avoid loss at large U.
    xs = [p['inner_V']-points[-1]['inner_V'] for p in points]
    ys = [p['current_A_per_um']-points[-1]['current_A_per_um'] for p in points]
    mx = sum(w*x for w,x in zip((1,2,3),xs))/6
    my = sum(w*y for w,y in zip((1,2,3),ys))/6
    g = sum(w*(x-mx)*(y-my) for w,x,y in zip((1,2,3),xs,ys))/sum(
        w*(x-mx)**2 for w,x in zip((1,2,3),xs))
    secants = [(ys[i+1]-ys[i])/(xs[i+1]-xs[i]) for i in (0,1)]
    low, high = min(g,*secants), max(g,*secants)
    if not all(math.isfinite(v) for v in (g,low,high)):
        raise ValueError('nonfinite slope estimate')
    trusted = low > 0 and high/low <= 2
    return dict(trusted=trusted, reason='positive_consistent' if trusted else 'nonpositive_or_inconsistent',
                slope_A_per_um_V=g, secants=secants, low=low, high=high)


def propose(records, policy, requested_span, remaining, resistance=1e12):
    if any(not math.isfinite(v) or v <= 0 for v in (requested_span, remaining, resistance)):
        raise ValueError('invalid step proposal input')
    if policy.predictor == 'limited_rational':
        return propose_hybrid(records, policy, requested_span, remaining, resistance)
    points = samples(records)
    fit = estimate(points)
    p = points[-1]
    allowance = policy.current_fraction*max(abs(p['current_A_per_um']), policy.current_floor_A_per_um)
    limits = dict(requested_V=requested_span, outer_V=policy.outer_max_V,
                  remaining_V=remaining, inner_V=policy.accepted_inner_max_V)
    if fit['trusted']:
        # Smaller positive slope is conservative for internal motion; the larger
        # one is conservative for predicted current growth. Neither is a bound
        # beyond the sampled interval, so actual increments are checked too.
        limits['inner_V'] *= 1+resistance*fit['low']
        limits['current_V'] = allowance*(1+resistance*fit['high'])/fit['high']
    h = min(limits.values())
    g = fit['slope_A_per_um_V'] if fit['trusted'] else 0.
    return dict(span_V=h, limiting=min(limits,key=limits.get), limits=limits, fit=fit,
                parent=p, predicted_inner_V=p['inner_V']+h/(1+resistance*g),
                predicted_current_A_per_um=p['current_A_per_um']+h*g/(1+resistance*g),
                current_allowance_A_per_um=allowance, samples=points)


def port_samples(records, count=6):
    """Causal holdout history; no interpolation across frames or reversals.

    Unlike I(v) fitting, v(U) can resolve an internal plateau. A negative
    observed increment is never clipped to manufacture a monotone history.
    """
    selected = []
    for index in range(len(records)-1, -1, -1):
        r = records[index]
        if r['base'] != records[-1]['base']:
            break
        origin = r['base'].get('potential_origin_V', 0.)
        p = dict(inner_V=r['point']['inner_V']+origin,
                 outer_V=r['outer_V']+origin,
                 current_A_per_um=r['point']['current_A_per_um'], record_index=index)
        if not all(math.isfinite(v) for v in p.values()):
            raise ValueError('nonfinite qualified port coordinate')
        if selected:
            du = selected[-1]['outer_V']-p['outer_V']
            if du < 0:
                break
            if du <= 1e-6:
                continue
            if (selected[-1]['inner_V'] < p['inner_V'] or
                    selected[-1]['current_A_per_um'] < p['current_A_per_um']):
                break
        selected.append(p)
        if len(selected) == count:
            break
    return selected[::-1]


def port_model(points, kind, horizon):
    """Shifted/scaled three-point fit, valid only on the checked interval."""
    if len(points) < 3:
        return dict(valid=False, reason='need_three_monotone_same_frame_points', kind=kind)
    p0,p1,p2 = points[-3:]
    h1=p1['outer_V']-p0['outer_V']; h2=p2['outer_V']-p1['outer_V']
    if min(h1,h2) <= 1e-6:
        return dict(valid=False, reason='unresolved_external_spacing', kind=kind)
    s1=(p1['inner_V']-p0['inner_V'])/h1
    s2=(p2['inner_V']-p1['inner_V'])/h2
    if min(s1,s2) < 0 or max(s1,s2) > 1+1e-8:
        return dict(valid=False, reason='observed_nonmonotone_port', kind=kind)
    if kind == 'limited_endpoint':
        a=max(0.,min(((2*h2+h1)*s2-h2*s1)/(h1+h2),3*s2,1.))
        return dict(valid=True,kind=kind,a=a,b=0.)
    # y=a*x/(1+b*x). Solve [1,-y] [a,b]^T = y/x after scaling.
    x0=-(h1+h2); x1=-h2
    y0=p0['inner_V']-p2['inner_V']; y1=p1['inner_V']-p2['inner_V']
    su=abs(x0); sv=max(abs(y0),abs(y1))
    if sv < 1e-15:
        return dict(valid=False,reason='unresolved_internal_change',kind=kind)
    z0=x0/su; z1=x1/su; q0=y0/sv; q1=y1/sv
    det=q0-q1
    if det == 0:
        return dict(valid=False,reason='singular_fit',kind=kind)
    cond=max(1+abs(q0),1+abs(q1))*max(abs(q0)+abs(q1),2)/abs(det)
    bb=(q1/z1-q0/z0)/det
    aa=q0/z0+q0*bb
    a=aa*sv/su; b=bb/su
    ends=(1+b*x0,1+b*horizon)
    if (not all(math.isfinite(v) for v in (a,b,cond,*ends)) or
            cond>1e8 or min(ends)<1e-3):
        return dict(valid=False,reason='pole_or_condition',kind=kind)
    derivatives=[a/(v*v) for v in ends]
    # Include x=0 too (it lies inside the checked interval).
    if min(*derivatives,a) < 0 or max(*derivatives,a) > 1+1e-12:
        return dict(valid=False,reason='derivative_shape',kind=kind)
    return dict(valid=True,kind=kind,a=a,b=b,condition=cond)


def port_increment(model, span):
    return model['a']*span/(1+model['b']*span)


def select_port_model(points, policy, horizon, resistance):
    baseline=port_model(points,'limited_endpoint',horizon)
    candidate=port_model(points,'rational3',horizon)
    scores=[]
    # Three completed one-step holdouts; each fit uses only its own past.
    for k in range(max(3,len(points)-3),len(points)):
        past=points[:k]; actual=points[k]; parent=past[-1]
        h=actual['outer_V']-parent['outer_V']
        models=[port_model(past,name,h) for name in ('limited_endpoint','rational3')]
        if not all(m['valid'] for m in models):
            scores.append(dict(valid=False, record_index=actual['record_index']))
            continue
        allowance=policy.current_fraction*max(abs(parent['current_A_per_um']),policy.current_floor_A_per_um)
        errors=[]
        for m in models:
            dv=port_increment(m,h)
            errors.append(max(abs(parent['inner_V']+dv-actual['inner_V'])/(.1*policy.accepted_inner_max_V),
                abs(parent['current_A_per_um']+(h-dv)/resistance-actual['current_A_per_um'])/(.1*allowance)))
        scores.append(dict(valid=True,record_index=actual['record_index'],limited_error=errors[0],rational_error=errors[1]))
    reason='insufficient_valid_holdouts'; chosen=baseline
    if len(scores)==3 and all(s['valid'] for s in scores):
        eb=statistics.median(s['limited_error'] for s in scores)
        er=statistics.median(s['rational_error'] for s in scores)
        # A meaningful median improvement plus no degradation on the latest
        # completed step. Near floating-point ties do not promote a model.
        if candidate['valid'] and baseline['valid'] and eb>1e-10 and er<.8*eb and scores[-1]['rational_error']<=scores[-1]['limited_error']:
            chosen=candidate;reason='rational_recent_holdout_improvement'
        else:
            reason='limited_preferred_or_candidate_ineligible'
    if not chosen['valid']:
        chosen=dict(valid=False,kind='conservative_fallback',a=1.,b=0.,reason=baseline['reason'])
        reason=baseline['reason']
    return chosen,dict(reason=reason,baseline=baseline,candidate=candidate,holdouts=scores)


def propose_hybrid(records, policy, requested_span, remaining, resistance):
    points=port_samples(records); p=points[-1]
    hard=min(requested_span,remaining,policy.outer_max_V)
    model,selection=select_port_model(points,policy,hard,resistance)
    allowance=policy.current_fraction*max(abs(p['current_A_per_um']),policy.current_floor_A_per_um)
    limits=dict(requested_V=requested_span,outer_V=policy.outer_max_V,remaining_V=remaining)
    for name,cap,fn in (
        ('inner_V',policy.accepted_inner_max_V,lambda h:port_increment(model,h)),
        ('current_V',resistance*allowance,lambda h:h-port_increment(model,h))):
        lo,hi=0.,hard
        if fn(hi)>cap:
            for _ in range(64):
                mid=(lo+hi)/2
                if fn(mid)<=cap:lo=mid
                else:hi=mid
            limits[name]=lo
        else:limits[name]=hard
    if not model['valid']:
        # No measured monotone law: both worst-case port changes are capped.
        limits['current_V']=min(limits['current_V'],resistance*allowance)
    h=min(limits.values());dv=port_increment(model,h)
    return dict(span_V=h,limiting=min(limits,key=limits.get),limits=limits,
        fit=dict(trusted=model['valid'],reason=selection['reason'],model=model,selection=selection),
        parent=p,predicted_inner_V=p['inner_V']+dv,
        predicted_current_A_per_um=p['current_A_per_um']+(h-dv)/resistance,
        current_allowance_A_per_um=allowance,samples=points,
        reference_origin_V=records[-1]['base'].get('potential_origin_V',0.),
        resistance_ohm_um=resistance)


def check_actual(proposal, inner, outer, current, policy, tolerance=1e-6):
    if not all(math.isfinite(v) for v in (inner,outer,current)):
        raise ValueError('nonfinite corrected coordinate')
    p = proposal['parent']
    if (inner < p['inner_V']-policy.accepted_inner_decrease_max_V-tolerance or
            inner-p['inner_V'] > policy.accepted_inner_max_V+tolerance):
        return 'accepted internal increment outside permitted positive branch'
    if outer <= p['outer_V']:
        return 'accepted external coordinate did not increase'
    if outer-p['outer_V'] > proposal.get('accepted_outer_max_V',policy.outer_max_V)+tolerance:
        return 'accepted external increment exceeds cap'
    if abs(current-p['current_A_per_um']) > proposal['current_allowance_A_per_um']*(1+1e-10):
        return 'accepted current change exceeds growth allowance'
    # Preserve the former campaign positive-branch requirement without its
    # separate 2x growth cap. Keep the existing increment-error precedence.
    if p['current_A_per_um'] > 1e-13 and current <= p['current_A_per_um']:
        return 'accepted current did not increase on the positive branch'
    return None


def matched_outer_prediction(proposal, actual_span, actual_reference_origin):
    """Evaluate the pre-step model at observed supply, without selecting a new fit.

    Domain checks may use the observed span but fitting uses only the saved past.
    Missing history/frame metadata must not produce a spurious clean observation.
    """
    try:
        origin=proposal['reference_origin_V']
        if (actual_reference_origin is None or not math.isfinite(origin) or
                not math.isfinite(actual_reference_origin) or origin != actual_reference_origin):
            return dict(valid=False,reason='reference_frame_mismatch')
        resistance=proposal['resistance_ohm_um']
        if not math.isfinite(resistance) or resistance <= 0:
            return dict(valid=False,reason='invalid_resistance')
        model=proposal['fit']['model']; points=proposal['samples']; parent=proposal['parent']
        if (not proposal['fit']['trusted'] or not model['valid'] or
                model['kind'] not in ('limited_endpoint','rational3') or
                len(points)<3 or points[-1] != parent):
            return dict(valid=False,reason='untrusted_model_or_history')
        domain=port_model(points,model['kind'],actual_span)
        if not domain['valid']:
            return dict(valid=False,reason=domain['reason'])
        if any(domain[k] != model[k] for k in ('a','b')):
            return dict(valid=False,reason='frozen_model_history_mismatch')
        dv=port_increment(model,actual_span)
        inner=parent['inner_V']+dv
        current=parent['current_A_per_um']+(actual_span-dv)/resistance
        if not all(math.isfinite(v) for v in (inner,current)):
            return dict(valid=False,reason='nonfinite_prediction')
        return dict(valid=True,reason='same_frozen_model_at_actual_outer',
            predicted_inner_V=inner,predicted_current_A_per_um=current,
            checked_span_V=actual_span,model=dict(model))
    except (KeyError,TypeError,ValueError,ZeroDivisionError,OverflowError):
        return dict(valid=False,reason='missing_or_invalid_model_context')


def feedback(proposal, inner, current, actual_span, policy, cost, recovered=False, *,
             minimum_actual_span=0., mode='proposed_outer', actual_reference_origin=None):
    if mode not in ('proposed_outer','matched_actual_outer'):
        raise ValueError('unknown port feedback mode')
    if (not all(math.isfinite(v) for v in (inner,current,actual_span,cost['iteration_rows'],cost['elapsed_s']))
            or actual_span <= 0 or cost['iteration_rows'] < 0 or cost['elapsed_s'] < 0):
        raise ValueError('invalid feedback coordinate or cost')
    if not math.isfinite(minimum_actual_span) or minimum_actual_span < 0 or actual_span < minimum_actual_span:
        raise ValueError('actual continuation progress below configured resolution')
    error = max(abs(inner-proposal['predicted_inner_V'])/(.1*policy.accepted_inner_max_V),
                abs(current-proposal['predicted_current_A_per_um'])/
                (.1*proposal.get('prediction_current_allowance_A_per_um',proposal['current_allowance_A_per_um'])))
    diagnostics={}
    current_target = proposal.get('proposal_protocol') == 'current_target_v1'
    if current_target:
        if actual_reference_origin != proposal['reference_origin_V']:
            raise ValueError('current-target feedback reference frame mismatch')
        # Compare at the accepted supply, but do not treat an unavailable
        # monotone history fit as a reason to suppress current-target growth.
        p=proposal['parent']
        predicted_current=p['current_A_per_um']+actual_span/proposal['resistance_ohm_um']
        diagnostics=dict(mode='current_target_v1', promotion_eligible=True,
                         original_normalized_port_prediction_error=error,
                         outer_translation_V=actual_span-proposal['span_V'])
        error=max(abs(inner-p['inner_V'])/(.1*policy.accepted_inner_max_V),
                  abs(current-predicted_current)/(.1*proposal['prediction_current_allowance_A_per_um']))
    elif mode == 'matched_actual_outer':
        matched=matched_outer_prediction(proposal,actual_span,actual_reference_origin)
        diagnostics=dict(mode=mode,original_normalized_port_prediction_error=error,
            promotion_eligible=matched['valid'],matched_prediction=matched,
            proposed_outer_V=proposal['parent']['outer_V']+proposal['span_V'],
            actual_outer_V=proposal['parent']['outer_V']+actual_span,
            outer_translation_V=actual_span-proposal['span_V'])
        if matched['valid']:
            error=max(abs(inner-matched['predicted_inner_V'])/(.1*policy.accepted_inner_max_V),
                abs(current-matched['predicted_current_A_per_um'])/
                (.1*proposal.get('prediction_current_allowance_A_per_um',proposal['current_allowance_A_per_um'])))
    factors = dict(growth=policy.growth, port_error=.8/math.sqrt(max(error,1e-30)),
                   iterations=math.sqrt(policy.target_iterations/max(1,cost['iteration_rows'])),
                   seconds=math.sqrt(policy.target_seconds/max(1e-6,cost['elapsed_s'])))
    if recovered:
        factors['recovery'] = .5
    if diagnostics and not diagnostics['promotion_eligible']:
        factors['invalid_model_no_growth'] = 1.
    factor = max(.5, min(factors.values()))
    if current_target:
        correction=min(1.,factor)
        return dict(next_span_V=proposal['nominal_target_span_V']*correction,
                    proposal_feedback_factor=correction, factor=factor, factors=factors,
                    normalized_port_prediction_error=error, cost=cost, **diagnostics)
    return dict(next_span_V=actual_span*factor, factor=factor, factors=factors,
                normalized_port_prediction_error=error, cost=cost, **diagnostics)


def measure_cost(folder, elapsed_s):
    """Count all saved nonlinear history rows, including rejected trials.

    This is a work indicator, not a count of linear solves or proof of complete
    tracing. Wall time includes corrections even if they have no history file.
    """
    files = sorted(Path(folder).rglob('iv_newton_history.csv'))
    count = 0
    for path in files:
        with path.open(encoding='utf-8-sig',newline='') as handle:
            count += sum(1 for row in csv.DictReader(handle) if float(row['iteration']) > 0)
    return dict(iteration_rows=count, history_files=len(files), elapsed_s=elapsed_s)
