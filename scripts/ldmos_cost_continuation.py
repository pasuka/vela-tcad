"""Experimental continuation feedback, independent of absolute device bias."""
import math
import time
from collections import defaultdict
from contextlib import contextmanager


def work_units(counters):
    """Deterministic work proxy; not wall seconds or hardware calibration.

    One factorization + Jacobian + residual is one nominal Newton unit.
    All attempts, including rejected candidates and recovery, enter counters.
    """
    values=[counters.get(k,0) for k in
            ('linear.factorize_calls','dd.jacobian_calls','dd.residual_calls')]
    if any(not math.isfinite(v) or v<0 for v in values):
        raise ValueError('Invalid work counters')
    return (values[0]+values[1]+.25*values[2])/2.25


def next_proposal(proposed,actual,work,min_alpha,recovery,*,minimum=.0025,maximum=.8,budget=6.):
    if not all(math.isfinite(x) and x>0 for x in (proposed,actual,minimum,maximum,budget)):
        raise ValueError('Invalid continuation limits')
    if not math.isfinite(work) or work<0 or not math.isfinite(min_alpha) or not 0<=min_alpha<=1:
        raise ValueError('Invalid progress')
    if minimum>maximum or actual>proposed+1e-10:raise ValueError('Inconsistent step')
    factor=max(.5,min(1.35,math.sqrt(budget/max(work,1.))))
    reason='cost_feedback'
    if recovery or min_alpha<.01:
        factor=.5;reason='recovery' if recovery else 'severe_damping'
    elif min_alpha<.25:
        factor=min(1.,factor);reason='limited_progress'
    # An easy tiny output fragment does not prove a larger step is affordable.
    # Preserve the pre-truncation proposal, but never grow it from such a fragment.
    if actual<.5*proposed and factor>=1.:
        value=proposed;reason='output_truncation_hold'
    else:value=actual*factor
    return min(maximum,max(minimum,value)),reason


class Timings:
    """Nested inclusive and exclusive durations; exclusive sums do not overlap."""
    def __init__(self,clock=time.perf_counter):
        self.clock=clock;self.stack=[];self.data=defaultdict(lambda:dict(calls=0,inclusive=0.,exclusive=0.))
    @contextmanager
    def stage(self,name):
        start=self.clock();frame=[0.];self.stack.append(frame)
        try:yield
        finally:
            elapsed=self.clock()-start;self.stack.pop()
            row=self.data[name];row['calls']+=1;row['inclusive']+=elapsed;row['exclusive']+=elapsed-frame[0]
            if self.stack:self.stack[-1][0]+=elapsed
    def wrap(self,name,fn):
        def measured(*args,**kwargs):
            with self.stage(name):return fn(*args,**kwargs)
        return measured
