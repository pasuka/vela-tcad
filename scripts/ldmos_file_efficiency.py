"""Bounded experimental file reuse and cost-per-volt continuation policy."""
from collections import Counter, OrderedDict
from pathlib import Path
import math


class CompletedFileCache:
    """CSV is eligible only after its worker request finishes; JSON is not cached.

    Metadata detects ordinary writes/replacements. Final evidence verification must
    disable this cache and reread bytes; metadata is not a cryptographic guarantee.
    Values are copied because the secant predictor mutates row dictionaries.
    """
    def __init__(self, rows, digest, entries=8):
        self.raw_rows=rows;self.raw_digest=digest;self.entries=entries
        self.completed=set();self.csv=OrderedDict();self.hashes=OrderedDict()
        self.stats=Counter();self.enabled=True

    @staticmethod
    def key(path):
        p=Path(path).resolve();s=p.stat()
        return p,(s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns)

    def seal(self, directory):
        self.completed.add(Path(directory).resolve())

    def disable(self):
        self.enabled=False;self.csv.clear();self.hashes.clear()

    def get(self,path,cache,read,label,eligible,limit):
        if not self.enabled or not eligible:
            self.stats[label+'_bypass']+=1;return read(path)
        p,key=self.key(path)
        item=cache.pop(p,None)
        if item is not None and item[0]==key:
            self.stats[label+'_hit']+=1;cache[p]=item;return item[1]
        self.stats[label+'_miss']+=1
        if item is not None:self.stats[label+'_invalidated']+=1
        value=read(path)
        if self.key(p)[1]!=key:raise RuntimeError('File changed during read: '+str(p))
        cache[p]=(key,value)
        while len(cache)>limit:cache.popitem(last=False)
        return value

    def rows(self,path):
        eligible=Path(path).resolve().parent in self.completed
        value=self.get(path,self.csv,self.raw_rows,'csv',eligible,self.entries)
        return [dict(row) for row in value] if self.enabled and eligible else value

    def digest(self,path):
        return self.get(path,self.hashes,self.raw_digest,'digest',True,128)


class EfficiencyPolicy:
    """Observe relative cost/volt, with no absolute-bias split or Newton budget.

    Compare increasing/decreasing steps only when step sizes differ by >=10%.
    A >20% deterioration reverses the last size change. Recovery/damping guards
    take precedence. Output fragments do not train the trend estimate.
    """
    def __init__(self,minimum=.0025,maximum=.8):
        if not 0<minimum<=maximum:raise ValueError('Invalid limits')
        self.minimum=minimum;self.maximum=maximum;self.previous=None

    def next(self,proposed,actual,cost,alpha,recovery,frame):
        if not all(math.isfinite(x) and x>0 for x in (proposed,actual,cost)):
            raise ValueError('Invalid cost or step')
        if actual>proposed+1e-10 or not math.isfinite(alpha) or not 0<=alpha<=1:
            raise ValueError('Invalid step or damping')
        rate=cost/actual;old=self.previous
        if old is not None and old['frame']!=frame:old=None;self.previous=None
        reason='efficiency_explore';value=actual*1.35
        if recovery or alpha<.01:
            value=actual*.5;reason='recovery' if recovery else 'severe_damping';self.previous=None
        elif alpha<.25:
            value=actual;reason='limited_progress';self.previous=None
        elif actual<.8*proposed:
            value=proposed;reason='output_fragment_hold';self.previous=None
        else:
            if old is not None and rate>1.2*old['rate']:
                if actual>1.1*old['step']:
                    value=old['step'];reason='larger_step_cost_regression'
                elif actual<.9*old['step']:
                    value=old['step'];reason='smaller_step_cost_regression'
            self.previous=dict(step=actual,rate=rate,frame=frame)
        return min(self.maximum,max(self.minimum,value)),dict(reason=reason,cost_per_V=rate,previous=old)
