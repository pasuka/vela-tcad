"""Bounded immutable VDS1 predictor seeds; solver outputs remain checkpoints."""
from collections import Counter, OrderedDict
from pathlib import Path
import hashlib
import dd_state_binary as binary


def replay_prediction(metadata,rows,digest):
    """Rebuild the frozen guarded-secant seed using its persisted parent states."""
    if metadata['mode']!='outer_secant_guarded_v1':raise ValueError('Unsupported prediction')
    for prefix in ('previous','current'):
        if digest(Path(metadata[prefix+'_state']))!=metadata[prefix+'_sha256']:
            raise ValueError('Prediction parent changed')
    current=rows(Path(metadata['current_state']));previous=rows(Path(metadata['previous_state']))
    if len(current)!=len(previous):raise ValueError('Prediction node mismatch')
    ratio=metadata['ratio']
    for c,p in zip(current,previous):
        if c['node_id']!=p['node_id']:raise ValueError('Prediction node order changed')
        for key in ('psi','phin','phip'):
            c[key]=format(float(c[key])+ratio*(float(c[key])-float(p[key])),'.17g')
        for field,ref,inc in [('phin','electron_qf_reference_V','electron_qf_increment_V'),
                              ('phip','hole_qf_reference_V','hole_qf_increment_V')]:
            c[ref]=c[field];c[inc]='0'
    raw=binary.encode(current)
    if hashlib.sha256(raw).hexdigest()!=metadata['predicted_sha256']:
        raise ValueError('Replayed prediction changed')
    return raw


class MemorySeeds:
    def __init__(self, limit=64*1024*1024):
        if limit<=0:raise ValueError('Positive memory budget required')
        self.limit=limit;self.data=OrderedDict();self.used=0;self.stats=Counter();self.manifest={}

    def contains(self,path):return Path(path).resolve() in self.data

    def get(self,path):return self.data[Path(path).resolve()]

    def add(self,path,raw):
        path=Path(path).resolve()
        if path.exists() or str(path) in self.manifest:raise ValueError('Seed is immutable: '+str(path))
        if not isinstance(raw,bytes):raise TypeError('Immutable bytes required')
        self.manifest[str(path)]=dict(sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw),storage='memory')
        self.data[path]=raw;self.used+=len(raw);self.stats['created']+=1
        while self.used>self.limit:self.persist(next(iter(self.data)),reason='budget_spill')
        self.stats['peak_retained_bytes']=max(self.stats['peak_retained_bytes'],self.used)

    def persist(self,path,reason='failure_evidence'):
        path=Path(path).resolve()
        if path not in self.data:return
        raw=self.data[path]
        with path.open('xb') as f:f.write(raw)
        del self.data[path];self.used-=len(raw)
        self.manifest[str(path)]['storage']=reason;self.stats[reason]+=1

    def flush(self):
        for path in list(self.data):self.persist(path)

    def audit(self):
        for name,item in self.manifest.items():
            path=Path(name);raw=self.get(path) if self.contains(path) else path.read_bytes()
            if len(raw)!=item['bytes'] or hashlib.sha256(raw).hexdigest()!=item['sha256']:
                raise ValueError('Predictor seed changed: '+name)
        return dict(stats=dict(self.stats),memory_seeds=len(self.data),memory_bytes=self.used,seeds=self.manifest)
