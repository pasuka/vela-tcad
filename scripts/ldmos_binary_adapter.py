"""Experimental VDS1 adapter for the frozen linked D5 driver; CSV default stays."""
from pathlib import Path
import math
import hashlib
import dd_state_binary as binary


def translated(data,subtract,active):
    if not math.isfinite(subtract):raise ValueError('Invalid frame')
    result=[dict(r) for r in data]
    for r in result:
        live=int(r['node_id']) in active
        for carrier in ('electron','hole'):
            ref=f'{carrier}_qf_reference_V';inc=f'{carrier}_qf_increment_V'
            if ref in r:
                a,b=float(r[ref]),float(r[inc]);new=a-subtract
                r[ref]=new if live else 0.
                r[inc]=math.fsum((a,-subtract,-new,b)) if live else 0.
        for key in ('psi','phin','phip'):
            r[key]=float(r[key])-subtract if key=='psi' or live else 0.
    return result


def install_binary(b,store=None,extension='.vds',codec=binary):
    if store is not None and extension!='.vds':raise ValueError('Memory transport requires VDS1')
    original_rows,original_digest,original_prepare,original_csv=b.rows,b.digest,b.prepare_config,b.csv_out
    def actual(path):
        p=Path(path)
        if p.suffix=='.csv' and p.resolve().is_relative_to(b.HERE.resolve()):
            alternate=p.with_suffix(extension)
            if alternate.exists() or store is not None and store.contains(alternate):return alternate
        return p
    def rows(path):
        path=actual(path)
        if store is not None and store.contains(path):return binary.decode(store.get(path))
        return codec.read(path) if path.suffix==extension else original_rows(path)
    def digest(path):
        path=actual(path)
        if store is not None and store.contains(path):return hashlib.sha256(store.get(path)).hexdigest()
        return original_digest(path)
    def prepare(*args,**kwargs):
        cfg=original_prepare(*args,**kwargs)
        cfg['sweep']['initial_state_file']=str(actual(cfg['sweep']['initial_state_file']))
        if cfg['sweep'].get('write_state_file'):
            cfg['sweep']['write_state_file']=str(Path(cfg['sweep']['write_state_file']).with_suffix(extension))
        return cfg
    def csv_out(path,data):
        if Path(path).parent.name=='predictor_inputs':
            target=Path(path).with_suffix(extension)
            if store is not None:store.add(target,binary.encode(data))
            else:codec.write(target,data)
        else:original_csv(path,data)
    def translate(source,destination,subtract):
        mesh=b.read(Path(b.BASE['mesh_file'].replace('@workspace',str(b.ROOT))))
        si={int(r['id']) for r in mesh['regions'] if r['material'].lower() in ('si','silicon')}
        active={int(n) for t in mesh['triangles'] if int(t['region_id']) in si for n in t['node_ids']}
        codec.write(Path(destination).with_suffix(extension),translated(b.rows(source),subtract,active))
    b.rows=rows;b.digest=digest;b.prepare_config=prepare;b.csv_out=csv_out;b.translate=translate
    original_point=b.CheckpointSweep.point
    def point(self,dest,bias):
        result=original_point(self,dest,bias);result['state']=str(actual(result['state']));return result
    b.CheckpointSweep.point=point
    original_child=b.FrameSweep.child
    def child(self,parent,*args,**kwargs):
        parent=actual(parent)
        result,dest=original_child(self,parent,*args,**kwargs)
        if store is not None and not b.good(result,dest,args[1]):store.persist(parent)
        return result,dest
    b.FrameSweep.child=child
    original_predict=b.Sweep.child
    def predict(self,*args,**kwargs):
        result=original_predict(self,*args,**kwargs)
        metadata=self.ledger['runs'][-1].get('outer_predictor')
        if metadata:metadata['predicted_state']=str(actual(metadata['predicted_state']))
        return result
    b.Sweep.child=predict
    original_frame=b.FrameSweep.change_frame
    def frame(self):
        original_frame(self);self.ledger['accepted_state']=str(actual(self.ledger['accepted_state']))
        for event in self.ledger['frame_events']:
            for key in ('translated_seed','canonical_for_comparison'):
                if key in event:event[key]=str(actual(event[key]))
        self.save()
    b.FrameSweep.change_frame=frame
    return actual
