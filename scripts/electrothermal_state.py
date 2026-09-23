"""HDF5 state references for production four-equation control records.

Inline dictionaries remain the in-memory solver API. Files passed to production
sweeps use this codec; historical inline JSON needs explicit conversion once.
"""
import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np
try:
    from . import state_archive
except ImportError:
    import state_archive

KEYS = ('state_interleaved', 'referenced_state_interleaved',
        'electron_qf_reference_V', 'hole_qf_reference_V', 'temperature_K')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pack(path, value, metadata):
    path = Path(path)
    sequence = 0
    def visit(record):
        nonlocal sequence
        out = dict(record)
        if 'state_archive' in out:
            raise ValueError('State record already packed')
        if 'state_interleaved' in out or 'referenced_state_interleaved' in out:
            referenced = 'referenced_state_interleaved' in out
            x = np.asarray(out['referenced_state_interleaved' if referenced else 'state_interleaved'], dtype='<f8')
            if x.ndim != 1 or not x.size or x.size % 4:
                raise ValueError('Invalid four-equation layout')
            x = x.reshape((-1, 4))
            physical = np.asarray(out.get('state_interleaved', x.ravel()), dtype='<f8').reshape(x.shape).copy()
            fields = dict(psi=x[:, 0], temperature_K=x[:, 3])
            if not np.array_equal(physical[:, (0, 3)], x[:, (0, 3)]):
                raise ValueError('Inconsistent thermal psi or temperature representations')
            for k, carrier, field in ((1, 'electron', 'phin'), (2, 'hole', 'phip')):
                ref = carrier+'_qf_reference_V'
                fields[ref] = np.asarray(out[ref], dtype='<f8')
                fields[carrier+'_qf_increment_V'] = x[:, k] if referenced else x[:, k]-fields[ref]
                fields[field] = physical[:, k] if 'state_interleaved' in out else x[:, k]+fields[ref]
            if 'temperature_K' in out and not np.array_equal(out['temperature_K'],fields['temperature_K']):
                raise ValueError('Inconsistent temperature field')
            meta = dict(metadata, mode='electrothermal',
                        potential_origin_V=out.get('potential_origin_V',metadata.get('potential_origin_V',0.)),
                        record_fields=[k for k in KEYS if k in out])
            target=path.with_name(path.stem+f'_state_{sequence}.h5'); sequence+=1
            if target.exists():
                raise FileExistsError(target)
            state_archive.write(target, fields, meta)
            for k in KEYS:
                out.pop(k, None)
            out['state_archive'] = dict(file=target.name, sha256=digest(target))
        if 'diagnostic_tangent_predictor' in out:
            out['diagnostic_tangent_predictor'] = visit(out['diagnostic_tangent_predictor'])
        if 'diagnostic_predictor_candidates' in out:
            out['diagnostic_predictor_candidates'] = [visit(v) for v in out['diagnostic_predictor_candidates']]
        return out
    return visit(value)


def unpack(path, value, nodes, mesh_sha256, potential_origin_V):
    path=Path(path)
    def visit(record):
        out=dict(record)
        if 'state_archive' in out:
            if any(k in out for k in KEYS):
                raise ValueError('Mixed inline and HDF5 thermal state')
            ref=out.pop('state_archive'); target=path.parent/Path(ref['file'])
            if digest(target)!=ref['sha256']:
                raise ValueError('Thermal state file digest mismatch')
            fields,meta=state_archive.read(target,nodes,mesh_sha256)
            if meta['mode']!='electrothermal' or meta['potential_origin_V']!=potential_origin_V:
                raise ValueError('Thermal state mode or reference origin mismatch')
            physical=np.column_stack([fields[k] for k in ('psi','phin','phip','temperature_K')]).ravel().tolist()
            increment=np.column_stack([fields[k] for k in ('psi','electron_qf_increment_V','hole_qf_increment_V','temperature_K')]).ravel().tolist()
            for k in meta['record_fields']:
                if k=='state_interleaved': out[k]=physical
                elif k=='referenced_state_interleaved': out[k]=increment
                elif k in KEYS: out[k]=fields[k].tolist()
                else: raise ValueError('Unknown saved thermal field')
        elif 'state_interleaved' in out or 'referenced_state_interleaved' in out:
            raise ValueError('Production thermal state requires HDF5')
        if 'diagnostic_tangent_predictor' in out:
            out['diagnostic_tangent_predictor']=visit(out['diagnostic_tangent_predictor'])
        if 'diagnostic_predictor_candidates' in out:
            out['diagnostic_predictor_candidates']=[visit(v) for v in out['diagnostic_predictor_candidates']]
        return out
    return visit(value)


def write(path, record, metadata):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    packed=pack(path,record,metadata)
    fd,temporary=tempfile.mkstemp(prefix=path.name+'.tmp.',dir=path.parent)
    try:
        with os.fdopen(fd,'w',encoding='utf-8') as stream:
            json.dump(packed,stream,allow_nan=False,indent=2)
        os.replace(temporary,path)
    finally:
        if os.path.exists(temporary): os.unlink(temporary)


def read(path, nodes, mesh_sha256, potential_origin_V):
    path=Path(path)
    return unpack(path,json.loads(path.read_text(encoding='utf-8')),nodes,mesh_sha256,potential_origin_V)


def read_bound_record(path):
    """Read current control/output using its input mesh, never its own hash alone."""
    path=Path(path)
    value=json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(value,dict) or 'state_archive' not in value:
        return value  # Scalar reports and historical audit records are not seeds.
    cfgpath=path.parent/Path(value['state_archive'].get('input_file','input.json'))
    cfg=value if 'mesh_file' in value else json.loads(cfgpath.read_text(encoding='utf-8'))
    meshpath=Path(cfg['mesh_file'])
    if not meshpath.is_absolute(): meshpath=(path.parent if 'mesh_file' in value else cfgpath.parent)/meshpath
    mesh=json.loads(meshpath.read_text(encoding='utf-8'))
    return unpack(path,value,len(mesh['nodes']),
        state_archive.mesh_identity(mesh,cfg.get('coordinate_to_metres',1.)),cfg.get('potential_origin_V',0.))
