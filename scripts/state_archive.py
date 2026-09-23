"""Vela state/2 HDF5 contract. No legacy file sniffing or format fallback."""
import json
import math
import os
from pathlib import Path
import re
import tempfile
import hashlib
import struct

import h5py
import numpy as np

SCHEMA = 'vela.state/2'
MAX_BYTES = 2**30
UNITS = {name:'V' for name in ('psi', 'phin', 'phip', 'electron_quantum_potential_V',
    'electron_quantum_potential_like_V', 'electron_qf_increment_V', 'hole_qf_increment_V',
    'electron_qf_reference_V', 'hole_qf_reference_V')}
UNITS.update(electrons_m3='m^-3', holes_m3='m^-3', temperature_K='K')
REQUIRED = {'psi', 'phin', 'phip', 'electrons_m3', 'holes_m3'}
SPLIT = {'electron_qf_increment_V', 'hole_qf_increment_V',
         'electron_qf_reference_V', 'hole_qf_reference_V'}


def mesh_identity(mesh, length_m_per_internal):
    """Canonical encoding shared with stateMeshIdentity, independent of JSON spacing."""
    if not math.isfinite(length_m_per_internal) or length_m_per_internal <= 0:
        raise ValueError('State mesh length unit must be positive and finite')
    data = bytearray(b'vela.mesh/1')
    def integer(v): data.extend(struct.pack('<Q', v))
    def real(v): data.extend(struct.pack('<d', v))
    def text(v):
        b = v.encode('utf-8'); integer(len(b)); data.extend(b)
    def indices(values):
        integer(len(values))
        for v in values: integer(v)
    real(length_m_per_internal)
    integer(len(mesh['nodes']))
    for n in mesh['nodes']: integer(n['id']); real(n['x']); real(n['y'])
    integer(len(mesh['triangles']))
    for c in mesh['triangles']: integer(c['id']); integer(c['region_id']); indices(c['node_ids'])
    integer(len(mesh['regions']))
    for r in mesh['regions']: integer(r['id']); text(r['name']); text(r['material']); indices(r['cell_ids'])
    integer(len(mesh['contacts']))
    for c in mesh['contacts']:
        integer(c['id']); text(c['name']); integer(c['region_id']); indices(c['node_ids'])
        edges = c.get('edge_node_ids', [])
        integer(len(edges))
        for a,b in edges: integer(a); integer(b)
    return hashlib.sha256(data).hexdigest()


def fields_to_rows(fields):
    # Keep binary64 values numeric. Formatting every element as decimal here
    # only makes the predictor parse it again; CSV export owns text formatting.
    names=list(fields)
    columns=[np.asarray(fields[name],dtype='<f8').tolist() for name in names]
    return [dict(zip(['node_id']+names,[str(i),*values])) for i,values in enumerate(zip(*columns))]


def inspect(path):
    """Read a self-described archive for reporting, without qualifying its mesh.

    Solver/seed callers must use read() with an independently known mesh identity.
    This entry is for table inspectors and conversion reports only.
    """
    path = Path(path)
    if path.suffix != '.h5' or path.stat().st_size > MAX_BYTES + 1024*1024:
        raise ValueError('Invalid state archive path/size')
    with h5py.File(path,'r') as f:
        metadata = json.loads(f.attrs['metadata_json'])
        count = int(f.attrs['node_count'])
    return read(path,count,metadata['mesh_sha256'])


def read_rows(path):
    return fields_to_rows(inspect(path)[0])


def rows_to_fields(rows):
    if not rows: raise ValueError('Empty state')
    names = set(rows[0])-{'node_id'}
    if not REQUIRED <= names or not names <= UNITS.keys(): raise ValueError('Unknown state columns')
    fields = {name:np.empty(len(rows),dtype='<f8') for name in names}
    seen = set()
    for row in rows:
        node = int(row['node_id'])
        if node in seen or node < 0 or node >= len(rows) or set(row)-{'node_id'} != names:
            raise ValueError('Invalid state node order or columns')
        seen.add(node)
        for name in names: fields[name][node] = float(row[name])
    return fields


def translate(source, destination, subtract, active_nodes, expected_nodes, expected_mesh_sha256):
    if not math.isfinite(subtract): raise ValueError('Invalid frame offset')
    if Path(destination).exists(): raise FileExistsError(destination)
    fields, metadata = read(source, expected_nodes, expected_mesh_sha256)
    active = set(active_nodes)
    for i in range(expected_nodes):
        live = i in active
        for carrier in ('electron','hole'):
            ref,inc = carrier+'_qf_reference_V',carrier+'_qf_increment_V'
            if ref in fields:
                a,b=float(fields[ref][i]),float(fields[inc][i]);new=a-subtract
                fields[ref][i]=new if live else 0.
                fields[inc][i]=math.fsum((a,-subtract,-new,b)) if live else 0.
        for key in ('psi','phin','phip'):
            fields[key][i]=float(fields[key][i])-subtract if key=='psi' or live else 0.
    metadata['potential_origin_V'] += subtract
    # Saved DD bias values are expressed in the solver reference frame.
    # Physical terminal bias = stored bias + potential_origin_V.
    if 'bias_V' in metadata: metadata['bias_V'] -= subtract
    if 'contact_biases_V' in metadata:
        metadata['contact_biases_V']={k:v-subtract for k,v in metadata['contact_biases_V'].items()}
    write(destination, fields, metadata)


def validate(fields, metadata):
    required = REQUIRED if metadata.get('mode') == 'dd' else {'psi', 'phin', 'phip', 'temperature_K'}
    if not required <= fields.keys() or not fields.keys() <= UNITS.keys():
        raise ValueError('Missing or unknown state field')
    if ('electrons_m3' in fields) != ('holes_m3' in fields):
        raise ValueError('Partial density pair')
    count = len(fields['psi'])
    if count < 1 or count * len(fields) * 8 > MAX_BYTES:
        raise ValueError('Invalid or oversized dimensions')
    if metadata.get('mode') not in ('dd', 'electrothermal'):
        raise ValueError('Unknown equation mode')
    if not re.fullmatch('[0-9a-f]{64}', metadata.get('mesh_sha256', '')):
        raise ValueError('Mesh identity required')
    if type(metadata['potential_origin_V']) not in (int, float) or not math.isfinite(metadata['potential_origin_V']):
        raise ValueError('Invalid potential origin')
    if len(json.dumps(metadata, allow_nan=False).encode()) > 65536:
        raise ValueError('Metadata too large')
    if ('temperature_K' in fields) != (metadata['mode'] == 'electrothermal'):
        raise ValueError('Temperature/mode mismatch')
    if 'electron_quantum_potential_like_V' in fields and 'electron_quantum_potential_V' not in fields:
        raise ValueError('Partial quantum state')
    if len(SPLIT.intersection(fields)) not in (0, 4):
        raise ValueError('Partial split QF state')
    for name, values in fields.items():
        a = np.asarray(values, dtype='<f8')
        if a.shape != (count,) or not np.all(np.isfinite(a)):
            raise ValueError('Invalid field '+name)
    if SPLIT <= fields.keys():
        for carrier, qf in (('electron', 'phin'), ('hole', 'phip')):
            combined = np.asarray(fields[carrier+'_qf_reference_V']) + np.asarray(fields[carrier+'_qf_increment_V'])
            value = np.asarray(fields[qf])
            ceiling = 32*np.finfo(float).eps*np.maximum(1., np.maximum(abs(value), abs(combined)))
            if not np.all(np.isfinite(combined)) or np.any(abs(value-combined) > ceiling):
                raise ValueError('Inconsistent split QF state')
    return count


def write(path, fields, metadata):
    path = Path(path)
    if path.suffix != '.h5':
        raise ValueError('State path must end in .h5')
    count = validate(fields, metadata)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=path.name+'.tmp.', dir=path.parent)
    os.close(fd)
    try:
        with h5py.File(temporary, 'w') as f:
            f.attrs['schema'] = SCHEMA
            f.attrs['node_count'] = np.uint64(count)
            f.attrs['metadata_json'] = json.dumps(metadata, allow_nan=False, separators=(',', ':'))
            group = f.create_group('fields')
            for name in sorted(fields):
                d = group.create_dataset(name, data=np.asarray(fields[name], dtype='<f8'))
                d.attrs['unit'] = UNITS[name]
            f.flush()
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def read(path, expected_nodes, expected_mesh_sha256):
    path = Path(path)
    if path.suffix != '.h5' or path.stat().st_size > MAX_BYTES+2**20:
        raise ValueError('Invalid path or oversized state')
    if not re.fullmatch('[0-9a-f]{64}', expected_mesh_sha256):
        raise ValueError('Expected mesh identity required')
    def text(value):
        if isinstance(value, bytes): return value.decode('utf-8')
        if not isinstance(value, str): raise ValueError('Expected scalar string attribute')
        return value
    with h5py.File(path, 'r') as f:
        if text(f.attrs['schema']) != SCHEMA:
            raise ValueError('Unsupported schema')
        raw_count = np.asarray(f.attrs['node_count'])
        if raw_count.shape != () or raw_count.dtype.kind != 'u' or raw_count.dtype.itemsize != 8:
            raise ValueError('Node count must be scalar uint64')
        count = int(raw_count)
        if count != expected_nodes or count < 1:
            raise ValueError('Node count mismatch')
        raw = text(f.attrs['metadata_json'])
        if len(raw.encode()) > 65536: raise ValueError('Metadata too large')
        metadata = json.loads(raw)
        if metadata['mesh_sha256'] != expected_mesh_sha256:
            raise ValueError('Mesh identity mismatch')
        group = f['fields']
        if not 4 <= len(group) <= len(UNITS) or count*len(group)*8 > MAX_BYTES:
            raise ValueError('Invalid or oversized dimensions')
        fields = {}
        for name in group:
            if name not in UNITS: raise ValueError('Unknown field')
            d = group[name]
            if d.shape != (count,) or d.dtype != np.dtype('<f8') or text(d.attrs['unit']) != UNITS[name]:
                raise ValueError('Field type, dimensions or units mismatch')
            fields[name] = d[...]
    validate(fields, metadata)
    return fields, metadata
