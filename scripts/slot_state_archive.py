"""SLOT DD restart checks against independent deck/mesh identity."""
from pathlib import Path
import hashlib
import json
import math

import numpy as np
import state_archive as archive


def identity(document, directory):
    mesh_path = (Path(directory) / document['mesh_file']).resolve()
    mesh = json.loads(mesh_path.read_text(encoding='utf-8-sig'))
    mode = document.get('scaling', {}).get('mode', 'legacy_si')
    if mode not in ('legacy_si', 'unit_scaling'):
        raise ValueError('Unsupported SLOT state scaling')
    origin = document.get('potential_origin_V', 0.)
    if type(origin) not in (int, float) or not math.isfinite(origin):
        raise ValueError('Invalid SLOT potential origin')
    return len(mesh['nodes']), archive.mesh_identity(mesh, 1e-6 if mode == 'unit_scaling' else 1.), origin


def read(path, document, directory):
    count, mesh_hash, origin = identity(document, directory)
    fields, metadata = archive.read(path, count, mesh_hash)
    if metadata['mode'] != 'dd' or metadata['potential_origin_V'] != origin:
        raise ValueError('SLOT restart mode/potential origin mismatch')
    if any(np.any(fields[key] < 0) for key in ('electrons_m3', 'holes_m3')):
        raise ValueError('Negative SLOT carrier density')
    return fields, metadata


def digest(path, document, directory):
    """Exact decoded state identity, independent of container/provenance bytes."""
    fields, metadata = read(path, document, directory)
    metadata['potential_origin_V'] = float(metadata['potential_origin_V'])
    result = hashlib.sha256(json.dumps({key: metadata[key] for key in
        ('mode', 'mesh_sha256', 'potential_origin_V')}, sort_keys=True).encode())
    for name in sorted(fields):
        result.update(name.encode() + b'\0')
        result.update(np.asarray(fields[name], dtype='<f8').tobytes())
    return result.hexdigest()


def write_candidate(path, fields, metadata, document, parents, bias, role):
    path = Path(path).resolve()
    if path.exists():
        raise FileExistsError(path)
    if not math.isfinite(bias):
        raise ValueError('Nonfinite candidate bias')
    metadata = dict(metadata)
    # Parent provenance remains recorded, but old accepted-state labels must not
    # imply that a predictor or equilibrium candidate has passed device gates.
    metadata.update(state_role=role, qualified=False, bias_V=bias,
                    parent_state_sha256=[hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in parents],
                    source_config_sha256=hashlib.sha256(json.dumps(document, sort_keys=True,
                        allow_nan=False).encode()).hexdigest())
    driven = document.get('sweep', {}).get('contact')
    if driven:
        metadata['bias_contact'] = driven
        metadata['contact_biases_V'] = {c['name']: bias if c['name'] == driven else c.get('bias', 0.)
                                        for c in document.get('contacts', [])}
    archive.write(path, fields, metadata)
