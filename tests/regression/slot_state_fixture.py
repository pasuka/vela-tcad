"""Independent tiny mesh/state fixtures for SLOT restart regressions."""
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
import state_archive


def mesh(count=3):
    return dict(nodes=[dict(id=i, x=float(i), y=0.) for i in range(count)],
                triangles=[], regions=[], contacts=[])


def complete_mesh(value):
    for region in value['regions']:
        region.setdefault('name', region['material'])
        region.setdefault('cell_ids', [c['id'] for c in value['triangles'] if c['region_id'] == region['id']])
    for index, contact in enumerate(value['contacts']):
        contact.setdefault('id', index)
        contact.setdefault('region_id', 0)
    return value


def write(path, rows, geometry=None, origin=0., **metadata):
    path = Path(path)
    geometry = complete_mesh(geometry if geometry is not None else mesh(len(rows)))
    path.parent.mkdir(parents=True, exist_ok=True)
    (path.parent / 'mesh.json').write_text(json.dumps(geometry))
    state_archive.write(path, state_archive.rows_to_fields(rows), dict(mode='dd',
        mesh_sha256=state_archive.mesh_identity(geometry, 1.), potential_origin_V=origin, **metadata))


def read(path):
    return state_archive.read_rows(path)
