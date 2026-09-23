"""One-time repository migration, deliberately separate from production readers."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import state_archive


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def migrate(source, mesh_path, output, scaling, origin):
    if output.exists(): raise FileExistsError(output)
    if source.suffix == '.csv':
        with source.open(encoding='utf-8-sig',newline='') as stream:
            data = list(csv.DictReader(stream))
    elif source.suffix == '.vds':
        import dd_state_binary  # Migration-only legacy decoder.
        data = dd_state_binary.read(source)
    else:
        raise ValueError('Migration input must be CSV or VDS1')
    mesh = json.loads(mesh_path.read_text(encoding='utf-8-sig'))
    fields = state_archive.rows_to_fields(data)
    if len(data) != len(mesh['nodes']): raise ValueError('Mesh/state node count differs')
    metadata = dict(mode='dd',mesh_sha256=state_archive.mesh_identity(mesh,1e-6 if scaling=='unit_scaling' else 1.),
                    potential_origin_V=origin,source_state_sha256=sha(source))
    state_archive.write(output,fields,metadata)
    restored,actual_metadata = state_archive.read(output,len(data),metadata['mesh_sha256'])
    if actual_metadata != metadata: raise ValueError('Metadata changed')
    hashes = {}
    for name in fields:
        before = np.asarray(fields[name],dtype='<f8').tobytes()
        after = np.asarray(restored[name],dtype='<f8').tobytes()
        if before != after: raise ValueError('Field changed: '+name)
        hashes[name] = hashlib.sha256(after).hexdigest()
    return dict(schema='vela.state_migration.conversion/1',source=str(source),output=str(output),
                source_sha256=sha(source),output_sha256=sha(output),mesh_file_sha256=sha(mesh_path),
                metadata=metadata,node_count=len(data),field_sha256=hashes,all_decoded_fields_equal=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--mesh',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--report',type=Path,required=True)
    p.add_argument('--scaling',choices=('legacy_si','unit_scaling'),required=True)
    p.add_argument('--potential-origin',type=float,required=True)
    a=p.parse_args()
    if a.report.exists(): raise FileExistsError(a.report)
    result=migrate(a.source,a.mesh,a.output,a.scaling,a.potential_origin)
    a.report.parent.mkdir(parents=True,exist_ok=True)
    with a.report.open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2)
    print(json.dumps(dict(node_count=result['node_count'],fields=len(result['field_sha256']),equal=True)))
