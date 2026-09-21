"""Extract D0 nodal temperatures using an exactly verified existing TDR mapping.

Requires h5py and a C++ importer export of the baseline TDR. Every input must
have identical non-state geometry datasets; no positional mapping is guessed.
"""
import argparse,csv,hashlib,json,math
from pathlib import Path
import h5py
import numpy as np
from run_templates_ldmos_sentaurus_vm import normalize_plt_files


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def geometry_hash(handle):
    result=hashlib.sha256()
    def visit(name,item):
        if isinstance(item,h5py.Dataset) and not any(p.startswith('state_') for p in name.split('/')):
            array=item[()]
            if array.dtype.hasobject:raise ValueError('Unsupported variable-length geometry dataset')
            result.update(name.encode());result.update(str(array.dtype).encode());result.update(str(array.shape).encode());result.update(array.tobytes())
    handle['collection/geometry_0'].visititems(visit)
    return result.hexdigest()
def values(handle,name,region):
    found=[]
    for item in handle['collection/geometry_0/state_0'].values():
        attr=item.attrs;label=attr.get('name',b'')
        if isinstance(label,bytes):label=label.decode()
        if label==name and int(attr.get('region',-1))==region:
            if int(attr['location type'])!=0 or float(attr['conversion factor'])!=1.:raise ValueError('Unsupported nodal field metadata')
            found.append(item['values'][()].astype(float))
    if len(found)!=1:raise ValueError(f'Expected one {name} region {region}')
    return found[0]
def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('baseline-tdr','baseline-export','mesh','raw','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    mesh=json.loads(a.mesh.read_text(encoding='utf-8'));count=len(mesh['nodes']);mapping={};sources={str(a.baseline_tdr):sha(a.baseline_tdr),str(a.mesh):sha(a.mesh)}
    with h5py.File(a.baseline_tdr,'r') as f:
        signature=geometry_hash(f)
        for region in (0,1,2):
            path=a.baseline_export/'fields'/f'LatticeTemperature_region{region}.csv';rows=list(csv.DictReader(path.open(encoding='utf-8')))
            ids=[int(r['node_id']) for r in rows];expected=np.array([float(r['component0']) for r in rows])
            if len(set(ids))!=len(ids) or any(not 0<=i<count for i in ids):raise ValueError('Invalid baseline node map')
            if not np.array_equal(expected,values(f,'LatticeTemperature',region)):raise ValueError('Baseline importer mapping is not exact')
            mapping[region]=ids;sources[str(path)]=sha(path)
    contacts={c['name']:c['node_ids'] for c in mesh['contacts']};index=[]
    for gate in (4,8):
        files=sorted(a.raw.glob(f'field_vg{gate}_*_des.tdr'))
        if len(files)!=31:raise ValueError('Expected 31 exact native fields per gate')
        for point,path in enumerate(files):
            merged={}
            with h5py.File(path,'r') as f:
                if geometry_hash(f)!=signature:raise ValueError(f'Geometry changed: {path}')
                for region,ids in mapping.items():
                    data=values(f,'LatticeTemperature',region)
                    if len(ids)!=len(data):raise ValueError('Field length mismatch')
                    for i,t in zip(ids,data):
                        if i in merged and abs(merged[i]-t)>1e-9:raise ValueError('Shared temperature mismatch')
                        merged[i]=float(t)
                fn=dict(zip(mapping[0],values(f,'eQuasiFermiPotential',0)))
                bias=float(fn[contacts['drain'][0]]-fn[contacts['source'][0]])
                if abs(bias-point*40/30)>1e-9:raise ValueError(f'Unexpected exact field bias {bias}')
            temperature=[merged[i] for i in range(count)]
            if not all(math.isfinite(t) and t>0. for t in temperature):raise ValueError('Invalid temperature')
            output=a.output/f'vg{gate}_{point:02d}.json'
            output.write_text(json.dumps(dict(node_id=list(range(count)),temperature_K=temperature)),encoding='utf-8')
            index.append(dict(gate_V=gate,point_index=point,bias_V=bias,peak_K=max(temperature),temperature_file=str(output.resolve()),tdr_sha256=sha(path),temperature_sha256=sha(output)))
    normalized=normalize_plt_files(a.raw,a.output/'normalized')
    (a.output/'manifest.json').write_text(json.dumps(dict(scope='62 native D0 temperatures; exact baseline importer mapping and geometry checks',geometry_sha256=signature,sources_sha256=sources,h5py_version=h5py.__version__,hdf5_runtime=h5py.version.hdf5_version,fields=index,normalized=normalized),indent=2),encoding='utf-8')
    print(json.dumps(dict(fields=len(index),nodes=count,geometry_sha256=signature)))
if __name__=='__main__':main()
