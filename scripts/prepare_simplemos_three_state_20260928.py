"""Prepare (a) exact split, (b) physical projection, (c) M60-export projection."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'build/simplemos_ep_20260928/three_state'
BASE=ROOT/'build/outlier_analysis_20260928'
EXPORT=ROOT/'build-release/reference_tcad/simplemos_sentaurus2022/m60_tight_convergence_port_burst/state_exports/default'
REMOTE='/workspaces/simplemos-hfs-merged-20260926/ep_review_20260928'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def readcsv(p):
    with p.open(newline='',encoding='utf-8-sig') as f: return list(csv.DictReader(f))
def scalar(folder,name):
    paths=sorted((folder/'fields').glob(name+'_region*.csv')) if name=='ElectrostaticPotential' else [folder/'fields'/f'{name}_region0.csv']
    out={}
    for p in paths:
        for r in readcsv(p):
            i,v=int(r['node_id']),float(r['component0'])
            if i in out and out[i]!=v: raise ValueError(f'Conflicting regional value {p}:{i}')
            out[i]=v
    return out

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    executable=ROOT/'build/simplemos_ep_20260928/repack.exe'
    manifest=dict(binary_expected_sha256='87b934c79451110467cee69b62f6ec2e312a8aea51ab85c07bf2c0577ad05ea3',
        packer_source_sha256=sha(ROOT/'scripts/diagnostics/simplemos_split_repack_20260928.cpp'),
        packer_binary_sha256=sha(executable), conversion_script_sha256=sha(Path(__file__)),states=[])
    env=os.environ.copy();env['PATH']='D:/msys64/ucrt64/bin;'+env.get('PATH','')
    for device in ('n23','n24'):
        parent=BASE/'detailed'/f'{device}_vd_0p05_vg_001.h5'
        exports=EXPORT/f'{device}_vd_0p05_vg_0p05'
        mesh=json.loads((BASE/'inputs'/device/'mesh.json').read_text())
        native_nodes=readcsv(exports/'nodes.csv')
        assert len(native_nodes)==len(mesh['nodes'])
        for r in native_nodes:
            node=mesh['nodes'][int(r['id'])]
            assert node['id']==int(r['id']) and node['x']==float(r['x_um']) and node['y']==float(r['y_um'])
        native={field:scalar(exports,name) for field,name in (
            ('psi','ElectrostaticPotential'),('phin','eQuasiFermiPotential'),('phip','hQuasiFermiPotential'),
            ('electrons_m3','eDensity'),('holes_m3','hDensity'))}
        assert len(native['psi'])==len(mesh['nodes'])
        assert all(set(native[k])==set(native['phin']) for k in ('phip','electrons_m3','holes_m3'))
        with h5py.File(parent) as f:
            fields={k:v[()] for k,v in f['fields'].items()};meta=json.loads(f.attrs['metadata_json'])
        for arm in ('a','b','c'):
            path=OUT/f'{device}_{arm}.h5'
            if path.exists(): raise FileExistsError(path)
            shutil.copy2(parent,path)
            delta={}
            if arm!='a':
                values={k:v.copy() for k,v in fields.items()}
                if arm=='c':
                    for k,rows in native.items():
                        for i,v in rows.items(): values[k][i]=v*(1e6 if k.endswith('_m3') else 1)
                request=dict(scale=meta['packed_potential_scale_V'],rows=[
                    [float(values[k][i]) for k in ('psi','phin','phip','electron_qf_reference_V','hole_qf_reference_V')]
                    for i in range(len(values['psi']))])
                result=json.loads(subprocess.run([str(executable)],input=json.dumps(request),text=True,
                    capture_output=True,env=env,check=True).stdout)
                for block,(field,packed) in enumerate(zip(('psi','phin','phip'),('packed_psi','packed_electron_qf_increment','packed_hole_qf_increment'))):
                    physical=np.array([r['physical'][block] for r in result])
                    delta[field]=float(np.max(np.abs(physical-values[field])))
                    values[field]=physical
                    values[packed]=np.array([r['hi'][block] for r in result])
                    values[packed+'_low']=np.array([r['lo'][block] for r in result])
                    if block: values[('electron' if block==1 else 'hole')+'_qf_increment_V']=np.array([r['increment'][block] for r in result])
                with h5py.File(path,'r+') as f:
                    for k,v in values.items(): f['fields'][k][...]=v
                    derived=dict(meta)
                    derived['diagnostic_provenance']=dict(role=arm,parent_sha256=sha(parent),
                        conversion_script_sha256=manifest['conversion_script_sha256'],
                        packer_source_sha256=manifest['packer_source_sha256'],
                        source_config_sha256_is_parent_config=True,
                        native_export_hashes={p.relative_to(exports).as_posix():sha(p) for p in exports.rglob('*.csv')} if arm=='c' else {})
                    f.attrs['metadata_json']=json.dumps(derived)
            manifest['states'].append(dict(device=device,arm=arm,file=path.name,sha256=sha(path),
                parent_sha256=sha(parent),input_reconstruction_max_abs_V=delta,
                split_restore_verified=False))
            for contact in ('drain','source','substrate'):
                cfg=json.loads((BASE/'sweeps'/f'{device}_vd_0p05'/'config.json').read_text())
                for k in ('sweep','initial_state_file','output_csv'):cfg.pop(k,None)
                cfg.update(simulation_type='terminal_current_functional_probe',state_format='hdf5',
                    state_file=f'{REMOTE}/{path.name}',contact=contact)
                for c in cfg['contacts']:c['bias']=.05 if c['name'] in ('gate','drain') else 0.
                if contact=='drain':cfg['residual_output_csv']=f'{REMOTE}/{device}_{arm}_residual.csv'
                (OUT/f'{device}_{arm}_{contact}.json').write_text(json.dumps(cfg,indent=2)+'\n')
    manifest['files']={p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(states=len(manifest['states']),reconstruction=[s['input_reconstruction_max_abs_V'] for s in manifest['states']])))

if __name__=='__main__':main()
