"""Verify native Load/Plot identity, then audit fixed and A-reclosed states.

Field errors are diagnostic; device acceptance remains in the point verifier.
Requires h5py for geometry/state identity and the C++ importer for node mapping.
"""
import argparse,hashlib,json,os,subprocess,sys
from pathlib import Path
import h5py
import numpy as np
from extract_templates_ldmos_d0_temperature import geometry_hash,values

def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('native-raw','native-bundle','point-results','importer','probe','parameters','output'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--gates',nargs='+',type=int,choices=(4,8),default=[4,8])
    p.add_argument('--indices',nargs='+',type=int,choices=(0,1,10,30),default=[0,1,10,30])
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    plan=read(a.native_bundle/'manifest.json');points=read(a.point_results)
    if points['status']!='pass':raise ValueError('Complete qualified A representative points required')
    for name,digest in plan['bundle_sha256'].items():
        if sha(a.native_bundle/name)!=digest:raise ValueError('Changed native bundle '+name)
    result=dict(status='running',scope=__doc__,cases=[],h5py_version=h5py.__version__,hdf5_runtime=h5py.version.hdf5_version,
                source_sha256={str(p):sha(p) for p in (a.point_results,a.importer,a.probe,a.parameters,a.native_bundle/'manifest.json',Path(__file__))})
    def save():(a.output/'summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    env=dict(os.environ)
    if os.name=='nt':env['PATH']='D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+env['PATH']
    for case in plan['cases']:
        if case['gate_V'] not in a.gates or case['index'] not in a.indices:continue
        name=case['name'];seed=a.native_bundle/case['seed'];tdr=a.native_raw/f'{name}_des.tdr';identity={}
        with h5py.File(seed,'r') as old,h5py.File(tdr,'r') as new:
            if geometry_hash(old)!=geometry_hash(new):raise ValueError('Native geometry changed')
            for field,regions,limit in [('LatticeTemperature',(0,1,2),1e-9),('ElectrostaticPotential',(0,1,2),1e-12),('eQuasiFermiPotential',(0,),1e-12),('hQuasiFermiPotential',(0,),1e-12)]:
                maximum=0.
                for region in regions:
                    x,y=values(old,field,region),values(new,field,region)
                    if x.shape!=y.shape or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):raise ValueError('Invalid native state')
                    maximum=max(maximum,float(np.max(np.abs(x-y))))
                if maximum>limit:raise ValueError(f'Load/Plot changed {name} {field}: {maximum}')
                identity[field]=maximum
        directory=a.output/name;directory.mkdir();export=directory/'export'
        with (directory/'import.log').open('w',encoding='utf-8') as log:
            subprocess.run([str(a.importer.resolve()),'--tdr',str(tdr.resolve()),'--export-dir',str(export.resolve())],stdout=log,stderr=subprocess.STDOUT,env=env,check=True)
        point=next(row for row in points['points'] if row['gate']==case['gate_V'] and row['index']==case['index'])
        item=dict(**case,native_state_max_difference=identity,A_updates=point['updates'],A_pass=point['pass_gate'])
        for mode in ('fixed_native','reclosed_A'):
            target=directory/mode
            command=[sys.executable,str(Path(__file__).with_name('audit_templates_ldmos_thermal_physics.py')),
                     '--export',str(export),'--probe',str(a.probe),'--output',str(target),'--parameters',str(a.parameters),'--require-auger']
            if mode=='reclosed_A':command+=['--vela-state',point['result']]
            with (directory/f'{mode}.log').open('w',encoding='utf-8') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env,check=True)
            item[mode]=read(target/'summary.json')['results']
        result['cases'].append(item)
        for path in (seed,tdr,Path(point['result'])):result['source_sha256'][str(path)]=sha(path)
        save();print(json.dumps(dict(case=name,status='completed')),flush=True)
    result['status']='completed';save()

if __name__=='__main__':main()
