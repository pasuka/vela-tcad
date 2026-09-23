"""Real Python -> C++ -> Python restart, format conversion and failure isolation."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import dd_state_public as public
import dd_state_binary as binary
import numpy as np
import h5py
RUNNER=None


class PublicStateTest(unittest.TestCase):
    def test_cross_language_roundtrip_and_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            mesh=dict(nodes=[dict(id=i,x=x*1e-6,y=y*1e-6) for i,(x,y) in enumerate(((0,0),(1,0),(1,1),(0,1)))],
                triangles=[dict(id=0,region_id=0,node_ids=[0,1,2]),dict(id=1,region_id=1,node_ids=[0,2,3])],
                regions=[dict(id=0,name='n',material='Si',cell_ids=[0]),dict(id=1,name='p',material='Si',cell_ids=[1])],
                contacts=[dict(id=0,name='anode',region_id=1,node_ids=[0,3]),dict(id=1,name='cathode',region_id=0,node_ids=[1,2])])
            (root/'mesh.json').write_text(json.dumps(mesh))
            state=[dict(node_id=str(i),psi=.11*(i+1),phin=1e-17,phip=-1e-17,electrons_m3=1e10,holes_m3=1e10,
                electron_quantum_potential_V=0.,electron_qf_increment_V=1e-17,hole_qf_increment_V=-1e-17,
                electron_qf_reference_V=0.,hole_qf_reference_V=0.) for i in range(4)]
            extensions=('.h5','.vfb','.npy')
            requests=[];expected=[]
            def request(seed,extension,valid):
                i=len(requests)
                cfg=dict(simulation_type='dc_sweep',mesh_file='mesh.json',output_csv=f'curve{i}.csv',
                    doping=[dict(region='n',donors=1e23,acceptors=0),dict(region='p',donors=0,acceptors=1e23)],
                    contacts=[dict(name='anode',bias=0),dict(name='cathode',bias=0)],
                    solver=dict(method='frozen_state',impact_ionization='none'),
                    sweep=dict(contact='anode',current_contact='anode',start=0.,stop=0.,step=.1,write_vtk=False,
                        initial_state_file=seed.name,write_state_file=f'out{i}{extension}'))
                path=root/f'config{i}.json';path.write_text(json.dumps(cfg))
                requests.append(json.dumps(dict(id=i,config=str(path))));expected.append((valid,root/f'out{i}{extension}'))
            for ext in extensions:
                seed=root/('seed'+ext);public.write(seed,state)
                self.assertEqual(public.read(seed),state)
                for out in extensions:request(seed,out,True)
                short=root/('short'+ext);public.write(short,state[:3]);request(short,ext,False)
                bad=root/('truncated'+ext);bad.write_bytes(seed.read_bytes()[:8]);request(bad,ext,False)
                with self.assertRaises(Exception):public.read(bad)
                # Nonfinite values must be rejected after decoding, not silently accepted.
                nonfinite=root/('nan'+ext)
                if ext=='.h5':
                    public.write(nonfinite,state)
                    with h5py.File(nonfinite,'r+') as f:f['psi'][0]=float('nan')
                elif ext=='.npy':
                    data=np.load(seed,allow_pickle=False);data['psi'][0]=float('nan');np.save(nonfinite,data)
                else:
                    data=bytearray(seed.read_bytes())
                    obj=public.fb.DDState.GetRootAs(data);offset=obj._tab.Vector(obj._tab.Offset(10))
                    import struct
                    struct.pack_into('<d',data,offset,float('nan'));nonfinite.write_bytes(data)
                with self.assertRaises(ValueError):public.read(nonfinite)
                request(nonfinite,ext,False)
                request(seed,ext,True)  # Failure must not poison next worker request.
            wrong=root/'units.h5';public.write(wrong,state)
            with h5py.File(wrong,'r+') as f:f.attrs['units']='cm^-3'
            with self.assertRaises(ValueError):public.read(wrong)
            request(wrong,'.h5',False)
            records=np.load(root/'seed.npy',allow_pickle=False)
            wrong=root/'float32.npy';np.save(wrong,records.astype([(n,'<f4') for n in records.dtype.names]))
            with self.assertRaises(ValueError):public.read(wrong)
            request(wrong,'.npy',False)
            # A C++ writer's output is then used as a C++ input (restart, not only export).
            request(root/'out0.h5','.npy',True)
            request(root/'out1.vfb','.h5',True)
            request(root/'out2.npy','.vfb',True)
            requests.append('{"shutdown":true}')
            r=subprocess.run([str(RUNNER),'--dc-worker'],input='\n'.join(requests)+'\n',capture_output=True,text=True,timeout=90)
            self.assertEqual(r.returncode,0,r.stderr)
            responses=[json.loads(line) for line in r.stdout.splitlines()]
            self.assertEqual(len(responses),len(expected),r.stdout)
            for result,(valid,path) in zip(responses,expected):
                with self.subTest(response=result,path=path):
                    if valid:
                        self.assertEqual(result.get('returncode'),0,result)
                        self.assertEqual(public.read(path),state)
                    else:self.assertTrue('error' in result or result.get('returncode',0)!=0,result)

    def test_python_optional_fields_and_subnormals(self):
        with tempfile.TemporaryDirectory() as tmp:
            for flags in (0,1,3,4,5,7):
                rows=[dict(node_id='0',**{n:0. for n in binary.fields(flags)})]
                rows[0]['psi']=5e-324
                if flags&4:
                    rows[0]['electron_qf_increment_V']=1e-17
                    rows[0]['phin']=1e-17
                for ext in ('.h5','.vfb','.npy'):
                    p=Path(tmp)/(str(flags)+ext);public.write(p,rows)
                    self.assertEqual(public.read(p),binary.decode(binary.encode(rows)))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--runner',type=Path,required=True)
    a,rest=p.parse_known_args();RUNNER=a.runner.resolve();unittest.main(argv=[__file__]+rest)
