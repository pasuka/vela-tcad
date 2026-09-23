"""Python/C++ VDS1 interoperability through the real sequential worker."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import dd_state_binary as binary
RUNNER=None

class BinaryWorkerTest(unittest.TestCase):
    def test_python_state_crosses_cpp_worker_and_bad_input_is_isolated(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            mesh=dict(nodes=[dict(id=i,x=x*1e-6,y=y*1e-6) for i,(x,y) in enumerate(((0,0),(1,0),(1,1),(0,1)))],
                triangles=[dict(id=0,region_id=0,node_ids=[0,1,2]),dict(id=1,region_id=1,node_ids=[0,2,3])],
                regions=[dict(id=0,name='n',material='Si',cell_ids=[0]),dict(id=1,name='p',material='Si',cell_ids=[1])],
                contacts=[dict(id=0,name='anode',region_id=1,node_ids=[0,3]),dict(id=1,name='cathode',region_id=0,node_ids=[1,2])])
            (root/'mesh.json').write_text(json.dumps(mesh))
            state=[dict(node_id=str(i),psi=.11*(i+1),phin=1e-17,phip=-1e-17,electrons_m3=1e10,holes_m3=1e10,
                electron_qf_increment_V=1e-17,hole_qf_increment_V=-1e-17,
                electron_qf_reference_V=0.,hole_qf_reference_V=0.) for i in range(4)]
            binary.write(root/'seed.vds',state);(root/'bad.vds').write_bytes(b'VELADS01bad')
            configs=[]
            for i,seed in enumerate(('seed.vds','bad.vds','seed.vds')):
                cfg=dict(simulation_type='dc_sweep',mesh_file='mesh.json',output_csv=f'curve{i}.csv',
                    doping=[dict(region='n',donors=1e23,acceptors=0),dict(region='p',donors=0,acceptors=1e23)],
                    contacts=[dict(name='anode',bias=0),dict(name='cathode',bias=0)],
                    solver=dict(method='frozen_state',impact_ionization='none'),
                    sweep=dict(contact='anode',current_contact='anode',start=0.,stop=0.,step=.1,write_vtk=False,
                        initial_state_file=seed,write_state_file=f'out{i}.vds'))
                path=root/f'config{i}.json';path.write_text(json.dumps(cfg));configs.append(path)
            requests=[json.dumps(dict(id=i,config=str(p))) for i,p in enumerate(configs)]+['{"shutdown":true}']
            requests.pop()
            cfg['sweep'].update(initial_state_file=str(root/'virtual.vds'),write_state_file=str(root/'memory_out.vds'))
            virtual=root/'virtual.json';virtual.write_text(json.dumps(cfg))
            payload=binary.encode(state).hex()
            extras=[dict(initial_state_vds_hex=payload),{},dict(initial_state_vds_hex='zz'),
                    dict(initial_state_vds_hex=payload[:-2]),dict(initial_state_vds_hex=payload)]
            requests.extend(json.dumps(dict(id=i+3,config=str(virtual),**extra)) for i,extra in enumerate(extras))
            requests.append('{"shutdown":true}')
            r=subprocess.run([str(RUNNER),'--dc-worker'],input='\n'.join(requests)+'\n',capture_output=True,text=True,timeout=60)
            self.assertEqual(r.returncode,0,r.stderr);responses=[json.loads(line) for line in r.stdout.splitlines()]
            self.assertEqual(responses[0]['returncode'],0,responses[0]);self.assertIn('error',responses[1])
            self.assertEqual(responses[2]['returncode'],0,responses[2])
            for i in (3,7):self.assertEqual(responses[i]['returncode'],0,responses[i])
            for i in (4,5,6):self.assertIn('error',responses[i])
            self.assertFalse((root/'virtual.vds').exists())
            loaded=binary.read(root/'memory_out.vds')
            for a,c in zip(state,loaded):
                for key in a:self.assertEqual(a[key],c[key],key)
            for i in (0,2):
                loaded=binary.read(root/f'out{i}.vds')
                for a,c in zip(state,loaded):
                    for key in a:self.assertEqual(a[key],c[key],key)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--runner',type=Path,required=True)
    a,rest=p.parse_known_args();RUNNER=a.runner.resolve();unittest.main(argv=[__file__]+rest)
