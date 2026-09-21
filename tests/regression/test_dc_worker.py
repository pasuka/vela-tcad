"""CLI isolation: reject a request, then solve from the next explicit input."""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

RUNNER = None


class DCWorkerTest(unittest.TestCase):
    def test_linear_context_reuses_analysis_and_invalidates_on_failure_and_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cells = [[0,1,4],[1,2,4],[2,3,4],[3,0,4]]
            mesh = dict(
                nodes=[dict(id=i,x=x*1e-6,y=y*1e-6) for i,(x,y) in enumerate(
                    ((0,0),(1,0),(1,1),(0,1),(.5,.5)))],
                triangles=[dict(id=i,region_id=i//2,node_ids=n) for i,n in enumerate(cells)],
                regions=[dict(id=0,name='n',material='Si',cell_ids=[0,1]),
                         dict(id=1,name='p',material='Si',cell_ids=[2,3])],
                contacts=[dict(id=0,name='anode',region_id=1,node_ids=[0,3]),
                          dict(id=1,name='cathode',region_id=0,node_ids=[1,2])])
            (root/'mesh.json').write_text(json.dumps(mesh))
            configs=[]
            for i in range(4):
                density=1e21 if i<3 else 5e20
                cfg=dict(simulation_type='dc_sweep',mesh_file='mesh.json',output_csv=f'curve{i}.csv',
                    doping=[dict(region='n',donors=density,acceptors=0),dict(region='p',donors=0,acceptors=density)],
                    contacts=[dict(name='anode',bias=0),dict(name='cathode',bias=0)],
                    solver=dict(method='newton',max_iter=30,reltol=1e-7,abstol=1e-20,
                                performance_profiling=dict(enabled=True,json_file=f'profile{i}.json')),
                    sweep=dict(contact='anode',current_contact='anode',start=.05,stop=.05,step=.1,write_vtk=False))
                path=root/f'config{i}.json';path.write_text(json.dumps(cfg));configs.append(path)
            expected=[]
            for path in configs:
                subprocess.run([str(RUNNER),'--config',str(path)],check=True,capture_output=True,timeout=30)
                expected.append((root/json.loads(path.read_text())['output_csv']).read_text())
            requests=[json.dumps(dict(id=i,config=str(path))) for i,path in enumerate(configs)]
            requests.insert(2,'malformed')
            requests.append('{"shutdown":true}')
            result=subprocess.run([str(RUNNER),'--dc-worker'],input='\n'.join(requests)+'\n',
                                  capture_output=True,text=True,timeout=60)
            self.assertEqual(result.returncode,0,result.stderr)
            responses=[json.loads(line) for line in result.stdout.splitlines()]
            self.assertIn('error',responses[2])
            self.assertEqual([r['returncode'] for i,r in enumerate(responses) if i!=2],[0]*4)
            for i in range(4):
                self.assertEqual((root/f'curve{i}.csv').read_text(),expected[i])
                counters=json.loads((root/f'profile{i}.json').read_text())['counters']
                self.assertGreater(counters.get('linear.factorize_calls',0),0)
                key='hits' if i==1 else 'misses'
                self.assertEqual(counters.get('dc.linear_context.'+key),1)
                if i==1:self.assertEqual(counters.get('linear.analyze_calls',0),0)
            # The explicit control must still rebuild across requests after
            # enabling reuse by default; a baseline must not silently reuse.
            result=subprocess.run([str(RUNNER),'--dc-worker-no-linear-reuse'],input='\n'.join(requests)+'\n',
                                  capture_output=True,text=True,timeout=60)
            self.assertEqual(result.returncode,0,result.stderr)
            counters=json.loads((root/'profile1.json').read_text())['counters']
            self.assertGreater(counters.get('linear.analyze_calls',0),0)
            self.assertNotIn('dc.linear_context.hits',counters)

    def test_failure_does_not_contaminate_following_requests(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mesh = dict(
                nodes=[dict(id=i, x=x*1e-6, y=y*1e-6) for i,(x,y) in enumerate(((0,0),(1,0),(1,1),(0,1)))],
                triangles=[dict(id=0,region_id=0,node_ids=[0,1,2]),dict(id=1,region_id=1,node_ids=[0,2,3])],
                regions=[dict(id=0,name='n',material='Si',cell_ids=[0]),dict(id=1,name='p',material='Si',cell_ids=[1])],
                contacts=[dict(id=0,name='anode',region_id=1,node_ids=[0,3]),dict(id=1,name='cathode',region_id=0,node_ids=[1,2])])
            (root/'mesh.json').write_text(json.dumps(mesh))
            configs=[];expected=[]
            for i,bias in enumerate((0.,.1,0.)):
                cfg=dict(simulation_type='dc_sweep',mesh_file='mesh.json',output_csv=f'curve{i}.csv',
                         doping=[dict(region='n',donors=1e23,acceptors=0),dict(region='p',donors=0,acceptors=1e23)],
                         contacts=[dict(name='anode',bias=0),dict(name='cathode',bias=0)],
                         solver=dict(max_iter=80,reltol=1e-5,damping_psi=.5,performance_profiling=dict(enabled=True,json_file=f'profile{i}.json')),
                         sweep=dict(contact='anode',current_contact='anode',start=bias,stop=bias,step=.1,write_vtk=False))
                path=root/f'config{i}.json';path.write_text(json.dumps(cfg));configs.append(path)
                subprocess.run([str(RUNNER),'--config',str(path)],check=True,capture_output=True,timeout=30)
                expected.append((root/f'curve{i}.csv').read_text())
            bad=json.loads(configs[0].read_text());bad['sweep']['initial_state_file']='missing.csv'
            bad['output_csv']='bad.csv'
            bad['solver']['performance_profiling']['json_file']='bad_profile.json'
            bad_path=root/'bad.json';bad_path.write_text(json.dumps(bad))
            inputs=['malformed',json.dumps(dict(id=0,config=str(configs[0]))),json.dumps(dict(id=1,config=str(bad_path)))]
            inputs += [json.dumps(dict(id=i+2,config=str(path))) for i,path in enumerate(configs[1:])]
            inputs.append('{"shutdown":true}')
            result=subprocess.run([str(RUNNER),'--dc-worker'],input='\n'.join(inputs)+'\n',capture_output=True,text=True,timeout=60)
            self.assertEqual(result.returncode,0,result.stderr)
            responses=[json.loads(line) for line in result.stdout.splitlines()]
            self.assertEqual(len(responses),5)
            self.assertIn('error',responses[0])
            self.assertEqual([r.get('id') for r in responses[1:]],[0,1,2,3])
            self.assertEqual([r['returncode'] for r in responses[1:]],[0,1,0,0])
            for i in range(3):
                self.assertEqual((root/f'curve{i}.csv').read_text(),expected[i])
            counters=json.loads((root/'profile2.json').read_text())['counters']
            self.assertEqual(counters['dc.prepared_inputs.hits'],1)


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--runner',type=Path,required=True)
    args,remaining=parser.parse_known_args();RUNNER=args.runner.resolve()
    unittest.main(argv=[__file__]+remaining)
