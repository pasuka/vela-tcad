"""Production electrothermal checkpoint, path and gate regression."""
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


@unittest.skipUnless(os.environ.get('VELA_RUNNER'), 'Matching build runner required')
class ElectrothermalRunnerTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.runner = os.environ['VELA_RUNNER']
        mesh = dict(regions=[dict(id=0,name='silicon',material='Silicon',cell_ids=[0,1])],nodes=[dict(id=i, x=x, y=y) for i, (x,y) in enumerate(((0,0),(1,0),(1,1),(0,1)))],
                    triangles=[dict(id=0, region_id=0, node_ids=[0,1,2]), dict(id=1, region_id=0, node_ids=[0,2,3])],
                    contacts=[dict(id=i,region_id=0,name=name, node_ids=[i]) for i,name in enumerate(('source','drain','substrate','gate'))])
        self.write('mesh.json', mesh)
        self.input = dict(mesh_file='mesh.json', coordinate_to_metres=1.,
            region_conductivity=[dict(region_id=0,model='constant',value_W_per_m_K=1.)], thermodes=[],
            silicon_area_m2=[1/3,1/6,1/3,1/6], fixed_charge_C_per_m=[0.]*4,
            edge_geometry=[dict(nodes=[a,b],poisson_F_per_m=1e-10,transport_weight=0.) for a,b in ((0,1),(1,2),(0,2),(2,3),(0,3))],
            donors_m3=[0.]*4, acceptors_m3=[0.]*4, mobility_SI=dict(model='constant'),
            boundaries=[dict(node=i,kind=k,value=v) for i in range(4) for k,v in (('psi',0.),('fn',0.),('fp',0.),('temperature',300.))],
            state_interleaved=[0.,0.,0.,300.]*4, referenced_state_interleaved=[0.,0.,0.,300.]*4,
            electron_qf_reference_V=[0.]*4, hole_qf_reference_V=[0.]*4,
            electrical_gate_solver=dict(carrier_row_convergence=dict(mode='enforce',eps_row=1e-8),
                block_absolute_convergence=dict(mode='enforce',psi_residual_ceiling=1e-8,electron_residual_ceiling=1e-11,hole_residual_ceiling=3e-10)))
        self.write('input.json',self.input)
        self.deck=dict(simulation_type='electrothermal_dc_sweep', input_file='input.json',output_directory='output',runtime_log=dict(enabled=False),
                       sweep=dict(bias_points_V=[0.],initial_step_V=.1,minimum_step_V=.025,maximum_step_V=.1,max_newton=2,growth_newton=2,predictor='none'))
    def write(self,name,data):
        (self.root/name).write_text(json.dumps(data),encoding='utf-8')
    def run_deck(self):
        self.write('deck.json',self.deck)
        return subprocess.run([self.runner,'--config',str(self.root/'deck.json')],cwd=self.root.parent,
                              stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    def test_zero_point_and_completed_restart_preserve_result(self):
        run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
        result=self.root/'output/step_0000/output.json';before=result.read_bytes()
        self.deck['resume']=True
        run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
        self.assertEqual(before,result.read_bytes())
        self.assertTrue((self.root/'output/curve.csv').exists())
        with (self.root/'output/curve.csv').open() as stream:
            rows=list(csv.DictReader(stream))
        self.assertEqual(len(rows),1)
        self.assertEqual(float(rows[0]['current_total_A_per_um']),0.)
        self.assertEqual(float(rows[0]['peak_temperature_K']),300.)
        # Volume accumulation can round the uniform mean by one ULP.
        self.assertAlmostEqual(float(rows[0]['mean_temperature_K']),300.,delta=1e-12)
    def test_resume_rejects_modified_geometry(self):
        self.assertEqual(self.run_deck().returncode,0)
        mesh=json.loads((self.root/'mesh.json').read_text());mesh['nodes'][0]['x']=.1;self.write('mesh.json',mesh)
        self.deck['resume']=True
        run=self.run_deck();self.assertNotEqual(run.returncode,0);self.assertIn('differs from checkpoint',run.stderr)
    def test_pause_then_failed_heat_gate_keeps_accepted_zero(self):
        self.deck['sweep']['bias_points_V']=[0.,.1];self.deck['pause_after_attempts']=1
        run=self.run_deck();self.assertEqual(run.returncode,1,run.stderr)
        ledger=json.loads((self.root/'output/ledger.json').read_text())
        self.assertEqual(ledger['status'],'stopped_at_checkpoint')
        accepted=Path(ledger['accepted_result']).read_bytes()
        self.deck.update(resume=True,pause_after_attempts=0)
        run=self.run_deck();self.assertEqual(run.returncode,1,run.stderr)
        ledger=json.loads((self.root/'output/ledger.json').read_text())
        self.assertEqual(ledger['status'],'failed');self.assertEqual(ledger['accepted_bias_V'],0.)
        self.assertEqual(Path(ledger['accepted_result']).read_bytes(),accepted)
        self.assertIn('heat_balance',ledger['runs'][-1]['gate']['reasons'])
    def test_disabling_row_gate_is_rejected_before_output(self):
        self.input['electrical_gate_solver']['carrier_row_convergence']['mode']='off';self.write('input.json',self.input)
        run=self.run_deck();self.assertNotEqual(run.returncode,0)
        self.assertFalse((self.root/'output').exists())
    def test_voltage_update_limit_changes_path_but_preserves_converged_state(self):
        self.deck['sweep']['max_newton']=10
        for boundary in self.input['boundaries']:
            if boundary['kind']=='psi':boundary['value']=.5
        counts=[];states=[]
        for limit in (.2,.4):
            self.deck['output_directory']='limit_'+str(limit)
            self.input['diagnostic_voltage_update_limit_V']=limit
            self.write('input.json',self.input)
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            ledger=json.loads((self.root/self.deck['output_directory']/'ledger.json').read_text())
            result=json.loads(Path(ledger['accepted_result']).read_text())
            self.assertTrue(ledger['runs'][0]['gate']['pass_gate'])
            counts.append(result['newton_updates']);states.append(result['state_interleaved'])
        self.assertLess(counts[1],counts[0])
        for a,b in zip(*states):self.assertAlmostEqual(a,b,delta=1e-12)
    def test_voltage_update_limit_rejects_nonpositive_values(self):
        self.input['diagnostic_voltage_update_limit_V']=0.
        self.write('input.json',self.input)
        run=self.run_deck();self.assertNotEqual(run.returncode,0)
    def test_density_bias_ceiling_preserves_high_bias_qf_path(self):
        self.input['diagnostic_density_update_iterations']=60
        self.write('input.json',self.input)
        self.deck['sweep'].update(bias_points_V=[0.,.1],density_update_maximum_bias_V=0.)
        self.deck['pause_after_attempts']=2
        run=self.run_deck();self.assertEqual(run.returncode,1,run.stderr)
        ledger=json.loads((self.root/'output/ledger.json').read_text())
        self.assertEqual(len(ledger['runs']),2)
        configs=[json.loads((Path(r['directory'])/'input.json').read_text()) for r in ledger['runs']]
        self.assertEqual(configs[0]['diagnostic_density_update_iterations'],60)
        self.assertEqual(configs[1]['diagnostic_density_update_iterations'],0)
        self.assertEqual(configs[1]['electrical_gate_solver'],self.input['electrical_gate_solver'])
        self.assertNotIn('solver_error',ledger['runs'][1]['gate']['reasons'])
        self.assertIn('heat_balance',ledger['runs'][1]['gate']['reasons'])
    def test_density_bias_ceiling_rejects_negative_value(self):
        self.deck['sweep']['density_update_maximum_bias_V']=-1.
        run=self.run_deck();self.assertNotEqual(run.returncode,0)
        self.assertFalse((self.root/'output').exists())
    def test_density_prediction_guard_requires_an_accepted_history(self):
        self.input['diagnostic_density_update_iterations']=60
        self.write('input.json',self.input)
        self.deck['sweep']['density_update_requires_prediction']=True
        run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
        ledger=json.loads((self.root/'output/ledger.json').read_text())
        self.assertFalse(ledger['runs'][0]['prediction']['used'])
        cfg=json.loads((Path(ledger['runs'][0]['directory'])/'input.json').read_text())
        self.assertEqual(cfg['diagnostic_density_update_iterations'],0)
    def test_density_prediction_guard_preserves_initialization_trajectory(self):
        self.deck['initialization']=dict(mode='neutral_300K',gate_voltage_V=.1,max_newton=10)
        self.deck['sweep']['density_update_requires_prediction']=True
        self.input['skip_equilibrium_poisson_transport']=True
        self.input['boundaries'][-4]['value']=.1
        trajectories=[]
        for count in (0,60):
            self.input['diagnostic_density_update_iterations']=count
            self.write('input.json',self.input)
            self.deck['output_directory']='initial_density_'+str(count)
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            ledger=json.loads((self.root/self.deck['output_directory']/'ledger.json').read_text())
            trajectory=[]
            self.assertGreater(len(ledger['initialization_runs']),3)
            for row in ledger['initialization_runs']:
                result_path=Path(row['result'])
                cfg=json.loads((result_path.parent/'input.json').read_text())
                self.assertEqual(cfg['diagnostic_density_update_iterations'],0)
                result=json.loads(result_path.read_text())
                self.assertTrue(row['gate']['pass_gate'])
                trajectory.append((row['gate_V'],row['solve_mode'],row['newton_updates'],
                                   result['state_interleaved']))
            trajectories.append(trajectory)
        self.assertEqual(trajectories[0],trajectories[1])
    def test_neutral_initialization_does_not_use_supplied_seed(self):
        self.deck['initialization']=dict(mode='neutral_300K',gate_voltage_V=.1,max_newton=10)
        self.input['skip_equilibrium_poisson_transport']=True
        self.input['boundaries'][-4]['value']=.1
        for key in ('state_interleaved','referenced_state_interleaved'):
            self.input[key]=[999.]*16
        self.write('input.json',self.input)
        run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
        ledger=json.loads((self.root/'output/ledger.json').read_text())
        self.assertGreater(len(ledger['initialization_runs']),3)
        result=json.loads(Path(ledger['accepted_result']).read_text())
        self.assertEqual(result['temperature_K'],[300.]*4)
        self.assertAlmostEqual(result['state_interleaved'][12],.1)
        self.assertTrue(all(abs(v)<1e-14 for v in result['state_interleaved'][1::4]))
        self.assertTrue(all(abs(v)<1e-14 for v in result['state_interleaved'][2::4]))


if __name__=='__main__':
    unittest.main()
