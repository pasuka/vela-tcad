"""Production electrothermal checkpoint, path and gate regression."""
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
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
    def test_ngmres_does_not_modify_a_converged_point(self):
        self.assertEqual(self.run_deck().returncode,0)
        baseline=json.loads((self.root/'output/step_0000/output.json').read_text())
        self.input['diagnostic_ngmres_recovery']=True;self.write('input.json',self.input)
        self.deck['output_directory']='ngmres'
        run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
        result=json.loads((self.root/'ngmres/step_0000/output.json').read_text())
        self.assertEqual(result['state_interleaved'],baseline['state_interleaved'])
        self.assertEqual(result['residual'],baseline['residual'])
        self.assertEqual(result['ngmres_recovery']['updates'],0)
        self.assertEqual(result['ngmres_recovery']['attempts'],[])
    def test_ngmres_recomputes_reference_history_and_preserves_failed_fallback(self):
        # Equal potential/QF shifts preserve carrier statistics while the large
        # Dirichlet target creates a deliberately cap-limited stagnating map.
        for boundary in self.input['boundaries']:
            if boundary['kind']!='temperature':boundary['value']=1e4
        self.input.update(performance_profiling=True,diagnostic_stagnation_window=3)
        self.deck['sweep']['max_newton']=20
        results=[]
        for enabled in (False,True):
            self.input['diagnostic_ngmres_recovery']=enabled;self.write('input.json',self.input)
            self.deck['output_directory']='stalled_'+str(enabled)
            run=self.run_deck();self.assertEqual(run.returncode,1)
            result=json.loads((self.root/self.deck['output_directory']/'step_0000/output.json').read_text())
            self.assertEqual(result['diagnostic_stop'],'diagnostic_stagnation_reject');results.append(result)
        for key in ('referenced_state_interleaved','electron_qf_reference_V','hole_qf_reference_V','residual','newton_updates'):
            self.assertEqual(results[0][key],results[1][key])
        recovery=results[1]['ngmres_recovery']
        self.assertGreater(recovery['reference_reassemblies'],0)
        self.assertEqual(recovery['updates'],0);self.assertEqual(len(recovery['attempts']),1)
        self.assertEqual(len(recovery['attempts'][0]['candidates']),4)
        self.assertTrue(all(not c['accepted'] for c in recovery['attempts'][0]['candidates']))
    def test_predictor_screening_selects_lower_residual_and_aligns_references(self):
        self.deck['sweep']['max_newton']=10
        self.input['performance_profiling']=True
        for boundary in self.input['boundaries']:
            if boundary['kind']=='psi':boundary['value']=.5
        results=[]
        for guarded in (False,True):
            self.deck['output_directory']='screen_'+str(guarded)
            if guarded:
                self.input['diagnostic_density_update_iterations']=60
                self.input['diagnostic_density_requires_primary_prediction']=True
                self.input['diagnostic_predictor_candidates']=[dict(label='exact',
                    referenced_state_interleaved=[.5,-.001,.002,300.]*4,
                    electron_qf_reference_V=[.001]*4,hole_qf_reference_V=[-.002]*4)]
            self.write('input.json',self.input)
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            ledger=json.loads((self.root/self.deck['output_directory']/'ledger.json').read_text())
            results.append(json.loads(Path(ledger['accepted_result']).read_text()))
        self.assertLess(results[1]['newton_updates'],results[0]['newton_updates'])
        self.assertEqual(results[1]['predictor_selection']['selected'],'exact')
        self.assertTrue(results[1]['predictor_selection']['density_update_disabled_after_fallback'])
        self.assertNotIn('density_coordinate_evaluations',results[1]['performance'])
        self.assertEqual(results[1]['performance']['residual_only_calls'],1)
        for a,b in zip(results[0]['state_interleaved'],results[1]['state_interleaved']):
            self.assertAlmostEqual(a,b,delta=1e-12)
    def test_predictor_screening_rejects_bad_residual_and_invalid_temperature(self):
        self.input['diagnostic_predictor_candidates']=[dict(label=label,
            referenced_state_interleaved=[psi,0.,0.,temp]*4,
            electron_qf_reference_V=[0.]*4,hole_qf_reference_V=[0.]*4)
            for label,psi,temp in (('worse',1.,300.),('invalid',0.,-1.))]
        self.write('input.json',self.input)
        run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
        result=json.loads((self.root/'output/step_0000/output.json').read_text())
        self.assertEqual(result['predictor_selection']['selected'],'primary')
        self.assertEqual(result['newton_updates'],0)
        self.assertIn('rejected_reason',result['predictor_selection']['candidates'][1])
    def test_tangent_prepares_from_source_bias_and_counts_extra_factorization(self):
        import copy
        self.input['boundaries']=[b for b in self.input['boundaries'] if not (b['node']==0 and b['kind']!='temperature')]
        self.input['boundaries'].append(dict(node=0,kind='neutral_contact',value=0.))
        self.input['performance_profiling']=True
        self.deck['sweep']['max_newton']=20
        self.write('input.json',self.input)
        run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
        zero=json.loads((self.root/'output/step_0000/output.json').read_text())
        keys=('state_interleaved','referenced_state_interleaved','electron_qf_reference_V','hole_qf_reference_V')
        source={key:copy.deepcopy(zero[key]) for key in keys}
        for key in ('state_interleaved','referenced_state_interleaved'):
            for k in range(3):source[key][k]-=.01
        source.update(source_bias_V=-.01,delta_bias_V=.01,moving_nodes=[0],potential_origin_V=0.)
        for key in keys:self.input[key]=copy.deepcopy(zero[key])
        self.input['state_interleaved'][0]+=.02;self.input['referenced_state_interleaved'][0]+=.02
        self.input['diagnostic_tangent_predictor']=source
        self.deck['output_directory']='tangent'
        self.write('input.json',self.input)
        run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
        result=json.loads((self.root/'tangent/step_0000/output.json').read_text())
        self.assertTrue(result['tangent_preparation']['prepared'])
        self.assertEqual(result['predictor_selection']['selected'],'tangent')
        self.assertEqual(result['performance']['factorizations'],result['newton_updates']+1)
        for a,b in zip(result['state_interleaved'],zero['state_interleaved']):self.assertAlmostEqual(a,b,delta=1e-10)
        source['moving_nodes']=[99];self.write('input.json',self.input)
        self.deck['output_directory']='tangent_fallback'
        run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
        result=json.loads((self.root/'tangent_fallback/step_0000/output.json').read_text())
        self.assertFalse(result['tangent_preparation']['prepared'])
        self.assertIn('fallback_reason',result['tangent_preparation'])
    def test_fixed_targets_attempts_exact_target_and_stops_on_failure(self):
        self.deck['sweep'].update(step_policy='fixed_targets',bias_points_V=[0.,.3,.4])
        run=self.run_deck();self.assertEqual(run.returncode,1,run.stderr)
        ledger=json.loads((self.root/'output/ledger.json').read_text())
        self.assertEqual([r['bias_V'] for r in ledger['runs']],[0.,.3])
        self.assertEqual(ledger['status'],'failed')
        self.assertEqual(ledger['accepted_bias_V'],0.)
    def test_predictor_study_checks_provenance_and_records_all_runs(self):
        self.assertEqual(self.run_deck().returncode,0)
        evidence=self.root/'evidence';(evidence/'cases/data').mkdir(parents=True)
        (evidence/'cases_r7/data').mkdir(parents=True)
        (evidence/'cases/data/mesh.json').write_bytes((self.root/'mesh.json').read_bytes())
        (evidence/'cases/data/ialmob_geometry.json').write_text('{}')
        for gate in (4,8):
            (evidence/f'cases_r7/data/input_vg{gate}.json').write_text(json.dumps(self.input))
            (evidence/f'cases_r7/vela_vg{gate}.json').write_text(json.dumps(self.deck))
            directory=evidence/f'results/r7_full_vg{gate}/step_0000';directory.mkdir(parents=True)
            (directory/'output.json').write_bytes((self.root/'output/step_0000/output.json').read_bytes())
            (directory.parent/'ledger.json').write_text(json.dumps(dict(runs=[dict(bias_V=0.,gate=dict(pass_gate=True),
                directory=f'/original/results/r7_full_vg{gate}/step_0000')])))
        sources=[p for p in evidence.rglob('*.json') if 'results' not in p.parts]
        profile=dict(files_sha256={p.relative_to(evidence).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                     dependencies=[dict(recorded_path='/original/cases/data/mesh.json')])
        self.write('profile.json',profile)
        script=Path(__file__).resolve().parents[2]/'scripts/run_templates_ldmos_predictor_study.py'
        argv=[sys.executable,str(script),'--evidence-root',str(evidence),'--profile',str(self.root/'profile.json'),
              '--runner',self.runner,'--output',str(self.root/'study'),'--maximum-bias','0']
        run=subprocess.run(argv,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        self.assertEqual(run.returncode,0,run.stderr)
        report=json.loads((self.root/'study/summary.json').read_text())
        self.assertEqual(len(report['runs']),6)
        self.assertTrue(all(r['status']=='complete' and r['attempts']==1 and r['newton_updates']==0 for r in report['runs']))
        (evidence/'cases/data/mesh.json').write_text('{}')
        argv[argv.index('--output')+1]=str(self.root/'changed_study')
        run=subprocess.run(argv,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        self.assertNotEqual(run.returncode,0)
        self.assertFalse((self.root/'changed_study').exists())
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


class PredictorComparisonTest(unittest.TestCase):
    def test_voltage_serialization_is_explicit_and_never_interpolates(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);rows=[]
            for variant,bias in (('baseline',.3),('actual_step',.30000000000000004)):
                directory=root/(variant+'_vg4')/'results';directory.mkdir(parents=True)
                point=directory/'point';point.mkdir()
                state=dict(potential_origin_V=0.,referenced_state_interleaved=[bias,0.,0.,300.],
                           electron_qf_reference_V=[0.],hole_qf_reference_V=[0.])
                (point/'output.json').write_text(json.dumps(state))
                (point/'input.json').write_text(json.dumps(dict(boundaries=[dict(node=0,kind='neutral_contact',value=bias)])))
                (directory/'ledger.json').write_text(json.dumps(dict(runs=[dict(bias_V=bias,directory=str(point),gate=dict(pass_gate=True))])))
                rows.append(dict(gate_V=4,variant=variant,status='complete',targets_V=[.3],failed_attempts=0,screening_seconds=0.))
            (root/'summary.json').write_text(json.dumps(dict(runs=rows)))
            script=Path(__file__).resolve().parents[2]/'scripts/analyze_templates_ldmos_predictor_study.py'
            argv=[sys.executable,str(script),'--matrix',str(root),'--output',str(root/'comparison.json')]
            run=subprocess.run(argv,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
            self.assertNotEqual(run.returncode,0)
            self.assertIn('missing requested',run.stderr)
            run=subprocess.run(argv+['--bias-digits','15'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
            self.assertEqual(run.returncode,0,run.stderr)
            result=json.loads((root/'comparison.json').read_text())['comparisons'][0]
            self.assertEqual(result['common_accepted_states'],1)
            self.assertGreater(result['state_comparison'][0]['max_physical_state_delta_V_V_V_K'][0],0.)
            self.assertFalse(result['all_compared_states_exact'])
    def test_comparison_retains_sub_ulp_qf_changes(self):
        sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
        from analyze_templates_ldmos_predictor_study import state_delta
        a=dict(potential_origin_V=0.,referenced_state_interleaved=[0.,1e-20,0.,300.],
               electron_qf_reference_V=[1e12],hole_qf_reference_V=[0.])
        b=dict(a,referenced_state_interleaved=[0.,0.,0.,300.])
        self.assertEqual(state_delta(a,b),[0.,1e-20,0.,0.])
        b['referenced_state_interleaved'][0]=float('nan')
        with self.assertRaisesRegex(ValueError,'Nonfinite'):state_delta(a,b)


if __name__=='__main__':
    unittest.main()
