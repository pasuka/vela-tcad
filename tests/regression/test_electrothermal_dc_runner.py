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

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import electrothermal_state
import state_archive


def without_trial_timing(value):
    if isinstance(value,dict):
        return {k:without_trial_timing(v) for k,v in value.items() if k!='trial_seconds'}
    if isinstance(value,list):
        return [without_trial_timing(v) for v in value]
    return value


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
        self.mesh_identity=state_archive.mesh_identity(mesh,1.)
        self.state_serial=0
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
        if 'state_interleaved' in data or 'referenced_state_interleaved' in data:
            self.state_serial+=1
            mesh=json.loads((self.root/Path(data['mesh_file'])).read_text())
            identity=state_archive.mesh_identity(mesh,data.get('coordinate_to_metres',1.))
            data=electrothermal_state.pack(self.root/f'input_{self.state_serial}.json',data,
                dict(mode='electrothermal',mesh_sha256=identity,potential_origin_V=data.get('potential_origin_V',0.)))
        (self.root/name).write_text(json.dumps(data),encoding='utf-8')
    def read(self,path):
        return electrothermal_state.read_bound_record(path)
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
    def test_screening_default_and_explicit_legacy_reach_point_service(self):
        self.input['performance_profiling']=True
        for method in (None,'legacy','halley'):
            with self.subTest(method=method):
                self.input.pop('diagnostic_ialmob_screening_method',None)
                if method is not None:
                    self.input['diagnostic_ialmob_screening_method']=method
                self.write('input.json',self.input)
                self.deck['output_directory']='screening_'+str(method)
                run=self.run_deck()
                self.assertEqual(run.returncode,0,run.stderr)
                result=self.read(self.root/self.deck['output_directory']/'step_0000/output.json')
                self.assertEqual(result['performance']['ialmob_screening_method'],method or 'halley')
    def test_completed_resume_rejects_corrupted_referenced_state(self):
        self.assertEqual(self.run_deck().returncode,0)
        path=self.root/'output/step_0000/output.json'
        record=json.loads(path.read_text())
        state=path.parent/record['state_archive']['file']
        with state.open('ab') as stream: stream.write(b'corruption')
        self.deck['resume']=True
        run=self.run_deck()
        self.assertNotEqual(run.returncode,0)
        self.assertIn('digest mismatch',run.stderr)
    def test_missing_temperature_is_not_replaced_by_ambient(self):
        import h5py
        path=self.root/'input.json';record=json.loads(path.read_text())
        state=path.parent/record['state_archive']['file']
        with h5py.File(state,'r+') as file: del file['fields/temperature_K']
        record['state_archive']['sha256']=hashlib.sha256(state.read_bytes()).hexdigest()
        path.write_text(json.dumps(record))
        run=self.run_deck()
        self.assertNotEqual(run.returncode,0)
        self.assertIn('temperature/mode mismatch',run.stderr)
        self.assertFalse((self.root/'output').exists())
    def test_resume_rejects_modified_geometry(self):
        self.assertEqual(self.run_deck().returncode,0)
        mesh=self.read((self.root/'mesh.json'));mesh['nodes'][0]['x']=.1;self.write('mesh.json',mesh)
        self.deck['resume']=True
        run=self.run_deck();self.assertNotEqual(run.returncode,0);self.assertIn('mesh identity mismatch',run.stderr)
    def test_resume_rejects_changed_external_input_bytes(self):
        # A declared source is fingerprinted even when this minimal fixture's
        # physics is specified inline. The path and main deck remain unchanged.
        (self.root/'material_source.json').write_text('{}')
        self.input['materials_file']=str(self.root/'material_source.json')
        self.write('input.json',self.input)
        self.assertEqual(self.run_deck().returncode,0)
        (self.root/'material_source.json').write_text('{"revision":2}')
        self.deck['resume']=True
        run=self.run_deck()
        self.assertNotEqual(run.returncode,0)
        self.assertIn('differs from checkpoint',run.stderr)
    def test_pause_then_failed_heat_gate_keeps_accepted_zero(self):
        self.deck['sweep']['bias_points_V']=[0.,.1];self.deck['pause_after_attempts']=1
        run=self.run_deck();self.assertEqual(run.returncode,1,run.stderr)
        ledger=self.read((self.root/'output/ledger.json'))
        self.assertEqual(ledger['status'],'stopped_at_checkpoint')
        accepted=Path(ledger['accepted_result']).read_bytes()
        self.deck.update(resume=True,pause_after_attempts=0)
        run=self.run_deck();self.assertEqual(run.returncode,1,run.stderr)
        ledger=self.read((self.root/'output/ledger.json'))
        self.assertEqual(ledger['status'],'failed');self.assertEqual(ledger['accepted_bias_V'],0.)
        self.assertEqual(Path(ledger['accepted_result']).read_bytes(),accepted)
        self.assertIn('heat_balance',ledger['runs'][-1]['gate']['reasons'])
    def test_disabling_row_gate_is_rejected_before_output(self):
        self.input['electrical_gate_solver']['carrier_row_convergence']['mode']='off';self.write('input.json',self.input)
        run=self.run_deck();self.assertNotEqual(run.returncode,0)
        self.assertFalse((self.root/'output').exists())
    def test_static_preparation_is_fresh_on_resume_and_shared_after_failure(self):
        self.deck['reuse_static_preparation']=True
        self.input['performance_profiling']=True
        self.write('input.json',self.input)
        self.test_pause_then_failed_heat_gate_keeps_accepted_zero()
        ledger=self.read((self.root/'output/ledger.json'))
        self.assertGreater(len(ledger['runs']),2)
        hits=[self.read((Path(r['directory'])/'output.json'))
              ['performance']['static_preparation_reused'] for r in ledger['runs']]
        self.assertEqual(hits[:2],[False,False])
        self.assertTrue(all(hits[2:]))
        linear_hits=[self.read(Path(r['directory'])/'output.json')['performance']['linear_object_reused']
                     for r in ledger['runs']]
        self.assertFalse(any(linear_hits))  # Fresh resume, then externally rejected heat gates.
    def test_structure_cache_is_fresh_on_resume_and_safe_after_failure(self):
        self.deck['reuse_jacobian_structure']=True
        self.input['performance_profiling']=True
        self.write('input.json',self.input)
        self.test_pause_then_failed_heat_gate_keeps_accepted_zero()
        ledger=self.read((self.root/'output/ledger.json'))
        self.assertGreater(len(ledger['runs']),2)
        hits=[self.read((Path(r['directory'])/'output.json'))
              ['performance']['static_preparation_reused'] for r in ledger['runs']]
        self.assertEqual(hits[:2],[False,False])
        self.assertTrue(all(hits[2:]))
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
            ledger=self.read((self.root/self.deck['output_directory']/'ledger.json'))
            result=self.read(Path(ledger['accepted_result']))
            self.assertTrue(ledger['runs'][0]['gate']['pass_gate'])
            counts.append(result['newton_updates']);states.append(result['state_interleaved'])
        self.assertLess(counts[1],counts[0])
        for a,b in zip(*states):self.assertAlmostEqual(a,b,delta=1e-12)
    def test_iteration_trace_preserves_numerical_trajectory(self):
        self.deck['sweep']['max_newton']=10
        for boundary in self.input['boundaries']:
            if boundary['kind']=='psi':boundary['value']=.5
        results=[]
        for trace in (False,True):
            self.input['diagnostic_iteration_trace']=trace
            self.input['performance_profiling']=True
            self.write('input.json',self.input)
            self.deck['output_directory']='trace_'+str(trace)
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            results.append(self.read((self.root/self.deck['output_directory']/'step_0000/output.json')))
        baseline,traced=results
        for key in ('state_interleaved','referenced_state_interleaved','residual','diagnostic_stop','newton_updates'):
            self.assertEqual(baseline[key],traced[key],key)
        self.assertGreater(traced['iteration_trace_seconds'],0.)
        self.assertGreater(len(traced['history']),1)
        for a,b in zip(baseline['history'],traced['history']):
            trace=b.pop('iteration_trace');self.assertEqual(a,b)
            self.assertEqual(len(trace['direction_maxima']),4)
            self.assertEqual(trace['direction_maxima'][0]['temperature_K'],300.)
            self.assertEqual(len(trace['line_search_candidates']),b['line_search_trials'])
            self.assertEqual(sum(x['accepted'] for x in trace['line_search_candidates']),int(b['accepted']))
            self.assertTrue(all(x['trial_seconds']>=0 for x in trace['line_search_candidates']))
            for stage in ('before_gates','after_gates'):
                self.assertEqual(len(trace[stage]['blocks']),3)
                self.assertIn('eps_row',trace[stage]['row'])

    def test_local_qf_limiter_does_not_reduce_the_potential_cap(self):
        self.deck['sweep']['max_newton']=15
        self.input['performance_profiling']=True
        for boundary in self.input['boundaries']:
            if boundary['kind']=='psi':boundary['value']=.5
            if boundary['kind']=='fn':boundary['value']=1.
        results=[]
        for enabled in (False,True):
            self.input['diagnostic_local_qf_limiter']=enabled
            self.write('input.json',self.input);self.deck['output_directory']='local_qf_'+str(enabled)
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            result=self.read((self.root/self.deck['output_directory']/'step_0000/output.json'))
            self.assertTrue(result['carrier_row_gate']['satisfied'])
            results.append(result)
        self.assertGreater(results[1]['history'][0]['initial_alpha'],results[0]['history'][0]['initial_alpha'])
        self.assertGreater(results[1]['history'][0]['local_qf_limiter']['clipped_electrons'],0)
        for a,b in zip(results[0]['state_interleaved'],results[1]['state_interleaved']):self.assertAlmostEqual(a,b,delta=1e-12)

    def test_local_qf_limiter_is_disabled_during_initialization(self):
        self.deck['initialization']=dict(mode='neutral_300K',gate_voltage_V=.1,max_newton=10)
        self.input['boundaries'][-4]['value']=.1
        self.input['diagnostic_local_qf_limiter']=True;self.write('input.json',self.input)
        run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
        ledger=self.read((self.root/'output/ledger.json'))
        for row in ledger['initialization_runs']:
            cfg=self.read((Path(row['result']).parent/'input.json'))
            self.assertFalse(cfg['diagnostic_local_qf_limiter'])
        cfg=self.read((Path(ledger['runs'][0]['directory'])/'input.json'))
        self.assertTrue(cfg['diagnostic_local_qf_limiter'])

    def test_projected_natural_updates_preserve_dirichlet_constraints(self):
        self.deck['sweep']['max_newton']=10
        for boundary in self.input['boundaries']:
            if boundary['kind']=='psi':boundary['value']=.5
        results=[]
        for mode in ('off','v1','v2'):
            self.input['diagnostic_density_projection']=mode
            self.input['diagnostic_natural_damping']=True
            self.input['performance_profiling']=True
            self.write('input.json',self.input)
            self.deck['output_directory']='natural_'+mode
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            result=self.read((self.root/self.deck['output_directory']/'step_0000/output.json'))
            self.assertGreater(result['natural_damping']['corrector_solves'],0)
            self.assertTrue(result['carrier_row_gate']['satisfied'])
            self.assertTrue(all(b['satisfied'] for b in result['electrical_block_gates']))
            results.append(result['state_interleaved'])
        self.assertEqual(results[0],results[1]);self.assertEqual(results[0],results[2])

    def test_adaptive_jacobian_closes_a_nonlinear_free_poisson_row(self):
        # Millimetre geometry keeps charge-subtraction roundoff below the absolute residual ceiling.
        self.input['coordinate_to_metres']=1e-3
        self.input['silicon_area_m2']=[a*1e-6 for a in self.input['silicon_area_m2']]
        self.deck['sweep']['max_newton']=30
        self.input['boundaries']=[b for b in self.input['boundaries'] if not (b['node']==1 and b['kind']=='psi')]
        for key in ('state_interleaved','referenced_state_interleaved'):self.input[key][4]=.1
        self.input['performance_profiling']=True
        results=[]
        for adaptive in (False,True):
            self.input['diagnostic_adaptive_jacobian']=adaptive;self.write('input.json',self.input)
            self.deck['output_directory']='adaptive_'+str(adaptive)
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            result=self.read((self.root/self.deck['output_directory']/'step_0000/output.json'))
            self.assertTrue(all(b['satisfied'] for b in result['electrical_block_gates']))
            results.append(result)
        self.assertGreater(results[1]['adaptive_jacobian']['lagged_solves'],0)
        for a,b in zip(results[0]['state_interleaved'],results[1]['state_interleaved']):self.assertAlmostEqual(a,b,delta=1e-8)
        self.assertLess(results[1]['performance']['factorizations'],results[1]['newton_updates'])

    def test_row_audit_decimal_sg_keeps_sub_ulp_difference(self):
        from decimal import Decimal, localcontext
        sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
        from analyze_templates_ldmos_row_roundoff import sg,D,contact
        row=dict(q=1.602176634e-19,kb=1.380649e-23)
        edge=dict(sa=dict(T=300.,rp=0.,fp=1.),sb=dict(T=300.,rp=1.,fp=1e-20),
                  pa=dict(p=1e10,Nv=1e25,ep=-35.),pb=dict(p=1e10,Nv=1e25,ep=-35.),
                  g=1.,logF=0.,mu=.03,weight=.7)
        with localcontext() as ctx:
            ctx.prec=100
            exact=sg(edge,row);rounded=sg(edge,row,rounded_qf=True)
            target=-D(row['q'])*D(edge['mu'])*D(edge['weight'])*D(1e10)*D(1e-20)
            self.assertLess(abs((exact-target)/target),Decimal('1e-65'))
            self.assertEqual(rounded,0)
            edge['sb']['fp']=0.
            self.assertEqual(sg(edge,row),0)
            boundary=dict(state=dict(rp=0.,fp=1e-20,psi=.5),contact=dict(bias=0.,neutral_potential=.5,
                vt=1.,coefficient=1.,Nv=1.,df=1.,ddf=0.))
            self.assertEqual(contact(boundary),D(1e-20))
            self.assertEqual(contact(boundary,D(1e-17)),D(1e-20)+D(1e-17))

    def test_hole_row_audit_is_read_only_and_rejects_updates(self):
        probe=Path(self.runner).parent/('electrothermal_probe.exe' if os.name=='nt' else 'electrothermal_probe')
        cfg=dict(self.input,mesh_file=str(self.root/'mesh.json'),solve_mode='coupled',initialization='provided_state',diagnostic_newton_max_iterations=0)
        cfg['boundaries']=[b for b in cfg['boundaries'] if not (b['node']==1 and b['kind']=='fp')]
        def run(value,label):
            self.write(label+'.json',value)
            proc=subprocess.run([str(probe),str(self.root/(label+'.json')),str(self.root/(label+'_out.json'))],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
            return proc,self.read((self.root/(label+'_out.json'))) if proc.returncode==0 else None
        p,base=run(cfg,'base');self.assertEqual(p.returncode,0,p.stderr)
        cfg['diagnostic_hole_row_audit_nodes']=[1]
        p,data=run(cfg,'audit');self.assertEqual(p.returncode,0,p.stderr)
        self.assertEqual(data['residual'],base['residual']);self.assertEqual(data['state_interleaved'],base['state_interleaved'])
        self.assertEqual(data['hole_row_audit'][0]['residual'],base['residual'][6])
        cfg['boundaries']=[b for b in cfg['boundaries'] if not (b['node']==1 and b['kind'] in ('psi','fn','fp'))]
        cfg['boundaries'].append(dict(node=1,kind='neutral_contact',value=0.,hole_recombination_velocity_m_per_s=1.,boundary_length_m=1e-6))
        for key in ('state_interleaved','referenced_state_interleaved'):cfg[key][4]=.1
        p,data=run(cfg,'raw_contact');self.assertEqual(p.returncode,0,p.stderr)
        self.assertEqual(data['state_interleaved'][4],.1)
        self.assertNotEqual(data['hole_row_audit'][0]['contact']['neutral_potential'],.1)
        cfg['diagnostic_newton_max_iterations']=1
        p,_=run(cfg,'invalid');self.assertNotEqual(p.returncode,0)
        self.assertIn('zero-update',p.stderr)

    def test_poisson_preparation_holds_qf_temperature_and_restores_coupled_rows(self):
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
        from run_templates_ldmos_poisson_initialization import restored_config
        from analyze_templates_ldmos_predictor_study import state_delta
        self.input.update(coordinate_to_metres=1e-3, potential_origin_V=0.,
                          diagnostic_newton_max_iterations=30, performance_profiling=True)
        self.input['silicon_area_m2']=[a*1e-6 for a in self.input['silicon_area_m2']]
        self.input['boundaries']=[b for b in self.input['boundaries'] if not (b['node']==1 and b['kind']=='psi')]
        for key in ('state_interleaved','referenced_state_interleaved'):
            self.input[key][4:8]=[.1,.02,.01,310.]
        original=json.loads(json.dumps(self.input))
        self.input['solve_mode']='poisson'
        probe=Path(self.runner).parent/'electrothermal_probe'
        if os.name=='nt':probe=probe.with_suffix('.exe')
        def point(cfg,name):
            cfg=dict(cfg,mesh_file=str(self.root/'mesh.json'))
            self.write(name+'_input.json',cfg)
            run=subprocess.run([str(probe),str(self.root/(name+'_input.json')),str(self.root/(name+'.json'))],
                               stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
            self.assertEqual(run.returncode,0,run.stderr)
            return self.read((self.root/(name+'.json')))
        prepared=point(self.input,'prepared')
        self.assertEqual(prepared['diagnostic_stop'],'diagnostic_scaled_residual')
        self.assertTrue(prepared['electrical_block_gates'][0]['satisfied'])
        delta=state_delta(original,prepared)
        self.assertGreater(delta[0],.01)
        for change in delta[1:]:self.assertLessEqual(change,1e-14)
        restored=restored_config(original,prepared)
        self.assertEqual(restored['boundaries'],original['boundaries'])
        restored['diagnostic_newton_max_iterations']=0
        audit=point(restored,'audit')
        self.assertEqual(audit['newton_updates'],0)
        self.assertEqual(state_delta(prepared,audit),[0.]*4)
        self.assertTrue(audit['electrical_block_gates'][0]['satisfied'])
        # Original carrier and thermal constraints must reappear after restoration.
        self.assertAlmostEqual(audit['residual'][5],.02)
        self.assertAlmostEqual(audit['residual'][6],.01)
        self.assertAlmostEqual(audit['residual'][7],10.)
        self.assertFalse(all(b['satisfied'] for b in audit['electrical_block_gates']))

    def test_poisson_preparation_gate_rejects_incomplete_or_changed_state(self):
        import copy
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
        from run_templates_ldmos_poisson_initialization import preparation_gate, restored_config, total_cost
        original=dict(self.input,potential_origin_V=0.)
        prepared=dict(copy.deepcopy(original),diagnostic_stop='diagnostic_scaled_residual',newton_updates=3,
                      electrical_block_gates=[dict(satisfied=True,weighted_l2=1e-9,limit=5e-8)])
        audit=dict(copy.deepcopy(prepared),newton_updates=0)
        self.assertTrue(preparation_gate(original,prepared,audit)['pass_gate'])
        for field,value in (('weighted_l2',1.),('satisfied',False)):
            broken=copy.deepcopy(audit);broken['electrical_block_gates'][0][field]=value
            self.assertFalse(preparation_gate(original,prepared,broken)['pass_gate'])
        broken=copy.deepcopy(prepared);broken['diagnostic_stop']='line_search_failed'
        self.assertFalse(preparation_gate(original,broken,audit)['pass_gate'])
        broken=copy.deepcopy(prepared);broken['electron_qf_reference_V'][1]=1e-6
        self.assertFalse(preparation_gate(original,broken,audit)['pass_gate'])
        broken=copy.deepcopy(prepared);broken['potential_origin_V']=1.
        with self.assertRaises(ValueError):restored_config(original,broken)
        original['boundaries'].append(dict(node=1,kind='neutral_contact',value=0.,
            hole_recombination_velocity_m_per_s=123.,boundary_length_m=1e-6))
        restored=restored_config(original,prepared)
        self.assertEqual(restored['boundaries'],original['boundaries'])
        restored['boundaries'][-1]['hole_recombination_velocity_m_per_s']=0.
        self.assertEqual(original['boundaries'][-1]['hole_recombination_velocity_m_per_s'],123.)
        keys=('newton_updates','attempts','trials','assemblies','factorizations','wall_seconds')
        self.assertEqual(total_cost(dict.fromkeys(keys,3),dict.fromkeys(keys,0),dict.fromkeys(keys,5)),dict.fromkeys(keys,8))

    def test_near_steady_rebase_preserves_sub_ulp_qf_without_copying_solution(self):
        import copy
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
        from run_templates_ldmos_near_steady_isolation import rebase_state, zero_small_qf_references
        from analyze_templates_ldmos_predictor_study import state_delta
        source=dict(copy.deepcopy(self.input),potential_origin_V=0.)
        source['electron_qf_reference_V']=[4.]*4
        source['hole_qf_reference_V']=[-2.]*4
        for i in range(4):
            source['referenced_state_interleaved'][4*i+1]=1e-17
            source['referenced_state_interleaved'][4*i+2]=-2e-17
        reference=copy.deepcopy(source)
        reference['electron_qf_reference_V']=[4.+2**-12]*4
        reference['hole_qf_reference_V']=[-2.-2**-12]*4
        reference['referenced_state_interleaved']=[9.,9.,9.,999.]*4
        result=rebase_state(source,reference)
        delta=state_delta(source,result)
        self.assertEqual(delta[0],0.);self.assertEqual(delta[3],0.)
        self.assertLess(delta[1],1e-19);self.assertLess(delta[2],1e-19)
        self.assertNotEqual(result['referenced_state_interleaved'][1],-2**-12)
        self.assertEqual(source['electron_qf_reference_V'],[4.]*4)
        source['electron_qf_reference_V'][0]=2**-12
        source['referenced_state_interleaved'][1]=-2**-12+1e-17
        local=zero_small_qf_references(source)
        self.assertEqual(local['electron_qf_reference_V'][0],0.)
        self.assertEqual(local['electron_qf_reference_V'][1:],source['electron_qf_reference_V'][1:])
        self.assertEqual(state_delta(source,local),[0.]*4)

    def test_pseudo_storage_preserves_algebraic_constraints(self):
        self.deck['sweep']['max_newton']=10
        for boundary in self.input['boundaries']:
            if boundary['kind']=='psi':boundary['value']=.5
        self.input['performance_profiling']=True
        results=[]
        for enabled in (False,True):
            self.input['diagnostic_pseudo_transient']=enabled;self.write('input.json',self.input)
            self.deck['output_directory']='pseudo_constraints_'+str(enabled)
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            results.append(self.read((self.root/self.deck['output_directory']/'step_0000/output.json')))
        for key in ('state_interleaved','residual','newton_updates','electrical_block_gates'):
            self.assertEqual(results[0][key],results[1][key])
        steps=results[1]['pseudo_transient']['steps'];self.assertGreater(len(steps),0)
        self.assertTrue(all(s['mass_nonzeros']==0 for s in steps))
        for mode in ('defect_ser','defect_model'):
            self.input['diagnostic_pseudo_acceptance']=mode;self.write('input.json',self.input)
            self.deck['output_directory']='algebraic_'+mode
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            result=self.read((self.root/self.deck['output_directory']/'step_0000/output.json'))
            for key in ('state_interleaved','residual','newton_updates','electrical_block_gates'):
                self.assertEqual(results[0][key],result[key])

    def test_pseudo_defect_requires_mass_and_preserves_original_terminal_gates(self):
        self.input['diagnostic_pseudo_acceptance']='defect_model';self.write('input.json',self.input)
        run=self.run_deck();self.assertNotEqual(run.returncode,0)
        self.assertIn('requires pseudo transient',(self.root/'output/step_0000/run.log').read_text())
        # A converged state needs no pseudo step; the original gates and
        # steady residual must still be emitted unchanged.
        self.input['diagnostic_pseudo_transient']=True;self.write('input.json',self.input)
        self.deck['output_directory']='converged_defect'
        run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
        result=self.read((self.root/'converged_defect/step_0000/output.json'))
        self.assertEqual(result['newton_updates'],0)
        self.assertEqual(result['pseudo_transient']['defect_evaluations'],0)
        self.assertEqual(result['carrier_row_gate']['eps_row'],1e-8)
        self.assertEqual([b['limit'] for b in result['electrical_block_gates']],[1e-8,1e-11,3e-10])

    def test_pseudo_storage_has_free_carrier_rows_and_rejects_mixed_experiments(self):
        self.input['boundaries']=[b for b in self.input['boundaries'] if not (b['node']==1 and b['kind'] in ('fn','fp'))]
        self.input['coordinate_to_metres']=1e-6
        self.input['silicon_area_m2']=[a*1e-12 for a in self.input['silicon_area_m2']]
        self.input['recombination_area_m2']=[a*2 for a in self.input['silicon_area_m2']]
        for edge in self.input['edge_geometry']:edge['transport_weight']=1.
        for key in ('state_interleaved','referenced_state_interleaved'):
            self.input[key][5]=.01;self.input[key][6]=-.01
        self.input.update(diagnostic_pseudo_transient=True,performance_profiling=True)
        self.deck['sweep']['max_newton']=1
        self.deck['sweep']['growth_newton']=1
        self.write('input.json',self.input)
        run=self.run_deck()
        self.assertTrue((self.root/'output/step_0000/output.json').exists(),run.stdout+run.stderr)
        result=self.read((self.root/'output/step_0000/output.json'))
        first=result['pseudo_transient']['steps'][0]
        self.assertGreater(first['mass_nonzeros'],0);self.assertLessEqual(first['mass_nonzeros'],8)
        self.assertGreater(first['tau_s'],0.)
        self.assertEqual(result['state_interleaved'][3::4],[300.]*4)
        for mode in ('defect_ser','defect_model'):
            self.input['diagnostic_pseudo_acceptance']=mode;self.write('input.json',self.input)
            self.deck['output_directory']=mode
            self.run_deck()
            trial=self.read((self.root/mode/'step_0000/output.json'))
            step=trial['pseudo_transient']['steps'][0]
            self.assertTrue(step['accepted']);self.assertTrue(step['defect_used'])
            self.assertLess(step['accepted_defect_norm'],step['fixed_residual_before'])
            self.assertEqual(trial['state_interleaved'][3::4],[300.]*4)
            self.assertEqual(trial['state_interleaved'][0::4],[0.]*4)
            self.assertGreater(trial['pseudo_transient']['defect_evaluations'],0)
        del self.input['diagnostic_pseudo_acceptance']
        self.input['diagnostic_pseudo_direction_audit']=True;self.write('input.json',self.input)
        self.deck['output_directory']='audit'
        self.run_deck()
        audited=self.read((self.root/'audit/step_0000/output.json'))
        for key in ('state_interleaved','residual','newton_updates'):
            self.assertEqual(result[key],audited[key])
        directions=audited['pseudo_direction_audit']['directions'];self.assertEqual(len(directions),8)
        self.assertTrue(all(d['steady_residual_exact'] for d in directions))
        self.assertTrue(all(d['relative_linear_residual']<1e-10 for d in directions))
        del self.input['diagnostic_pseudo_direction_audit']
        self.input['diagnostic_adaptive_jacobian']=True;self.write('input.json',self.input)
        self.deck['output_directory']='mixed'
        run=self.run_deck();self.assertNotEqual(run.returncode,0)
        self.assertIn('must be isolated',(self.root/'mixed/step_0000/run.log').read_text())

    def test_near_steady_switch_guards_and_retirement(self):
        self.check_near_steady_representation(standalone=False)

    def test_near_steady_rebase_isolated_r7(self):
        self.check_near_steady_representation(standalone=True)

    def test_contact_consistency_preserves_inactive_and_rejected_newton(self):
        self.input['boundaries']=[b for b in self.input['boundaries'] if not
            (b['node']==1 and b['kind'] in ('psi','fn','fp'))]
        self.input['boundaries'].append(dict(node=1,kind='neutral_contact',value=0.,
            hole_recombination_velocity_m_per_s=1.,boundary_length_m=1e-6))
        self.input['electrical_gate_solver']['carrier_row_convergence'].update(min_flux_scale=1e-100,scale_floor=1e-30)
        self.input.update(performance_profiling=True,diagnostic_iteration_trace=True)
        self.deck['sweep'].update(max_newton=8,growth_newton=1)
        for amplitude in (0.,5e-14):
            for key in ('state_interleaved','referenced_state_interleaved'):self.input[key][6]=amplitude
            pair=[]
            for enabled in (False,True):
                self.input['diagnostic_near_steady_contact_consistency']=enabled
                self.write('input.json',self.input)
                self.deck['output_directory']='contact_'+str(amplitude)+'_'+str(enabled)
                self.run_deck()
                pair.append(self.read((self.root/self.deck['output_directory']/'step_0000/output.json')))
            for key in ('history','state_interleaved','referenced_state_interleaved','residual','carrier_row_gate','newton_updates'):
                self.assertEqual(without_trial_timing(pair[0][key]),without_trial_timing(pair[1][key]),key)
            event=pair[1]['near_steady_contact_consistency']
            if amplitude==0.:
                self.assertFalse(event['triggered'])
            else:
                self.assertTrue(event['triggered']);self.assertEqual(len(event['events']),1)
                self.assertFalse(event['events'][0]['accepted'])
                self.assertEqual(event['events'][0]['reason'],'already_consistent')
        self.input['diagnostic_density_update_iterations']=1
        self.write('input.json',self.input);self.deck['output_directory']='invalid_contact'
        self.assertNotEqual(self.run_deck().returncode,0)
        self.assertIn('Contact consistency requires',(self.root/'invalid_contact/step_0000/run.log').read_text())

    def test_contact_consistency_rolls_back_a_reassembled_trial(self):
        # A temperature update makes the finite-contact potential inconsistent,
        # while another free hole row cannot be repaired by a contact projection.
        # Its rejection must preserve the complete ordinary Newton trajectory.
        self.input['boundaries']=[b for b in self.input['boundaries'] if not
            (b['node']==1 and b['kind'] in ('psi','fn','fp')) and not (b['node']==2 and b['kind']=='fp')]
        self.input['boundaries'].append(dict(node=1,kind='neutral_contact',value=0.,
            hole_recombination_velocity_m_per_s=1.,boundary_length_m=1e-6))
        self.input['electrical_gate_solver']['carrier_row_convergence'].update(min_flux_scale=1e-100,scale_floor=1e-30)
        self.input.update(performance_profiling=True,diagnostic_iteration_trace=True)
        for key in ('state_interleaved','referenced_state_interleaved'):
            self.input[key][7]=299.99;self.input[key][10]=1e-9
        self.deck['sweep'].update(max_newton=10,growth_newton=1)
        pair=[]
        for enabled in (False,True):
            self.input['diagnostic_near_steady_contact_consistency']=enabled;self.write('input.json',self.input)
            self.deck['output_directory']='rollback_contact_'+str(enabled);self.run_deck()
            pair.append(self.read((self.root/self.deck['output_directory']/'step_0000/output.json')))
        event=pair[1]['near_steady_contact_consistency']['events']
        self.assertEqual(len(event),1);self.assertTrue(event[0]['reassembled']);self.assertFalse(event[0]['accepted'])
        self.assertGreater(event[0]['max_potential_change_V'],0.)
        for key in ('history','state_interleaved','referenced_state_interleaved','electron_qf_reference_V',
                    'hole_qf_reference_V','residual','carrier_row_gate','electrical_block_gates','diagnostic_stop'):
            self.assertEqual(without_trial_timing(pair[0][key]),without_trial_timing(pair[1][key]),key)
        self.assertEqual(pair[1]['performance']['assembly_calls'],pair[0]['performance']['assembly_calls']+1)
        self.assertEqual(pair[1]['performance']['factorizations'],pair[0]['performance']['factorizations'])

    def check_near_steady_representation(self, standalone):
        # Keep contacts fixed and explicitly qualify the tiny nonzero flux
        # at a free interior carrier node in this controlled test.
        mesh=self.read((self.root/'mesh.json'))
        mesh['nodes'].append(dict(id=4,x=.5,y=.5))
        mesh['triangles']=[dict(id=i,region_id=0,node_ids=[i,(i+1)%4,4]) for i in range(4)]
        mesh['regions'][0]['cell_ids']=list(range(4));self.write('mesh.json',mesh)
        self.input['edge_geometry']=[dict(nodes=[a,b],poisson_F_per_m=1e-10,transport_weight=1.)
            for a,b in ((0,1),(1,2),(2,3),(0,3),(0,4),(1,4),(2,4),(3,4))]
        self.input['boundaries'] += [dict(node=4,kind='psi',value=0.),dict(node=4,kind='temperature',value=300.)]
        for key in ('donors_m3','acceptors_m3','fixed_charge_C_per_m','electron_qf_reference_V','hole_qf_reference_V'):
            self.input[key].append(0.)
        for key in ('state_interleaved','referenced_state_interleaved'):self.input[key] += [0.,0.,0.,300.]
        self.input['silicon_area_m2']=[1/6]*4+[1/3]
        self.input['electrical_gate_solver']['carrier_row_convergence'].update(min_flux_scale=1e-100,scale_floor=1e-30)
        self.input['coordinate_to_metres']=1e-6
        self.input['silicon_area_m2']=[a*1e-12 for a in self.input['silicon_area_m2']]
        for edge in self.input['edge_geometry']:edge['transport_weight']=1.
        self.input.update(diagnostic_pseudo_transient=True,diagnostic_density_projection='v1',
                          diagnostic_pseudo_acceptance='defect_model',diagnostic_near_steady_qf_switch=True,
                          performance_profiling=True)
        if standalone:
            for key in ('diagnostic_pseudo_transient','diagnostic_density_projection',
                        'diagnostic_pseudo_acceptance','diagnostic_near_steady_qf_switch'):
                self.input.pop(key)
            self.input['diagnostic_near_steady_qf_rebase']=True
        self.deck['sweep'].update(max_newton=8,growth_newton=1)
        for amplitude,label in ((0.,'converged'),(5e-14,'near'),(.01,'far')):
            for key in ('state_interleaved','referenced_state_interleaved'):
                self.input[key][17]=amplitude;self.input[key][18]=-amplitude
            self.deck['output_directory']=label;self.write('input.json',self.input)
            self.run_deck()
            data=self.read((self.root/label/'step_0000/output.json'))
            switch=data['near_steady_qf_rebase' if standalone else 'near_steady_qf_switch']
            if standalone:
                self.assertNotIn('pseudo_transient',data)
                self.assertNotIn('density_projection',data)
            if label=='near':
                self.assertTrue(switch['triggered'],str({k:data[k] for k in ('newton_updates','diagnostic_stop','electrical_block_gates','carrier_row_gate')}));self.assertEqual(len(switch['events']),1)
                self.assertLess(switch['events'][0]['merit_before'],1e-9)
                if not standalone:self.assertEqual(data['pseudo_transient']['mass_assemblies'],0)
                self.assertTrue(data['carrier_row_gate']['satisfied'])
                self.assertEqual(data['diagnostic_stop'],'diagnostic_scaled_residual')
                self.assertTrue(all(h['near_steady_qf_active'] and 'pseudo_transient' not in h and 'projection_trials' not in h for h in data['history']))
            elif label=='converged':
                self.assertFalse(switch['triggered']);self.assertEqual(data['newton_updates'],0)
            else:
                self.assertFalse(data['history'][0]['near_steady_qf_active'])
                if not standalone:self.assertGreater(data['pseudo_transient']['mass_assemblies'],0)
        if standalone:
            self.input['diagnostic_density_update_iterations']=1
            self.write('input.json',self.input);self.deck['output_directory']='invalid_rebase'
            self.assertNotEqual(self.run_deck().returncode,0)
            self.assertIn('Near-steady rebase requires',(self.root/'invalid_rebase/step_0000/run.log').read_text())
            return
        self.input['diagnostic_pseudo_transient']=False
        self.input['diagnostic_pseudo_acceptance']='steady';self.write('input.json',self.input)
        self.deck['output_directory']='invalid_switch'
        self.assertNotEqual(self.run_deck().returncode,0)
        self.assertIn('Near-steady switch requires',(self.root/'invalid_switch/step_0000/run.log').read_text())

    def test_voltage_update_limit_rejects_nonpositive_values(self):
        self.input['diagnostic_voltage_update_limit_V']=0.
        self.write('input.json',self.input)
        run=self.run_deck();self.assertNotEqual(run.returncode,0)
    def test_ngmres_does_not_modify_a_converged_point(self):
        self.assertEqual(self.run_deck().returncode,0)
        baseline=self.read((self.root/'output/step_0000/output.json'))
        self.input['diagnostic_ngmres_recovery']=True;self.write('input.json',self.input)
        self.deck['output_directory']='ngmres'
        run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
        result=self.read((self.root/'ngmres/step_0000/output.json'))
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
            result=self.read((self.root/self.deck['output_directory']/'step_0000/output.json'))
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
            ledger=self.read((self.root/self.deck['output_directory']/'ledger.json'))
            results.append(self.read(Path(ledger['accepted_result'])))
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
        result=self.read((self.root/'output/step_0000/output.json'))
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
        zero=self.read((self.root/'output/step_0000/output.json'))
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
        result=self.read((self.root/'tangent/step_0000/output.json'))
        self.assertTrue(result['tangent_preparation']['prepared'])
        self.assertEqual(result['predictor_selection']['selected'],'tangent')
        self.assertEqual(result['performance']['factorizations'],result['newton_updates']+1)
        for a,b in zip(result['state_interleaved'],zero['state_interleaved']):self.assertAlmostEqual(a,b,delta=1e-10)
        source['moving_nodes']=[99];self.write('input.json',self.input)
        self.deck['output_directory']='tangent_fallback'
        run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
        result=self.read((self.root/'tangent_fallback/step_0000/output.json'))
        self.assertFalse(result['tangent_preparation']['prepared'])
        self.assertIn('fallback_reason',result['tangent_preparation'])
    def test_fixed_targets_attempts_exact_target_and_stops_on_failure(self):
        self.deck['sweep'].update(step_policy='fixed_targets',bias_points_V=[0.,.3,.4])
        run=self.run_deck();self.assertEqual(run.returncode,1,run.stderr)
        ledger=self.read((self.root/'output/ledger.json'))
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
            for payload in (self.root/'output/step_0000').iterdir():
                if payload.suffix=='.h5' or payload.name=='input.json':
                    (directory/payload.name).write_bytes(payload.read_bytes())
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
        report=self.read((self.root/'study/summary.json'))
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
        ledger=self.read((self.root/'output/ledger.json'))
        self.assertEqual(len(ledger['runs']),2)
        configs=[self.read((Path(r['directory'])/'input.json')) for r in ledger['runs']]
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
        ledger=self.read((self.root/'output/ledger.json'))
        self.assertFalse(ledger['runs'][0]['prediction']['used'])
        cfg=self.read((Path(ledger['runs'][0]['directory'])/'input.json'))
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
            ledger=self.read((self.root/self.deck['output_directory']/'ledger.json'))
            trajectory=[]
            self.assertGreater(len(ledger['initialization_runs']),3)
            for row in ledger['initialization_runs']:
                result_path=Path(row['result'])
                cfg=self.read((result_path.parent/'input.json'))
                self.assertEqual(cfg['diagnostic_density_update_iterations'],0)
                result=self.read(result_path)
                self.assertTrue(row['gate']['pass_gate'])
                trajectory.append((row['gate_V'],row['solve_mode'],row['newton_updates'],
                                   result['state_interleaved']))
            trajectories.append(trajectory)
        self.assertEqual(trajectories[0],trajectories[1])

    def test_contact_consistency_preserves_neutral_initialization_and_gate_prebias(self):
        self.deck['initialization']=dict(mode='neutral_300K',gate_voltage_V=.1,max_newton=10)
        self.input['skip_equilibrium_poisson_transport']=True
        self.input['boundaries'][-4]['value']=.1
        trajectories=[]
        for enabled in (False,True):
            self.input['diagnostic_near_steady_contact_consistency']=enabled;self.write('input.json',self.input)
            self.deck['output_directory']='contact_init_'+str(enabled)
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            ledger=self.read((self.root/self.deck['output_directory']/'ledger.json'))
            trajectory=[]
            self.assertGreater(len(ledger['initialization_runs']),3)
            for row in ledger['initialization_runs']:
                result_path=Path(row['result']);cfg=self.read((result_path.parent/'input.json'))
                self.assertFalse(cfg['diagnostic_near_steady_contact_consistency'])
                result=self.read(result_path);self.assertNotIn('near_steady_contact_consistency',result)
                trajectory.append((row['gate_V'],row['solve_mode'],result['history'],result['state_interleaved'],result['residual']))
            trajectories.append(trajectory)
            cfg=self.read((Path(ledger['runs'][0]['directory'])/'input.json'))
            self.assertEqual(cfg['diagnostic_near_steady_contact_consistency'],enabled)
        self.assertEqual(trajectories[0],trajectories[1])

    def test_static_preparation_reuse_preserves_initialization_and_resume(self):
        self.deck['reuse_jacobian_structure']=False  # Isolate immutable input preparation.
        self.deck['reuse_linear_analysis']=False
        self.deck['initialization']=dict(mode='neutral_300K',gate_voltage_V=.1,max_newton=10)
        self.input.update(skip_equilibrium_poisson_transport=True,performance_profiling=True,
                          reuse_physics_preparation=True)
        self.input['boundaries'][-4]['value']=.1
        self.write('input.json',self.input)
        trajectories=[]
        for enabled in (False,True):
            self.deck.update(reuse_static_preparation=enabled,output_directory='static_'+str(enabled),resume=False)
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            ledger=self.read((self.root/self.deck['output_directory']/'ledger.json'))
            paths=[Path(r['result']) for r in ledger['initialization_runs']]
            paths += [Path(r['directory'])/'output.json' for r in ledger['runs']]
            trajectory=[]
            for i,path in enumerate(paths):
                result=self.read(path);perf=result.pop('performance')
                self.assertEqual(perf['static_preparation_reused'],enabled and i>0)
                trajectory.append(result)
            trajectories.append(trajectory)
            before=paths[-1].read_bytes();self.deck['resume']=True
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            self.assertEqual(paths[-1].read_bytes(),before)
        self.assertEqual(trajectories[0],trajectories[1])

    def test_linear_context_preserves_initialization_and_resume(self):
        self.deck['initialization']=dict(mode='neutral_300K',gate_voltage_V=.1,max_newton=10)
        self.input.update(skip_equilibrium_poisson_transport=True,performance_profiling=True,
                          reuse_physics_preparation=True)
        self.input['boundaries'][-4]['value']=.1
        self.write('input.json',self.input)
        trajectories=[];analyses=[]
        for enabled in (False,True):
            self.deck.update(output_directory='linear_'+str(enabled),resume=False,reuse_linear_analysis=enabled)
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            ledger=self.read(self.root/self.deck['output_directory']/'ledger.json')
            paths=[Path(r['result']) for r in ledger['initialization_runs']]
            paths += [Path(r['directory'])/'output.json' for r in ledger['runs']]
            trajectory=[];total=0;hits=0
            for path in paths:
                result=self.read(path);perf=result.pop('performance')
                total+=perf['symbolic_analyses'];hits+=int(perf['linear_object_reused'])
                trajectory.append(result)
            self.assertEqual(hits>0,enabled)
            analyses.append(total);trajectories.append(trajectory)
            before=paths[-1].read_bytes();self.deck['resume']=True
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            self.assertEqual(paths[-1].read_bytes(),before)
            self.deck['reuse_linear_analysis']=not enabled
            run=self.run_deck();self.assertNotEqual(run.returncode,0)
            self.assertIn('linear solver policy differs',run.stderr)
        self.assertEqual(trajectories[0],trajectories[1])
        self.assertLess(analyses[1],analyses[0])

    def test_structure_reuse_preserves_initialization_and_resume(self):
        self.deck['reuse_linear_analysis']=False
        self.deck['initialization']=dict(mode='neutral_300K',gate_voltage_V=.1,max_newton=10)
        self.input.update(skip_equilibrium_poisson_transport=True,performance_profiling=True,
                          reuse_physics_preparation=True)
        self.input['boundaries'][-4]['value']=.1
        self.write('input.json',self.input)
        trajectories=[]
        for requested in (None,False,True):
            enabled=True if requested is None else requested
            self.deck.update(output_directory='structure_'+str(requested),resume=False)
            self.deck.pop('reuse_jacobian_structure',None)
            if requested is not None:self.deck['reuse_jacobian_structure']=requested
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            ledger=self.read((self.root/self.deck['output_directory']/'ledger.json'))
            paths=[Path(r['result']) for r in ledger['initialization_runs']]
            paths += [Path(r['directory'])/'output.json' for r in ledger['runs']]
            trajectory=[]
            for i,path in enumerate(paths):
                result=self.read(path);perf=result.pop('performance')
                self.assertEqual(perf['static_preparation_reused'],enabled and i>0)
                self.assertEqual(perf['jacobian_structure_reuse_enabled'],enabled)
                if enabled:
                    self.assertGreater(perf['jacobian_structure_hits']+perf['jacobian_structure_builds'],0)
                trajectory.append(result)
            trajectories.append(trajectory)
            before=paths[-1].read_bytes();self.deck['resume']=True
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            self.assertEqual(paths[-1].read_bytes(),before)
        for trajectory in trajectories[1:]:self.assertEqual(trajectories[0],trajectory)

    def test_neutral_root_newton_preserves_initialization_and_reports_all_work(self):
        self.deck['initialization']=dict(mode='neutral_300K',gate_voltage_V=.1,max_newton=10)
        self.input.update(skip_equilibrium_poisson_transport=True,performance_profiling=True,
                          reuse_physics_preparation=True)
        self.input['boundaries'][-4]['value']=.1
        trajectories=[];iterations=[]
        for enabled in (False,True):
            self.input['diagnostic_neutral_root_newton']=enabled;self.write('input.json',self.input)
            self.deck['output_directory']='root_newton_'+str(enabled)
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            ledger=self.read((self.root/self.deck['output_directory']/'ledger.json'))
            trajectory=[];counts=[0]*5
            paths=[Path(r['result']) for r in ledger['initialization_runs']]
            paths += [Path(r['directory'])/'output.json' for r in ledger['runs']]
            for path in paths:
                cfg=self.read((path.parent/'input.json'))
                self.assertEqual(cfg['diagnostic_neutral_root_newton'],enabled)
                result=self.read(path)
                counts=[a+b for a,b in zip(counts,result['performance']['neutral_root_iteration_counts'])]
                trajectory.append((result['history'],result['state_interleaved'],result['residual']))
            trajectories.append(trajectory);iterations.append(counts)
        self.assertEqual(trajectories[0],trajectories[1])
        self.assertEqual(iterations[0][1],0)
        self.assertGreater(iterations[1][1],0)
        self.assertGreater(iterations[1][4],0)  # Intrinsic near-zero reference uses the legacy guard.

    def test_density_projection_preserves_poisson_prebias(self):
        self.deck['initialization']=dict(mode='neutral_300K',gate_voltage_V=.1,max_newton=10)
        self.input['skip_equilibrium_poisson_transport']=True
        self.input['boundaries'][-4]['value']=.1
        trajectories=[]
        for mode in ('off','v1','v2'):
            self.input['diagnostic_density_projection']=mode;self.write('input.json',self.input)
            self.deck['output_directory']='projection_init_'+mode
            run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
            ledger=self.read((self.root/self.deck['output_directory']/'ledger.json'))
            trajectory=[]
            for row in ledger['initialization_runs']:
                if row['solve_mode']!='poisson':continue
                result=self.read(Path(row['result']))
                self.assertNotIn('density_projection',result)
                self.assertTrue(row['gate']['pass_gate'])
                trajectory.append((row['gate_V'],result['state_interleaved']))
            self.assertGreater(len(trajectory),1);trajectories.append(trajectory)
        self.assertEqual(trajectories[0],trajectories[1]);self.assertEqual(trajectories[0],trajectories[2])

    def test_neutral_initialization_does_not_use_supplied_seed(self):
        self.deck['initialization']=dict(mode='neutral_300K',gate_voltage_V=.1,max_newton=10)
        self.input['skip_equilibrium_poisson_transport']=True
        self.input['boundaries'][-4]['value']=.1
        for key in ('state_interleaved','referenced_state_interleaved'):
            self.input[key]=[999.]*16
        self.write('input.json',self.input)
        run=self.run_deck();self.assertEqual(run.returncode,0,run.stderr)
        ledger=self.read((self.root/'output/ledger.json'))
        self.assertGreater(len(ledger['initialization_runs']),3)
        result=self.read(Path(ledger['accepted_result']))
        self.assertEqual(result['temperature_K'],[300.]*4)
        self.assertAlmostEqual(result['state_interleaved'][12],.1)
        self.assertTrue(all(abs(v)<1e-14 for v in result['state_interleaved'][1::4]))
        self.assertTrue(all(abs(v)<1e-14 for v in result['state_interleaved'][2::4]))


class PredictorComparisonTest(unittest.TestCase):
    def test_contact_sweep_changes_only_flag_and_exact_prefix(self):
        sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
        from run_templates_ldmos_contact_sweep import prepare, FLAG
        original=dict(sentinel=[1,2,3],electrothermal_linear_solver='umfpack')
        deck=dict(input_file='old',output_directory='old',initialization=dict(mode='neutral_300K'),
                  sweep=dict(bias_points_V=[i*40/30 for i in range(31)],growth_newton=12))
        base,bd=prepare(original,deck,Path('/test'),'baseline','first8')
        cand,cd=prepare(original,deck,Path('/test'),'contact','first8')
        self.assertEqual(base,original);self.assertEqual(cand,dict(original,**{FLAG:True}))
        self.assertEqual(cd,bd);self.assertEqual(cd['sweep']['bias_points_V'],deck['sweep']['bias_points_V'][:8])
        self.assertEqual(cd['initialization'],deck['initialization'])
        _,full=prepare(original,deck,Path('/test'),'contact','full')
        self.assertEqual(full['sweep'],deck['sweep']);self.assertNotIn(FLAG,original)
        with self.assertRaises(ValueError):prepare(cand,deck,Path('/test'),'contact','full')

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
