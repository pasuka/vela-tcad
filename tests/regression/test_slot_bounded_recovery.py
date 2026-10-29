"""Bounded recovery must stop wasted work without granting qualification."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
import tempfile
import json
import unittest
from types import SimpleNamespace
from unittest.mock import Mock, patch
import run_slot_ldmos_bv_campaign as c
import run_slot_ldmos_noavalanche_continuation as h
from tests.regression.test_run_slot_ldmos_loadline_continuation import deck


class BoundedRecoveryTests(unittest.TestCase):
    def test_corrector_stops_stagnation_without_accepting_state(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);document=dict(simulation_type='dc_sweep',sweep={'contact':'drain'},
                contacts=[dict(name='drain',bias=0)],solver={})
            def execute(runner,repo,config,run_dir,timeout):
                d=json.loads(config.read_text())
                if d['simulation_type']=='newton_block_step_probe':
                    (run_dir/'trial_carrier_only.h5').write_bytes(b'trial')
                else:
                    (run_dir/'iv.csv').write_text('converged,newton_iterations,failure_reason,carrier_row_violations,global_continuity_closure_satisfied,final_psi_residual_norm\n0,0,line_search_non_decrease,0,1,1.3e-9\n')
                    (run_dir/'rejected.h5').write_bytes(b'rejected')
                    (run_dir/'attempts.csv').write_text('rejected_final_state_file\n'+str(run_dir/'rejected.h5')+'\n')
                return 0,False,.1
            with patch.object(h,'execute_config',side_effect=execute),patch.object(h,'audit_local_candidate',return_value={'passed':True}):
                result=h.correct_local_candidate(document,root/'seed.h5',1.,root/'runner',root,
                    root/'attempt',10,root,8,stagnation_passes=2)
            self.assertFalse(result['passed']);self.assertEqual(result['state_file'],'')
            self.assertEqual(len(result['steps']),2)
            self.assertEqual(result['failure_reason'],'poisson_stagnation_with_qualified_carriers')

    def test_stagnation_does_not_ignore_carrier_or_global_failure_or_progress(self):
        good=dict(failure_reason='line_search_non_decrease',carrier_row_violations='0',
                  global_continuity_closure_satisfied='1',final_psi_residual_norm='1.3e-9')
        self.assertTrue(h.poisson_stagnation([good,good],2))
        for change in [dict(carrier_row_violations='1'), dict(global_continuity_closure_satisfied='0'),
                       dict(final_psi_residual_norm='1e-10'), dict(final_psi_residual_norm='nan'),
                       dict(failure_reason='different')]:
            self.assertFalse(h.poisson_stagnation([good,{**good,**change}],2))
        self.assertFalse(h.poisson_stagnation([good],2))
        self.assertFalse(h.poisson_stagnation([good,good],0))

    def backend(self):
        b=SimpleNamespace(base={},deadline=100.,repo=Path('.'))
        b.set_deadline=lambda x:setattr(b,'deadline',x)
        b.validate=Mock()
        b.solve=Mock(return_value=dict(inner_V=31.2,current_A_per_um=6.5e-15,state_sha256='new'))
        b.solve_predicted=Mock(side_effect=c.ll.NumericalFailure('stalled'))
        return b

    def test_only_one_prediction_and_one_parent_with_qualified_result(self):
        with tempfile.TemporaryDirectory() as tmp:
            b=self.backend();failed=[]
            p,proof=c.actual_outer_candidate(b,dict(inner_V=31.1,state_sha256='seed'),31.106,31.206,
                Path(tmp),previous={'point':dict(state_sha256='old')},predictor_ratio=1.6,
                max_inner_step=.1,predicted_inner=31.2,failed_trials=failed)
            self.assertEqual(b.solve_predicted.call_count,1);self.assertEqual(b.solve.call_count,1)
            b.validate.assert_called_once_with(p)
            self.assertTrue(proof['recovered']);self.assertEqual(len(failed),1)
            self.assertAlmostEqual(proof['actual_outer_V'],31.2065)
            saved=c.load(Path(tmp)/'qualified_device.json')
            self.assertEqual(saved,p)

    def test_failed_voltage_not_repeated_from_same_parent(self):
        with tempfile.TemporaryDirectory() as tmp:
            b=self.backend();b.solve.side_effect=c.ll.NumericalFailure('failure');failed=[]
            for n in range(2):
                with self.assertRaises(c.ll.NumericalFailure):
                    c.actual_outer_candidate(b,dict(inner_V=31.1,state_sha256='seed'),31.106,31.206,
                        Path(tmp)/str(n),previous=None,predictor_ratio=None,max_inner_step=.1,
                        predicted_inner=31.2,failed_trials=failed)
            self.assertEqual(b.solve.call_count,1)

    def test_environment_and_evidence_failures_never_become_numerical_retry(self):
        for error in [OSError('disk'),c.ll.EvidenceError('hash')]:
            with tempfile.TemporaryDirectory() as tmp:
                b=self.backend();b.solve_predicted.side_effect=error
                with self.assertRaises(type(error)):
                    c.actual_outer_candidate(b,dict(inner_V=31.1,state_sha256='seed'),31.106,31.206,
                        Path(tmp),previous={'point':dict(state_sha256='old')},predictor_ratio=1.6,
                        max_inner_step=.1,predicted_inner=31.2,failed_trials=[])
                b.solve.assert_not_called()

    def test_candidate_deadline_is_restored_and_global_budget_remains_terminal(self):
        for now,expected in [(15.,c.ll.NumericalFailure),(100.,c.ll.BudgetExhausted)]:
            b=self.backend()
            with patch.object(c.time,'monotonic',return_value=now), patch.object(c,'candidate',
                    side_effect=c.ll.BudgetExhausted('budget')):
                with self.assertRaises(expected):
                    c.bounded_candidate(b,{},0.,1.,Path('.'),candidate_seconds=10.)
            self.assertEqual(b.deadline,100.)

    def test_actual_supply_still_requires_all_gates_and_independent_confirmation(self):
        for valid,origin in [(True,0.),(False,0.),(True,35.)]:
            with tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);base=deck();base['potential_origin_V']=origin
                config=root/'base.json';seed=root/'seed.json';out=root/'out'
                c.ll.atomic_json(config,base)
                point=dict(id=0,inner_V=1.-origin,current_A_per_um=1e-15,state_sha256='seed')
                c.ll.atomic_json(seed,dict(completed=True,failure='',identity={'config_sha256':'physical'},
                    inner_ledger=[point],outer_ledger=[dict(device_id=0,outer_V=1.001-origin,confirmation={})]))
                args=SimpleNamespace(base_config=config,runner=root/'runner',output=out,repo_root=root,
                    max_seconds=20.,timeout_seconds=1.,max_attempts=1,criterion=1e-7,resume=False,
                    initial_loadline_checkpoint=seed,port_step_control=True,adaptive_initial_step=.1,
                    actual_outer_continuation=True,candidate_max_seconds=10.,outer_goal=2.,
                    reference_origin_floor=origin if origin else None)
                backend=SimpleNamespace(validate=lambda p:None,validate_confirmation=lambda *a:None)
                answer=dict(inner_V=(1.100007 if valid else 1.3)-origin,current_A_per_um=1e-15,state_sha256='new')
                with patch.object(c.ll,'freeze_inputs',return_value=(base,{'files':{},'config_sha256':'physical'})), \
                     patch.object(c.ll,'verify_hashes'),patch.object(c,'backend_for',return_value=backend), \
                     patch.object(c,'bounded_candidate',return_value=(answer,{'evidence':{}})), \
                     patch.object(c,'confirm',return_value=(answer['inner_V']+.001,{'proof':True})) as confirm:
                    c.run(args)
                d=c.load(out/'checkpoint.json')
                self.assertEqual(len(d['records']),2 if valid else 1)
                if valid:
                    confirm.assert_called_once()
                    self.assertAlmostEqual(d['attempts'][0]['proposal_residual_V'],7e-6)
                else:confirm.assert_not_called()


if __name__=='__main__':unittest.main()
