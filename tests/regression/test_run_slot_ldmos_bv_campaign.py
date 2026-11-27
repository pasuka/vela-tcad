"""Campaign safety properties; no heavy numerical solves."""
import importlib.util
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch, Mock

ROOT = Path(__file__).resolve().parents[2]
import sys
sys.path.insert(0, str(ROOT / 'scripts'))
import run_slot_ldmos_bv_campaign as campaign
from tests.regression.test_run_slot_ldmos_loadline_continuation import deck


class CampaignTests(unittest.TestCase):
    def test_repaired_current_is_rechecked_confirmed_and_used_in_feedback(self):
        for current_jump,accepted in ((1e-16,True),(1e-12,False)):
            with tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);config=root/'base.json';source=root/'seed.json';out=root/'out'
                audit=root/'audit.exe';audit.write_bytes(b'audit')
                base=deck();base['potential_origin_V']=0.;campaign.ll.atomic_json(config,base)
                us=(1.,1.1,1.2)
                points=[dict(id=k,inner_V=.4*u,current_A_per_um=.6*u/1e12,state_sha256=str(k)) for k,u in enumerate(us)]
                campaign.ll.atomic_json(source,dict(completed=True,failure='',identity={'config_sha256':'physical'},
                    inner_ledger=points,outer_ledger=[dict(device_id=k,outer_V=u,confirmation={}) for k,u in enumerate(us)]))
                args=SimpleNamespace(base_config=config,runner=root/'runner',output=out,repo_root=root,
                    max_seconds=10.,timeout_seconds=1.,max_attempts=1,criterion=1e-7,resume=False,
                    initial_loadline_checkpoint=source,outer_goal=2.,adaptive_initial_step=.002,
                    progress_guard=True,port_step_control=True,outer_max_step=.02,proposal_max_step=.002,
                    port_predictor='limited_rational',port_feedback_mode='matched_actual_outer',
                    current_growth_control=True,candidate_conservation_runner=audit,row_quality_repair=True)
                b=SimpleNamespace(validate=lambda p:None,validate_confirmation=lambda *a:None)
                def candidate(b,p,o,target,folder,**kw):
                    return dict(p,inner_V=.4*target,current_A_per_um=.6*target/1e12),{'evidence':{}}
                def repair(b,p,proof,folder,control,**kw):
                    self.assertTrue(kw['eligible'])
                    return dict(p,current_A_per_um=p['current_A_per_um']+current_jump),dict(status='improved',elapsed_s=40.)
                with patch.object(campaign.ll,'freeze_inputs',return_value=(base,{'files':{},'config_sha256':'physical'})), \
                     patch.object(campaign.ll,'verify_hashes'),patch.object(campaign,'backend_for',return_value=b), \
                     patch.object(campaign,'candidate',side_effect=candidate),patch.object(campaign,'repair_row_quality',side_effect=repair), \
                     patch.object(campaign.steps,'measure_cost',return_value=dict(iteration_rows=7,history_files=3,elapsed_s=42.)), \
                     patch.object(campaign,'growth_metrics',return_value=dict(candidate_error=.002,candidate_closure=.002,row_ratio=.0002)), \
                     patch.object(campaign,'confirm',side_effect=lambda b,p,f:(p['inner_V']+1e12*p['current_A_per_um'],{})) as confirm:
                    campaign.run(args)
                r=campaign.load(out/'checkpoint.json')
                self.assertEqual(len(r['records']),4 if accepted else 3)
                self.assertEqual(r['identity']['continuation_policy']['row_quality_repair']['max_passes'],1)
                if accepted:
                    last=r['records'][-1]
                    self.assertAlmostEqual(last['outer_V'],1.2021)
                    self.assertEqual(confirm.call_args.args[1]['current_A_per_um'],last['point']['current_A_per_um'])
                    fb=r['attempts'][-1]['step_feedback']
                    self.assertEqual(fb['cost']['elapsed_s'],42.)
                    self.assertAlmostEqual(fb['actual_outer_V'],last['outer_V'])
                else:
                    confirm.assert_not_called()
                    self.assertEqual(r['attempts'][-1]['status'],'numerical_failure')

    def test_row_repair_is_opt_in_and_identity_bound(self):
        self.assertIsNone(campaign.quality_repair_policy(SimpleNamespace(),None))
        with self.assertRaises(ValueError):
            campaign.quality_repair_policy(SimpleNamespace(row_quality_repair=True),None)
        policy=campaign.quality_repair_policy(SimpleNamespace(row_quality_repair=True),campaign.steps.GrowthPolicy())
        self.assertEqual(policy['max_passes'],1)
        old={'files':{str(Path(campaign.__file__).resolve()):'old'},'continuation_policy':{}}
        new={'files':{str(Path(campaign.__file__).resolve()):'new'},'continuation_policy':{'row_quality_repair':policy}}
        self.assertFalse(campaign.allow_driver_update(old,new))

    def repair_fixture(self, root):
        p=dict(inner_V=2.,state_file=str(root/'seed.h5'),state_sha256='seed',evidence={})
        q=dict(p,state_file=str(root/'result.h5'),state_sha256='result',evidence={})
        helper=SimpleNamespace(correct_local_candidate=Mock(return_value=dict(passed=True,state_file='guess.h5')))
        b=SimpleNamespace(base={},runner=root/'runner',repo=root,timeout=1.,restore_poly=True,
            corrector_passes=8,helper=helper,validate=Mock(),adopt_validation=Mock(return_value=q))
        before=dict(candidate_error=.002,candidate_closure=.002,row_ratio=.00099)
        after=dict(before,row_ratio=.0002)
        return p,q,b,before,after

    def test_row_repair_skips_clean_rows_bad_candidates_and_ineligible_predictions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p,q,b,before,after=self.repair_fixture(root)
            for metrics,eligible in [(after,True),(before,False),(dict(before,candidate_error=.006),True),
                                     (dict(before,candidate_closure=.006),True)]:
                with patch.object(campaign,'growth_metrics',return_value=metrics):
                    out,report=campaign.repair_row_quality(b,p,{},root/'repair',campaign.steps.GrowthPolicy(),eligible=eligible)
                    self.assertIs(out,p);self.assertFalse(report['attempted'])
            b.helper.correct_local_candidate.assert_not_called();b.adopt_validation.assert_not_called()

    def test_row_repair_once_validates_and_binds_improved_state(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p,q,b,before,after=self.repair_fixture(root)
            def solve(*args,**kwargs):
                self.assertEqual(b.corrector_passes,8);self.assertTrue(kwargs['same_bias']);return q
            b.adopt_validation.side_effect=solve
            with patch.object(campaign,'growth_metrics',side_effect=[before,after]):
                out,report=campaign.repair_row_quality(b,p,{},root/'repair',campaign.steps.GrowthPolicy(),eligible=True)
            self.assertIs(out,q);self.assertEqual(report['status'],'improved')
            self.assertEqual(b.corrector_passes,8);self.assertEqual(b.adopt_validation.call_count,1)
            self.assertEqual(b.helper.correct_local_candidate.call_args.kwargs['max_passes'],1)
            self.assertIn(str(root/'repair/summary.json'),q['evidence'])
            campaign.ll.verify_hashes(q['evidence'])
            self.assertGreaterEqual(b.validate.call_count,3)

    def test_row_repair_numerical_failure_keeps_original_without_extra_recovery(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p,q,b,before,after=self.repair_fixture(root)
            b.adopt_validation.side_effect=campaign.ll.NumericalFailure('failed full qualification')
            with patch.object(campaign,'growth_metrics',return_value=before):
                out,report=campaign.repair_row_quality(b,p,{},root/'repair',campaign.steps.GrowthPolicy(),eligible=True)
            self.assertIs(out,p);self.assertEqual(report['status'],'numerical_rejection')
            self.assertEqual(b.corrector_passes,8);self.assertEqual(b.adopt_validation.call_count,1)

    def test_row_repair_unqualified_block_does_not_seed_a_solve(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p,q,b,before,after=self.repair_fixture(root)
            b.helper.correct_local_candidate.return_value={'passed':False}
            with patch.object(campaign,'growth_metrics',return_value=before):
                out,report=campaign.repair_row_quality(b,p,{},root/'repair',campaign.steps.GrowthPolicy(),eligible=True)
            self.assertIs(out,p);b.adopt_validation.assert_not_called()
            self.assertEqual(report['status'],'numerical_rejection')

    def test_row_repair_budget_and_evidence_errors_are_not_numerical_fallback(self):
        for error in (campaign.ll.BudgetExhausted('total budget'),campaign.ll.EvidenceError('bad hash'),OSError('environment')):
            with tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);p,q,b,before,after=self.repair_fixture(root)
                b.adopt_validation.side_effect=error
                with patch.object(campaign,'growth_metrics',return_value=before):
                    with self.assertRaises(type(error)):
                        campaign.repair_row_quality(b,p,{},root/'repair',campaign.steps.GrowthPolicy(),eligible=True)
                self.assertEqual(b.corrector_passes,8)
                self.assertEqual(campaign.load(root/'repair/summary.json')['status'],'aborted')

    def test_row_repair_discards_worse_row_or_candidate_error(self):
        for changed in (dict(row_ratio=.001),dict(candidate_error=.006)):
            with tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);p,q,b,before,after=self.repair_fixture(root)
                after.update(changed)
                with patch.object(campaign,'growth_metrics',side_effect=[before,after]):
                    out,report=campaign.repair_row_quality(b,p,{},root/'repair',campaign.steps.GrowthPolicy(),eligible=True)
                self.assertIs(out,p);self.assertEqual(report['status'],'numerical_rejection')

    def test_matched_feedback_cli_parsing_and_invalid_combination(self):
        import contextlib
        import io
        argv=['campaign','--base-config','base.json','--runner','runner.exe',
              '--repo-root','.', '--output','out','--port-feedback-mode','matched_actual_outer']
        with patch.object(sys,'argv',argv),patch.object(campaign,'run') as run, \
             contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as error:campaign.main()
            self.assertEqual(error.exception.code,2);run.assert_not_called()
        with patch.object(sys,'argv',argv+['--port-step-control','--port-predictor','limited_rational']), \
             patch.object(campaign,'run',return_value=0) as run:
            self.assertEqual(campaign.main(),0)
            self.assertEqual(run.call_args.args[0].port_feedback_mode,'matched_actual_outer')

    def test_matched_feedback_requires_hybrid_controller(self):
        self.assertEqual(campaign.feedback_mode(SimpleNamespace(),None),'proposed_outer')
        args=SimpleNamespace(port_feedback_mode='matched_actual_outer')
        for policy in (None,campaign.steps.Policy()):
            with self.assertRaises(ValueError):campaign.feedback_mode(args,policy)
        self.assertEqual(campaign.feedback_mode(args,campaign.steps.Policy(predictor='limited_rational')),
                         'matched_actual_outer')

    def test_matched_feedback_campaign_wiring_identity_and_resume_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);config=root/'base.json';source=root/'seed.json';out=root/'out'
            audit=root/'audit.exe';audit.write_bytes(b'fake audit runner')
            base=deck();base['potential_origin_V']=0.;campaign.ll.atomic_json(config,base)
            points=[dict(id=k,inner_V=.4*u,current_A_per_um=.6*u/1e12,state_sha256=str(k))
                    for k,u in enumerate((1.,1.1,1.2))]
            campaign.ll.atomic_json(source,dict(completed=True,failure='',identity={'config_sha256':'physical'},
                inner_ledger=points,outer_ledger=[dict(device_id=k,outer_V=u,confirmation={})
                for k,u in enumerate((1.,1.1,1.2))]))
            args=SimpleNamespace(base_config=config,runner=root/'runner',output=out,repo_root=root,
                max_seconds=10.,timeout_seconds=1.,max_attempts=1,criterion=1e-7,resume=False,
                initial_loadline_checkpoint=source,outer_goal=2.,adaptive_initial_step=.002,
                progress_guard=True,port_step_control=True,outer_max_step=.02,proposal_max_step=.002,
                port_predictor='limited_rational',port_feedback_mode='matched_actual_outer',
                current_growth_control=True,candidate_conservation_runner=audit)
            b=SimpleNamespace(validate=lambda p:None,validate_confirmation=lambda *a:None)
            def solve(b,p,o,target,folder,**kw):
                actual=target+.01
                return dict(p,inner_V=.4*actual,current_A_per_um=.6*actual/1e12),{'evidence':{}}
            with patch.object(campaign.ll,'freeze_inputs',return_value=(base,{'files':{},'config_sha256':'physical'})), \
                 patch.object(campaign.ll,'verify_hashes'),patch.object(campaign,'backend_for',return_value=b), \
                 patch.object(campaign,'candidate',side_effect=solve) as solver, \
                 patch.object(campaign,'growth_metrics',return_value=dict(candidate_error=.002,candidate_closure=.002,row_ratio=.0005)), \
                 patch.object(campaign,'confirm',side_effect=lambda b,p,f:(p['inner_V']+1e12*p['current_A_per_um'],{})):
                campaign.run(args)
                r=campaign.load(out/'checkpoint.json');fb=r['attempts'][-1]['step_feedback']
                self.assertEqual(len(r['records']),4)
                self.assertGreater(fb['original_normalized_port_prediction_error'],1.)
                self.assertLess(fb['normalized_port_prediction_error'],1e-10)
                self.assertTrue(r['attempts'][-1]['growth_observation']['clean'])
                self.assertEqual(r['identity']['continuation_policy']['port_feedback']['version'],
                                 'matched_actual_outer_v1')
                args.resume=True;args.port_feedback_mode='proposed_outer'
                with self.assertRaisesRegex(campaign.ll.EvidenceError,'identity mismatch'):
                    campaign.run(args)
                self.assertEqual(solver.call_count,1)
                changed=campaign.ll.clone(r['identity'])
                changed['continuation_policy'].pop('port_feedback')
                self.assertFalse(campaign.allow_driver_update(r['identity'],changed))

    def test_growth_cli_state_and_identity_survive_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);config=root/'base.json';source=root/'seed.json';out=root/'out'
            audit=root/'audit.exe';audit.write_bytes(b'fake audit runner')
            base=deck();base['potential_origin_V']=0.;campaign.ll.atomic_json(config,base)
            p=dict(id=0,inner_V=.9,current_A_per_um=1e-13,state_sha256='seed')
            campaign.ll.atomic_json(source,dict(completed=True,failure='',identity={'config_sha256':'physical'},
                inner_ledger=[p],outer_ledger=[dict(device_id=0,outer_V=1.,confirmation={})]))
            args=SimpleNamespace(base_config=config,runner=root/'runner',output=out,repo_root=root,
                max_seconds=10.,timeout_seconds=1.,max_attempts=1,criterion=1e-7,resume=False,
                initial_loadline_checkpoint=source,outer_goal=2.,adaptive_initial_step=.002,
                progress_guard=True,port_step_control=True,outer_max_step=.01,proposal_max_step=.002,
                port_predictor='limited_rational',current_growth_control=True,candidate_conservation_runner=audit)
            b=SimpleNamespace(validate=lambda p:None,validate_confirmation=lambda *a:None)
            def solve(b,p,o,target,folder,**kw):
                i=p['current_A_per_um']+1e-16
                return dict(p,inner_V=target-1e12*i,current_A_per_um=i),{'evidence':{}}
            with patch.object(campaign.ll,'freeze_inputs',return_value=(base,{'files':{},'config_sha256':'physical'})), \
                 patch.object(campaign.ll,'verify_hashes'),patch.object(campaign,'backend_for',return_value=b), \
                 patch.object(campaign,'candidate',side_effect=solve), \
                 patch.object(campaign,'growth_metrics',return_value=dict(candidate_error=.002,candidate_closure=.002,row_ratio=.0005)), \
                 patch.object(campaign,'confirm',side_effect=lambda b,p,f:(p['inner_V']+1e12*p['current_A_per_um'],{})):
                campaign.run(args);args.resume=True;campaign.run(args)
            r=campaign.load(out/'checkpoint.json')
            self.assertEqual(len(r['records']),3)
            self.assertEqual(len(r['growth_state']['observations']),2)
            self.assertEqual(r['identity']['continuation_policy']['current_growth_control']['fractions'],[.04,.08,.20])
            self.assertEqual(r['attempts'][-1]['step_proposal']['growth_policy']['actual_current_fraction_guard'],.25)
            changed=campaign.ll.clone(r['identity']);changed['continuation_policy']['current_growth_control']['fractions'][0]=.05
            self.assertFalse(campaign.allow_driver_update(r['identity'],changed))

    def test_proposal_cap_does_not_replace_actual_acceptance_cap(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); config=root/'base.json'; source=root/'seed.json'; out=root/'out'
            base=deck(); base['potential_origin_V']=0.; campaign.ll.atomic_json(config,base)
            p=dict(id=0,inner_V=.999,current_A_per_um=1e-15,state_sha256='seed')
            campaign.ll.atomic_json(source,dict(completed=True,failure='',identity={'config_sha256':'physical'},
                inner_ledger=[p],outer_ledger=[dict(device_id=0,outer_V=1.,confirmation={})]))
            args=SimpleNamespace(base_config=config,runner=root/'runner',output=out,repo_root=root,
                max_seconds=10.,timeout_seconds=1.,max_attempts=2,criterion=1e-7,resume=False,
                initial_loadline_checkpoint=source,outer_goal=1.1,adaptive_initial_step=.009,
                progress_guard=True,port_step_control=True,outer_max_step=.01,proposal_max_step=.002)
            b=SimpleNamespace(validate=lambda p:None,validate_confirmation=lambda *a:None)
            targets=[]
            def solve(b,p,o,target,folder,**kw):
                targets.append(target-o)
                # Strict correction exceeds proposal cap, but respects actual cap.
                return dict(p,inner_V=target-.001+.0003),{'evidence':{}}
            with patch.object(campaign.ll,'freeze_inputs',return_value=(base,{'files':{},'config_sha256':'physical'})), \
                 patch.object(campaign.ll,'verify_hashes'),patch.object(campaign,'backend_for',return_value=b), \
                 patch.object(campaign,'candidate',side_effect=solve), \
                 patch.object(campaign,'confirm',side_effect=lambda b,p,f:(p['inner_V']+1e12*p['current_A_per_um'],{})):
                campaign.run(args)
            r=campaign.load(out/'checkpoint.json')
            self.assertEqual(len(r['records']),3)
            self.assertTrue(all(abs(h-.002)<1e-12 for h in targets))
            self.assertEqual(r['identity']['continuation_policy']['proposal_max_step_V'],.002)
            changed=campaign.ll.clone(r['identity']); changed['continuation_policy']['proposal_max_step_V']=.003
            self.assertFalse(campaign.allow_driver_update(r['identity'],changed))

    def test_progress_guard_rejects_archived_drift_but_accepts_vertical_current_growth(self):
        policy=campaign.progress_policy(SimpleNamespace(progress_guard=True))
        p=dict(inner_V=2.400112295021657,current_A_per_um=6.038377250647024e-14)
        q=dict(inner_V=2.400112295021745,current_A_per_um=6.038399986535761e-14)
        outer=lambda x:x['inner_V']+1e12*x['current_A_per_um']
        self.assertFalse(campaign.assess_progress(p,q,outer(p),outer(q),policy)['resolved'])
        q=dict(p,current_A_per_um=1.1*p['current_A_per_um'])
        self.assertTrue(campaign.assess_progress(p,q,outer(p),outer(q),policy)['resolved'])

    def test_stagnation_retains_proof_and_retries_proposal_without_history_pollution(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);config=root/'base.json';source=root/'seed.json';out=root/'out'
            base=deck();base['potential_origin_V']=0.;campaign.ll.atomic_json(config,base)
            p=dict(id=0,inner_V=.999,current_A_per_um=1e-15,state_sha256='seed')
            campaign.ll.atomic_json(source,dict(completed=True,failure='',identity={'config_sha256':'physical'},
                inner_ledger=[p],outer_ledger=[dict(device_id=0,outer_V=1.,confirmation={})]))
            args=SimpleNamespace(base_config=config,runner=root/'runner',output=out,repo_root=root,
                max_seconds=10.,timeout_seconds=1.,max_attempts=2,criterion=1e-7,resume=False,
                initial_loadline_checkpoint=source,outer_goal=1.03,adaptive_initial_step=.026,
                progress_guard=True,port_step_control=True)
            b=SimpleNamespace(validate=lambda p:None,validate_confirmation=lambda *a:None)
            targets=[]
            def solve(b,p,o,target,folder,**kw):
                targets.append(target)
                if len(targets)==1:return dict(p,current_A_per_um=1.0002e-15),{'evidence':{}}
                return dict(p,inner_V=target-.001),{'evidence':{}}
            with patch.object(campaign.ll,'freeze_inputs',return_value=(base,{'files':{},'config_sha256':'physical'})), \
                 patch.object(campaign.ll,'verify_hashes'),patch.object(campaign,'backend_for',return_value=b), \
                 patch.object(campaign,'candidate',side_effect=solve), \
                 patch.object(campaign,'confirm',side_effect=lambda b,p,f:(p['inner_V']+1e12*p['current_A_per_um'],{})):
                campaign.run(args)
            r=campaign.load(out/'checkpoint.json')
            self.assertEqual(r['attempts'][0]['status'],'stagnation')
            self.assertIn('stalled_point',r['attempts'][0])
            self.assertNotIn('step_feedback',r['attempts'][0])
            self.assertAlmostEqual(targets[1],1.013)
            self.assertEqual(len(r['records']),2)
            self.assertEqual(r['attempts'][1]['seed_sha256'],'seed')

    def test_zero_rhs_is_candidate_only_and_keeps_strict_tolerance(self):
        for enabled in (False, True):
            with self.subTest(enabled=enabled), tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);base=deck();original=campaign.ll.clone(base);captured=[]
                def render(*args):
                    return campaign.ll.clone(base),root/'candidate.json',root/'run'
                def write(rendered):
                    captured.append(rendered[0])
                    raise campaign.ll.NumericalFailure('stop before solver')
                backend=SimpleNamespace(base=base,helper=SimpleNamespace(render_attempt=render),_write_run=write)
                with self.assertRaises(campaign.ll.NumericalFailure):
                    campaign.candidate(backend,dict(inner_V=38.5,state_file='seed'),38.6,38.61,root,
                        zero_converged_device_rhs=enabled)
                circuit=captured[0]['sweep']['external_circuit']
                self.assertEqual(captured[0]['sweep']['bias_points'],[38.61])
                self.assertEqual(circuit['initial_inner_voltage_V'],38.5)
                self.assertEqual(circuit['coupled_equation_tolerance'],1e-8)
                if enabled:
                    self.assertTrue(circuit['coupled_inexact_device_forcing']['zero_converged_residual'])
                    self.assertEqual(circuit['coupled_inexact_device_forcing']['max_equation_tolerance'],1e-8)
                else:
                    self.assertNotIn('coupled_inexact_device_forcing',circuit)
                self.assertEqual(base,original)
                self.assertEqual(captured[0]['solver'],base['solver'])

    def test_tangent_starts_at_parent_and_preserves_strict_base(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);base=deck();original=campaign.ll.clone(base);captured=[]
            def render(*args):
                self.assertEqual(args[3], Path('parent.h5'))
                return campaign.ll.clone(base),root/'candidate.json',root/'run'
            def write(rendered):
                captured.append(rendered[0])
                raise campaign.ll.NumericalFailure('stop before solver')
            backend=SimpleNamespace(base=base,helper=SimpleNamespace(render_attempt=render),_write_run=write)
            with self.assertRaises(campaign.ll.NumericalFailure):
                campaign.candidate(backend,dict(inner_V=38.5,state_file='parent.h5'),38.6,38.7,root,
                    state_predictor='tangent',previous={'not':'a secant parent'},predictor_ratio=1.6)
            circuit=captured[0]['sweep']['external_circuit']
            self.assertEqual(circuit['coupled_state_predictor'],'tangent')
            self.assertEqual(circuit['initial_inner_voltage_V'],38.5)
            self.assertEqual(captured[0]['sweep']['bias_points'],[38.7])
            self.assertEqual(base,original)

    def test_candidate_policy_change_cannot_be_hidden_as_driver_update(self):
        path=str(Path(campaign.__file__).resolve())
        old={'files':{path:'a'},'continuation_policy':{}}
        for policy in ({'candidate_zero_converged_device_rhs':True},
                       {'progress_guard':campaign.progress_policy(SimpleNamespace(progress_guard=True))},
                       {'coupled_candidate_protocol':'matched_inner_direct_target_v1'},
                       {'candidate_conservation':'frozen_candidate_conservation_v1'}):
            new={'files':{path:'b'},'continuation_policy':policy}
            self.assertFalse(campaign.allow_driver_update(old,new))

    def test_port_controller_rejects_oversized_acceptance_before_ledger(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);config=root/'base.json';source=root/'seed.json';out=root/'out'
            base=deck();base['potential_origin_V']=0.
            campaign.ll.atomic_json(config,base)
            point=dict(id=0,inner_V=.999,current_A_per_um=1e-15,state_sha256='seed')
            campaign.ll.atomic_json(source,dict(completed=True,failure='',identity={'config_sha256':'physical'},
                inner_ledger=[point],outer_ledger=[dict(device_id=0,outer_V=1.,confirmation={})]))
            args=SimpleNamespace(base_config=config,runner=root/'runner',output=out,repo_root=root,
                max_seconds=10.,timeout_seconds=1.,max_attempts=1,criterion=1e-7,resume=False,
                initial_loadline_checkpoint=source,port_step_control=True,adaptive_initial_step=.1,
                accepted_inner_max_step=.2,outer_max_step=.5,trial_inner_max_step=.007,
                candidate_zero_converged_device_rhs=True,candidate_state_predictor='tangent')
            backend=SimpleNamespace(validate=lambda p:None,validate_confirmation=lambda *a:None)
            answer=dict(inner_V=1.299,current_A_per_um=1e-15,state_sha256='new')
            with patch.object(campaign.ll,'freeze_inputs',return_value=(base,{'files':{},'config_sha256':'physical'})), \
                 patch.object(campaign.ll,'verify_hashes'),patch.object(campaign,'backend_for',return_value=backend), \
                 patch.object(campaign,'candidate',return_value=(answer,{'evidence':{}})) as solve, \
                 patch.object(campaign,'confirm') as confirm:
                self.assertEqual(campaign.run(args),1)
                confirm.assert_not_called()
            self.assertEqual(solve.call_args.kwargs['max_inner_step'],.007)
            self.assertEqual(solve.call_args.kwargs['coupled_inner_step'],.007)
            self.assertTrue(solve.call_args.kwargs['zero_converged_device_rhs'])
            self.assertEqual(solve.call_args.kwargs['state_predictor'],'tangent')
            result=campaign.load(out/'checkpoint.json')
            self.assertTrue(result['identity']['continuation_policy']['candidate_zero_converged_device_rhs'])
            self.assertEqual(result['identity']['continuation_policy']['candidate_state_predictor'],'bordered_tangent_v1')
            self.assertEqual(len(result['records']),1)
            self.assertEqual(result['attempts'][0]['status'],'numerical_failure')
            self.assertEqual(result['next_span_V'],.05)

    def test_policy_is_opt_in_and_resume_bound(self):
        self.assertIsNone(campaign.step_policy(SimpleNamespace()))
        policy=campaign.step_policy(SimpleNamespace(port_step_control=True))
        self.assertEqual((policy.outer_max_V,policy.accepted_inner_max_V,policy.trial_inner_max_V),(.2,.2,.1))
        path=str(Path(campaign.__file__).resolve())
        old={'files':{path:'a'},'continuation_policy':{'outer_max_V':.1}}
        new={'files':{path:'b'},'continuation_policy':{'outer_max_V':.2}}
        self.assertFalse(campaign.allow_driver_update(old,new))

    def test_high_prediction_retries_numerical_failure_but_not_environment(self):
        for error in (campaign.ll.NumericalFailure('rejected'), OSError('disk error')):
            with self.subTest(error=type(error).__name__), tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp); previous_file=root/'previous.h5';seed=root/'seed.h5'
                previous_file.write_bytes(b'previous');seed.write_bytes(b'seed')
                point=dict(inner_V=38.51,state_file=str(seed))
                previous=dict(point=dict(inner_V=38.50,state_file=str(previous_file)),outer_V=100.)
                helper=SimpleNamespace(write_secant_restart=lambda a,b,c,*args,**kw:c.write_bytes(b'predicted'))
                backend=SimpleNamespace(base={},repo=root,helper=helper,validate=lambda p:None)
                execute=campaign.candidate
                with patch.object(campaign,'candidate',side_effect=[error,({'ok':True},{'evidence':{}})]) as nested:
                    if isinstance(error,OSError):
                        with self.assertRaises(OSError):
                            execute(backend,point,200.,300.,root,previous=previous,predictor_ratio=1.6,
                                    zero_converged_device_rhs=True)
                        self.assertEqual(nested.call_count,1)
                    else:
                        result,proof=execute(backend,point,200.,300.,root,previous=previous,predictor_ratio=1.6,
                                             zero_converged_device_rhs=True)
                        self.assertTrue(result['ok']);self.assertTrue(proof['recovered'])
                        self.assertEqual(nested.call_args_list[1].args[1]['state_file'],str(seed))
                        self.assertTrue(proof['evidence'])
                    self.assertTrue(all(call.kwargs['zero_converged_device_rhs'] for call in nested.call_args_list))
                    predicted=campaign.load(root/'predictor.json')
                    self.assertEqual(nested.call_args_list[0].args[1]['inner_V'],predicted['target_inner_V'])

    def test_stops_at_current_or_physical_outer_goal(self):
        self.assertEqual(campaign.stop_reason(-1e-7, 1e5, 1e-7, 2e5), 'current_criterion')
        self.assertEqual(campaign.stop_reason(1e-8, 2e5, 1e-7, 2e5), 'outer_voltage_goal')
        self.assertIsNone(campaign.stop_reason(1e-8, 199999., 1e-7, 2e5))
        with self.assertRaises(campaign.ll.EvidenceError):
            campaign.stop_reason(float('nan'), 1., 1e-7, 2e5)

    def test_predictor_never_crosses_reference_frames(self):
        records = [dict(base={'potential_origin_V': 0}, point={'inner_V': 2.9}),
                   dict(base={'potential_origin_V': 5}, point={'inner_V': -1.9})]
        self.assertIsNone(campaign.predictor_parent(records, {'potential_origin_V': 5},
            {'inner_V': -2.}, .1, 1.6))
        self.assertIs(campaign.predictor_parent(records, {'potential_origin_V': 0},
            {'inner_V': 3.}, .1, 1.6), records[0])

    def test_outer_history_selects_parent_across_internal_turn(self):
        base={'potential_origin_V':36.}
        records=[dict(base=base,outer_V=u,point=dict(inner_V=v,current_A_per_um=i))
                 for u,v,i in [(100.,2.50512,9.75e-11),(125.,2.50511,1.225e-10)]]
        now=dict(inner_V=2.50510,current_A_per_um=1.475e-10)
        self.assertIsNone(campaign.predictor_parent(records,base,now,20.,1.6))
        self.assertIs(campaign.predictor_parent(records,base,now,20.,1.6,outer=150.),records[-1])
        self.assertIsNone(campaign.predictor_parent(records,{'potential_origin_V':35.},now,20.,1.6,outer=150.))

    def test_current_target_cli_fold_identity_and_retry(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);config=root/'base.json';source=root/'seed.json';out=root/'out'
            audit=root/'audit.exe';audit.write_bytes(b'audit')
            base=deck();base['potential_origin_V']=0.;campaign.ll.atomic_json(config,base)
            p=dict(id=0,inner_V=.9,current_A_per_um=1e-10,state_sha256='seed')
            campaign.ll.atomic_json(source,dict(completed=True,failure='',identity={'config_sha256':'physical'},
                inner_ledger=[p],outer_ledger=[dict(device_id=0,outer_V=100.9,confirmation={})]))
            args=SimpleNamespace(base_config=config,runner=root/'runner',output=out,repo_root=root,
                max_seconds=10.,timeout_seconds=1.,max_attempts=2,criterion=1e-7,resume=False,
                initial_loadline_checkpoint=source,outer_goal=1000.,adaptive_initial_step=.001,
                progress_guard=True,port_step_control=True,outer_max_step=.02,
                port_predictor='limited_rational',current_growth_control=True,candidate_conservation_runner=audit,
                current_proposal_mode='current_target',accepted_inner_decrease_max_step=20e-6)
            b=SimpleNamespace(validate=lambda p:None,validate_confirmation=lambda *a:None)
            spans=[]
            def solve(b,p,o,target,folder,**kw):
                self.assertTrue(kw['predictor_outer_coordinate'])
                spans.append(target-o)
                if len(spans)==1:raise campaign.ll.NumericalFailure('synthetic candidate rejection')
                v=p['inner_V']-3e-6
                return dict(p,inner_V=v,current_A_per_um=(target-v)/1e12),{'evidence':{}}
            with patch.object(campaign.ll,'freeze_inputs',return_value=(base,{'files':{},'config_sha256':'physical'})), \
                 patch.object(campaign.ll,'verify_hashes'),patch.object(campaign,'backend_for',return_value=b), \
                 patch.object(campaign,'candidate',side_effect=solve), \
                 patch.object(campaign,'growth_metrics',return_value=dict(candidate_error=.002,candidate_closure=.002,row_ratio=.0005)), \
                 patch.object(campaign,'confirm',side_effect=lambda b,p,f:(p['inner_V']+1e12*p['current_A_per_um'],{})):
                campaign.run(args)
            r=campaign.load(out/'checkpoint.json');identity=r['identity']['continuation_policy']
            self.assertEqual(len(r['records']),2)
            self.assertAlmostEqual(spans[0],4.)
            self.assertAlmostEqual(spans[1],2.)
            self.assertEqual(identity['current_growth_control']['version'],'current_target_fold_v1')
            self.assertEqual(identity['port_step_control']['accepted_inner_decrease_max_V'],20e-6)
            self.assertEqual(identity['predictor_history_coordinate'],'outer_current_v1')
            changed=campaign.ll.clone(r['identity'])
            changed['continuation_policy']['port_step_control']['accepted_inner_decrease_max_V']=30e-6
            self.assertFalse(campaign.allow_driver_update(r['identity'],changed))

    def test_scalar_parent_import_requires_all_rows_and_circuit_confirmation(self):
        point = dict(id=5, inner_V=3., active_ids={'electron':[1], 'hole':[1]})
        previous = dict(id=4, inner_V=2.9, active_ids=point['active_ids'])
        calls = []
        control = SimpleNamespace(snapshot={'inner_ledger':[]}, save=lambda: None,
            accept_outer=lambda p,o: calls.append((p['id'],o)))
        backend = SimpleNamespace(validate=lambda p: calls.append(p['id']))
        campaign.seed_scalar(control, backend, point, previous, 3.001)
        self.assertEqual(calls,[4,5,(1,3.001)])
        self.assertEqual(point['id'],5)
        bad = {**previous, 'active_ids':{'electron':[], 'hole':[]}}
        with self.assertRaises(campaign.ll.EvidenceError):
            campaign.seed_scalar(control, backend, point, bad, 3.)

    def test_seeded_campaign_clips_goal_and_finishes_on_last_budgeted_attempt(self):
        for by_current in (False, True):
            with self.subTest(by_current=by_current), tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp); config=root/'base.json'; source=root/'seed.json'; out=root/'out'
                base=deck();base['potential_origin_V']=0.
                campaign.ll.atomic_json(config,base)
                point=dict(id=0,inner_V=2.999,current_A_per_um=1e-15,state_sha256='seed')
                campaign.ll.atomic_json(source,dict(completed=True,failure='',identity={'config_sha256':'physical'},
                    inner_ledger=[point],outer_ledger=[dict(device_id=0,outer_V=3.,confirmation={})]))
                args=SimpleNamespace(base_config=config,runner=root/'runner',output=out,repo_root=root,
                    max_seconds=10.,timeout_seconds=1.,max_attempts=1,criterion=1e-7,resume=False,
                    initial_loadline_checkpoint=source,outer_goal=3.01,state_predictor_max_ratio=1.6,
                    adaptive_initial_step=.025,adaptive_growth=1.4,low_voltage_max_step=.1)
                args.port_step_control = not by_current
                backend=SimpleNamespace(validate=lambda p:None,validate_confirmation=lambda *a:None)
                seen=[]
                def solve(b,p,o,target,folder,**kw):
                    seen.append((target,kw))
                    i=1e-7 if by_current else 1e-15
                    return dict(inner_V=target-1e12*i,current_A_per_um=i,state_sha256='new'),{'evidence':{}}
                def confirm(b,p,f):return p['inner_V']+1e12*p['current_A_per_um'],{}
                with patch.object(campaign.ll,'freeze_inputs',return_value=(base,{'files':{},'config_sha256':'physical'})), \
                     patch.object(campaign.ll,'verify_hashes'),patch.object(campaign,'backend_for',return_value=backend), \
                     patch.object(campaign,'candidate',side_effect=solve),patch.object(campaign,'confirm',side_effect=confirm):
                    self.assertEqual(campaign.run(args),0)
                self.assertAlmostEqual(seen[0][0],3.01)
                self.assertEqual(seen[0][1]['predictor_ratio'],1.6)
                result=campaign.load(out/'checkpoint.json')
                self.assertEqual(result['stop_reason'],'current_criterion' if by_current else 'outer_voltage_goal')
                self.assertEqual(len(result['attempts']),1)
                if not by_current:
                    self.assertAlmostEqual(result['next_span_V'],.014)
                    self.assertIn('port_step_control',result['identity']['continuation_policy'])

    def test_driver_update_never_allows_physical_resource_change(self):
        path = str(Path(campaign.__file__).resolve())
        old = {'files': {path: 'old', 'runner': 'same', 'mesh': 'same'}, 'criterion': 1e-7}
        new = campaign.ll.clone(old)
        new['files'][path] = 'new'
        self.assertTrue(campaign.allow_driver_update(old, new))
        new['files']['mesh'] = 'changed'
        self.assertFalse(campaign.allow_driver_update(old, new))

    def test_current_guard_does_not_divide_equilibrium_noise(self):
        campaign.check_step(0., 1e-16, .025, 0.)
        campaign.check_step(-1e-24, 1e-16, .025, 0.)
        with self.assertRaises(campaign.ll.NumericalFailure):
            campaign.check_step(1e-10, 2.1e-10, 200., 100.)
        with self.assertRaises(campaign.ll.NumericalFailure):
            campaign.check_step(1e-10, 1.1e-10, 100., 100.)

    def test_nonfinite_coordinate_is_evidence_failure(self):
        with self.assertRaises(campaign.ll.EvidenceError):
            campaign.check_step(0., float('nan'), 1., 0.)

    def test_small_resistor_drop_root_stays_bracketed_without_forced_bisection(self):
        # A linear device has an exact root very close to the high endpoint.
        samples = [(0., -.025), (.025, .00005)]
        value, bracketed = campaign.ll.proposal(samples, .025, .025, 0, 1e-6)
        self.assertTrue(bracketed)
        self.assertGreater(value, 0.)
        self.assertLess(value, .025)
        self.assertAlmostEqual(1.002 * value - .025, 0., places=15)
        default, _ = campaign.ll.proposal(samples, .025, .025, 0)
        self.assertEqual(default, .0125)
        for bad in (0., .5, float('nan')):
            with self.assertRaises(ValueError):
                campaign.ll.proposal(samples, .025, .025, 0, bad)

    def test_frame_change_preserves_physical_contacts_and_source(self):
        base = deck()
        base.update(contacts=[dict(name=n, bias=-36., type='ohmic') for n in
                              ('source', 'drain', 'gate', 'SLOT', 'substrate')],
                    potential_origin_V=36.)
        base['sweep']['contact'] = 'drain'
        base['solver']['impact_ionization'] = dict(coupling_mode='self_consistent')
        cold = campaign.frame.translate_document(base, 0.)
        self.assertTrue(all(c['bias'] == 0 for c in cold['contacts']))
        self.assertEqual(cold['solver'], base['solver'])

    def test_high_current_step_uses_resistor_scale(self):
        self.assertEqual(campaign.proposed_span(38.52, 1e-8), 4000.)
        self.assertLessEqual(campaign.proposed_span(0., 0.), .025)

    def test_low_and_high_candidates_use_valid_separate_solver_options(self):
        for bias, expected in ((38.4, 'coupled_newton'),):
            base = deck()
            rendered = []
            def stop_before_solve(request):
                rendered.append(request[0])
                raise StopIteration
            backend = SimpleNamespace(base=base,
                helper=SimpleNamespace(render_attempt=lambda *a: (campaign.ll.clone(base), Path('config'), Path('run'))),
                _write_run=stop_before_solve)
            with self.assertRaises(StopIteration):
                campaign.candidate(backend, dict(inner_V=bias, state_file='seed'), bias, bias+.025, Path('unused'),coupled_inner_step=.007)
            circuit = rendered[0]['sweep']['external_circuit']
            self.assertTrue(circuit['enabled'])
            self.assertEqual(circuit['solver'], expected)
            if bias == 0:
                self.assertFalse(any(k.startswith('coupled_') for k in circuit))
            else:
                self.assertEqual(circuit['coupled_linear_solver'], 'direct_bordered_qr')
                self.assertEqual(circuit['max_inner_voltage_step_V'], .007)

    def test_low_candidate_uses_audited_scalar_controller(self):
        backend = SimpleNamespace(base=deck(), repo=ROOT, runner=ROOT/'runner', deadline=1e20)
        with patch.object(campaign.ll, 'freeze_inputs', return_value=(deck(), {'files': {}})), \
             patch.object(campaign.ll, 'Campaign', side_effect=StopIteration) as controller:
            with self.assertRaises(StopIteration):
                campaign.candidate(backend, dict(inner_V=0., state_file='seed'), 0., .025, Path('unused'),max_inner_step=.007)
            self.assertEqual(controller.call_args.args[2], [.025])
            self.assertEqual(controller.call_args.kwargs['initial_outer'], 0.)
            self.assertEqual(controller.call_args.kwargs['max_inner_step'], .007)

    def test_circuit_confirmation_is_explicitly_validated(self):
        calls = []
        backend = SimpleNamespace(confirm=lambda p,o,f: {'marker': o},
            validate_confirmation=lambda p,o,c: calls.append((o,c)))
        outer, proof = campaign.confirm(backend,
            dict(inner_V=38., current_A_per_um=1e-14), Path('unused'))
        self.assertAlmostEqual(outer, 38.01)
        self.assertEqual(calls, [(outer, proof)])

    def test_resume_rejects_changed_resource_before_any_solve(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = root / 'base.json'
            campaign.ll.atomic_json(config, deck())
            campaign.ll.atomic_json(root / 'checkpoint.json', dict(identity={'wrong': True}))
            args = SimpleNamespace(base_config=config, runner=root/'runner', output=root,
                repo_root=root, max_seconds=1., criterion=1e-7, resume=True)
            with patch.object(campaign.ll, 'freeze_inputs', return_value=(deck(), {'files': {}})), \
                 patch.object(campaign, 'backend_for') as backend:
                with self.assertRaisesRegex(campaign.ll.EvidenceError, 'identity mismatch'):
                    campaign.run(args)
                backend.assert_not_called()


if __name__ == '__main__':
    unittest.main()
