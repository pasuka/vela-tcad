"""Physical load-line identities and bounded controller behavior; no solver."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import slot_loadline_step_control as s


def records(xs=(1.,1.08,1.2), slope=2e-12, origin=0.):
    return [dict(base={'potential_origin_V':origin}, outer_V=x+1e12*(1e-12+slope*x)-origin,
                 point=dict(inner_V=x-origin,current_A_per_um=1e-12+slope*x)) for x in xs]


class StepControlTests(unittest.TestCase):
    def test_current_target_does_not_recycle_shrunken_actual_span(self):
        rs=records();p=s.Policy(predictor='limited_rational',outer_max_V=.02)
        control=s.GrowthPolicy(proposal_mode='current_target')
        state=s.initial_growth_state();state['stage']=2
        d=s.growth_propose(rs,p,control,state,.001,100.)
        self.assertAlmostEqual(d['span_V'],.2*1e12*rs[-1]['point']['current_A_per_um'])
        actual=.68*d['span_V'];par=d['parent']
        fb=s.feedback(d,par['inner_V'],par['current_A_per_um']+actual/1e12,
            actual,p,dict(iteration_rows=2,elapsed_s=20.),actual_reference_origin=0.)
        self.assertEqual(fb['proposal_feedback_factor'],1.)
        self.assertAlmostEqual(fb['next_span_V'],d['span_V'])
        self.assertNotIn('invalid_model_no_growth',fb['factors'])
        state['proposal_feedback_factor']=fb['proposal_feedback_factor']
        self.assertEqual(s.growth_propose(rs,p,control,state,fb['next_span_V'],100.)['span_V'],d['span_V'])

    def test_current_target_fold_and_unchanged_final_guards(self):
        rs=[dict(base={'potential_origin_V':0.},outer_V=v+1e12*i,
                 point=dict(inner_V=v,current_A_per_um=i))
            for v,i in [(38.5051216,1.17e-10),(38.5051205,1.75e-10),(38.5051091,2.62e-10)]]
        p=s.Policy(predictor='limited_rational',outer_max_V=.02,accepted_inner_max_V=.1,
                   accepted_inner_decrease_max_V=20e-6)
        state=s.initial_growth_state();state['stage']=2
        d=s.growth_propose(rs,p,s.GrowthPolicy(proposal_mode='current_target'),state,1e-6,1e5)
        v=rs[-1]['point']['inner_V']-3e-6;i=rs[-1]['point']['current_A_per_um']*1.2
        self.assertIsNone(s.check_actual(d,v,v+1e12*i,i,p))
        self.assertIsNotNone(s.check_actual(d,v-.001,v+1e12*i,i,p))
        self.assertIsNotNone(s.check_actual(d,v,v+1e12*i,rs[-1]['point']['current_A_per_um']*1.251,p))
        old=s.replace(p,accepted_inner_decrease_max_V=0.)
        self.assertIsNotNone(s.check_actual(d,v,v+1e12*i,i,old))

    def test_current_target_preserves_retry_quality_cost_and_goal_caps(self):
        rs=records();p=s.Policy(predictor='limited_rational',outer_max_V=.02)
        control=s.GrowthPolicy(proposal_mode='current_target');state=s.initial_growth_state()
        for kwargs in ({'retry_span':.002},{'proposal_cap':.002}):
            self.assertAlmostEqual(s.growth_propose(rs,p,control,state,1.,100.,**kwargs)['span_V'],.002)
        self.assertEqual(s.growth_propose(rs,p,control,state,1.,.001)['span_V'],.001)
        state['observations']=[dict(candidate_error=.01,span_V=.002)]
        self.assertAlmostEqual(s.growth_propose(rs,p,control,state,1.,100.)['span_V'],.001)
        d=s.growth_propose(rs,p,control,s.initial_growth_state(),1.,100.);par=d['parent']
        fb=s.feedback(d,par['inner_V'],par['current_A_per_um']+d['span_V']/1e12,
            d['span_V'],p,dict(iteration_rows=100,elapsed_s=600.),actual_reference_origin=0.)
        self.assertEqual(fb['proposal_feedback_factor'],.5)
        with self.assertRaises(ValueError):
            s.Policy(accepted_inner_decrease_max_V=float('nan'))

    def matched_case(self, law=lambda u:.8*u/(1+.1*u)):
        rs=self.hybrid_records([1,1.2,1.4,1.6,1.8,2],law)
        policy=s.Policy(predictor='limited_rational',current_fraction=1.)
        proposal=s.propose(rs,policy,.05,10.)
        u=2.07
        return proposal,dict(proposal=proposal,inner=law(u),current=(u-law(u))/1e12,
            actual_span=.07,policy=policy,cost=dict(iteration_rows=2,elapsed_s=40.),
            actual_reference_origin=5.)

    def test_matched_feedback_separates_supply_migration_from_model_error(self):
        import copy
        d,kw=self.matched_case();saved=copy.deepcopy(d)
        old=s.feedback(**kw)
        result=s.feedback(**kw,mode='matched_actual_outer')
        self.assertGreater(old['normalized_port_prediction_error'],.1)
        self.assertLess(result['normalized_port_prediction_error'],1e-10)
        self.assertEqual(result['original_normalized_port_prediction_error'],old['normalized_port_prediction_error'])
        self.assertAlmostEqual(result['outer_translation_V'],.02)
        self.assertTrue(result['promotion_eligible'])
        self.assertEqual(d,saved)  # no fitting to the accepted point
        self.assertEqual(old,s.feedback(**kw,mode='proposed_outer'))

    def test_matched_feedback_detects_actual_model_error_and_preserves_cost_limits(self):
        d,kw=self.matched_case()
        # Preserve exact load line but move the device off the predicted curve.
        bad=dict(kw,inner=kw['inner']+.03,current=kw['current']-.03/1e12)
        result=s.feedback(**bad,mode='matched_actual_outer')
        self.assertGreater(result['normalized_port_prediction_error'],1)
        self.assertLess(result['factor'],1.)
        self.assertEqual(s.feedback(**kw,mode='matched_actual_outer',recovered=True)['factor'],.5)
        costly=dict(kw,cost=dict(iteration_rows=80,elapsed_s=400.))
        self.assertEqual(s.feedback(**costly,mode='matched_actual_outer')['factor'],.5)

    def test_matched_feedback_vertical_branch_and_real_growth_guard(self):
        d,kw=self.matched_case(lambda u:.5)
        result=s.feedback(**kw,mode='matched_actual_outer')
        self.assertTrue(result['promotion_eligible'])
        self.assertLess(result['normalized_port_prediction_error'],1e-10)
        self.assertIsNone(s.check_actual(d,.5,2.07,1.57e-12,kw['policy']))
        self.assertIn('growth',s.check_actual(d,.5,2.07,4e-12,kw['policy']))

    def test_matched_feedback_missing_history_and_frame_mismatch_never_promote(self):
        import copy
        d,kw=self.matched_case()
        missing=copy.deepcopy(d);missing.pop('samples')
        for values in (dict(kw,proposal=missing),dict(kw,actual_reference_origin=6.),
                       dict(kw,actual_reference_origin=None)):
            with self.subTest(values=values.get('actual_reference_origin')):
                result=s.feedback(**values,mode='matched_actual_outer')
                self.assertFalse(result['promotion_eligible'])
                self.assertLessEqual(result['factor'],1.)
                self.assertEqual(result['normalized_port_prediction_error'],
                    result['original_normalized_port_prediction_error'])
        state=s.initial_growth_state();state['clean_streak']=4
        d['growth_policy']=dict(proposed_current_fraction=.04)
        s.growth_observe(state,s.GrowthPolicy(),d,candidate_error=0.,candidate_closure=0.,
            row_ratio=0.,actual_current=kw['current'],recovered=False,
            prediction_error=0.,prediction_valid=False)
        self.assertEqual(state['stage'],0);self.assertEqual(state['clean_streak'],0)
        self.assertEqual(state['multiplier'],.5)

    def test_matched_feedback_rechecks_extended_domain_and_frozen_coefficients(self):
        import copy
        law=lambda u:.01*u/(1-u)
        rs=self.hybrid_records([0,.05,.1],law)
        p=s.Policy(predictor='limited_rational',current_fraction=1.)
        d=s.propose(rs,p,.001,10.)
        # Explicit valid frozen rational model before a future pole crossing.
        d['fit']['model']=s.port_model(d['samples'],'rational3',.001)
        self.assertTrue(d['fit']['model']['valid'])
        kw=dict(proposal=d,inner=.02,current=1e-12,actual_span=1.,policy=p,
                cost=dict(iteration_rows=1,elapsed_s=1.),actual_reference_origin=5.)
        result=s.feedback(**kw,mode='matched_actual_outer')
        self.assertFalse(result['promotion_eligible'])
        self.assertEqual(result['matched_prediction']['reason'],'pole_or_condition')
        d2=copy.deepcopy(d);d2['fit']['model']['a']*=1.01
        r=s.feedback(**dict(kw,proposal=d2,actual_span=.001),mode='matched_actual_outer')
        self.assertFalse(r['promotion_eligible'])
        self.assertEqual(r['matched_prediction']['reason'],'frozen_model_history_mismatch')

    def test_matched_feedback_does_not_bypass_growth_quality(self):
        d,kw=self.matched_case();d['growth_policy']=dict(proposed_current_fraction=.04)
        for changes in (dict(candidate_closure=.02),dict(candidate_error=.01),
                        dict(row_ratio=.00095),dict(recovered=True)):
            state=s.initial_growth_state();state['clean_streak']=4
            metrics=dict(candidate_error=0.,candidate_closure=0.,row_ratio=0.,
                recovered=False,actual_current=kw['current'],prediction_error=0.,prediction_valid=True)
            metrics.update(changes)
            s.growth_observe(state,s.GrowthPolicy(),d,**metrics)
            self.assertEqual(state['stage'],0);self.assertEqual(state['clean_streak'],0)

    def test_growth_reference_change_discards_calibration_and_promotion(self):
        rs=self.hybrid_records([40.,41.,42.],lambda u:38.5)
        policy=s.Policy(predictor='limited_rational');control=s.GrowthPolicy();state=s.initial_growth_state()
        s.growth_propose(rs,policy,control,state,.005,10.)
        state.update(stage=2,clean_streak=4,observations=[dict(span_V=.002,candidate_error=.002)])
        changed=[dict(base={'potential_origin_V':6.},outer_V=36.,point=dict(inner_V=32.5,current_A_per_um=3.5e-12))]
        d=s.growth_propose(changed,policy,control,state,.005,10.)
        self.assertEqual(state['stage'],0);self.assertEqual(state['clean_streak'],0)
        self.assertEqual(state['observations'],[])
        self.assertIsNone(d['growth_policy']['candidate_error_cap_V'])

    def test_unified_acceptance_preserves_positive_branch(self):
        policy=s.Policy()
        d=s.propose(records(),policy,.01,10.)
        for current in (d['parent']['current_A_per_um'],.99*d['parent']['current_A_per_um']):
            self.assertIn('positive branch',s.check_actual(d,d['parent']['inner_V']+.01,
                d['parent']['outer_V']+.01,current,policy))

    def test_growth_control_vertical_branch_and_final_guard_are_separate(self):
        rs=self.hybrid_records([40.,41.,42.],lambda u:38.5)
        policy=s.Policy(predictor='limited_rational',outer_max_V=.01,accepted_inner_max_V=.1)
        control=s.GrowthPolicy();state=s.initial_growth_state()
        d=s.growth_propose(rs,policy,control,state,10.,1e5)
        self.assertAlmostEqual(d['span_V'],.04*3.5)
        self.assertAlmostEqual(d['accepted_outer_max_V'],.25*3.5)
        self.assertAlmostEqual(d['current_allowance_A_per_um'],.25*3.5e-12)
        self.assertIsNone(s.check_actual(d,38.5,42.2,3.7e-12,policy))
        self.assertIn('growth',s.check_actual(d,38.5,42.2,4.5e-12,policy))

    def test_growth_candidate_cap_is_local_and_does_not_relax_guard(self):
        rs=self.hybrid_records([40.,41.,42.],lambda u:38.5)
        policy=s.Policy(predictor='limited_rational',outer_max_V=.01)
        state=s.initial_growth_state();state['observations']=[dict(span_V=.002,candidate_error=.002)]
        d=s.growth_propose(rs,policy,s.GrowthPolicy(),state,1.,10.)
        self.assertAlmostEqual(d['span_V'],.005)
        self.assertAlmostEqual(d['growth_policy']['candidate_error_cap_V'],.005)
        self.assertEqual(d['growth_policy']['actual_current_fraction_guard'],.25)

    def test_growth_promotes_only_clean_confirmed_observations_and_backs_off(self):
        rs=self.hybrid_records([40.,41.,42.],lambda u:38.5)
        policy=s.Policy(predictor='limited_rational');control=s.GrowthPolicy();state=s.initial_growth_state()
        d=s.growth_propose(rs,policy,control,state,.01,10.)
        kw=dict(candidate_error=.002,candidate_closure=.002,row_ratio=.0005,
                actual_current=d['predicted_current_A_per_um'],recovered=False,prediction_error=.01)
        for _ in range(4):s.growth_observe(state,control,d,**kw)
        self.assertEqual(state['stage'],0)
        s.growth_observe(state,control,d,**kw);self.assertEqual(state['stage'],1)
        s.growth_failed(state);self.assertEqual(state['multiplier'],.5)
        self.assertEqual(state['clean_streak'],0)
        self.assertEqual(len(state['observations']),5)  # failed trial is not a fit sample
        s.growth_observe(state,control,d,**dict(kw,row_ratio=.00095))
        self.assertEqual(state['clean_streak'],0);self.assertEqual(state['multiplier'],.25)

    def test_growth_headroom_and_resume_state(self):
        import json
        policy=s.Policy(predictor='limited_rational')
        state=s.initial_growth_state();s.growth_failed(state)
        copied=json.loads(json.dumps(state))
        rs=self.hybrid_records([40.,41.,42.],lambda u:38.5)
        self.assertEqual(s.growth_propose(rs,policy,s.GrowthPolicy(),state,.01,10.),
                         s.growth_propose(rs,policy,s.GrowthPolicy(),copied,.01,10.))
        with self.assertRaises(ValueError):s.GrowthPolicy().validate(s.Policy(predictor='limited_rational',current_fraction=.20))

    def hybrid_records(self, us, law):
        return [dict(base={'potential_origin_V':5.},outer_V=u-5.,
                     point=dict(inner_V=law(u)-5.,current_A_per_um=(u-law(u))/1e12)) for u in us]

    def test_hybrid_preserves_linear_loadline_and_independent_limits(self):
        rs=self.hybrid_records([1,1.1,1.25,1.4,1.6,1.9],lambda u:.4*u)
        p=s.Policy(predictor='limited_rational',outer_max_V=10,accepted_inner_max_V=.1,current_fraction=1.)
        d=s.propose(rs,p,2,10)
        self.assertEqual(d['fit']['model']['kind'],'limited_endpoint')
        self.assertAlmostEqual(d['span_V'],.25)
        dv=d['predicted_inner_V']-d['parent']['inner_V']
        di=d['predicted_current_A_per_um']-d['parent']['current_A_per_um']
        self.assertAlmostEqual(dv+1e12*di,d['span_V'])

    def test_rational_selected_by_completed_holdouts_and_predicts_unseen_point(self):
        law=lambda u:.8*u/(1+.1*u)
        rs=self.hybrid_records([1,1.2,1.4,1.6,1.8,2],law)
        p=s.Policy(predictor='limited_rational',current_fraction=1.)
        d=s.propose(rs,p,.1,10)
        self.assertEqual(d['fit']['model']['kind'],'rational3')
        self.assertEqual(len(d['fit']['selection']['holdouts']),3)
        self.assertAlmostEqual(d['predicted_inner_V'],law(2+d['span_V']),places=12)
        early=s.propose(rs[:5],p,.1,10)
        self.assertEqual(early['fit']['model']['kind'],'limited_endpoint')

    def test_rational_pole_or_derivative_violation_is_ineligible(self):
        rs=self.hybrid_records([0,.1,.2],lambda u:.1*u/(1-u))
        points=s.port_samples(rs)
        m=s.port_model(points,'rational3',1.)
        self.assertFalse(m['valid'])
        self.assertIn(m['reason'],('pole_or_condition','derivative_shape'))

    def test_reversal_is_not_hidden_and_frames_reset_holdouts(self):
        rs=self.hybrid_records([1,2,3,4,5,6],lambda u:.4*u)
        rs[-1]['point']['inner_V']=rs[-2]['point']['inner_V']-.0000001
        self.assertEqual(len(s.port_samples(rs)),1)
        d=s.propose(rs,s.Policy(predictor='limited_rational'),.2,10)
        self.assertFalse(d['fit']['trusted'])
        rs=self.hybrid_records([1,2,3,4,5,6],lambda u:.4*u)
        rs[-1]['base']={'potential_origin_V':15.}
        self.assertEqual(len(s.port_samples(rs)),1)

    def test_internal_plateau_can_be_predicted_with_current_cap(self):
        rs=self.hybrid_records([1,2,3,4,5,6],lambda u:.5)
        p=s.Policy(predictor='limited_rational',outer_max_V=100,current_fraction=.01)
        d=s.propose(rs,p,50,100)
        self.assertEqual(d['fit']['model']['a'],0.)
        self.assertAlmostEqual(d['span_V'],.055)
        self.assertAlmostEqual(d['predicted_inner_V'],.5)

    def test_hybrid_current_bound_and_actual_rejection_remain_independent(self):
        rs=self.hybrid_records([1,1.2,1.4,1.6,1.8,2],lambda u:.8*u/(1+.1*u))
        p=s.Policy(predictor='limited_rational',outer_max_V=100,current_fraction=.01)
        d=s.propose(rs,p,.5,100)
        change=d['predicted_current_A_per_um']-d['parent']['current_A_per_um']
        self.assertAlmostEqual(change/d['current_allowance_A_per_um'],1.,places=10)
        self.assertIn('current',s.check_actual(d,d['predicted_inner_V'],d['parent']['outer_V']+d['span_V'],
            d['parent']['current_A_per_um']+2*change,p))

    def test_nonuniform_linear_device_exact_port_prediction(self):
        p=s.Policy(outer_max_V=10.,accepted_inner_max_V=.2,current_fraction=1.)
        d=s.propose(records(),p,5.,100.)
        self.assertAlmostEqual(d['fit']['slope_A_per_um_V'],2e-12)
        self.assertAlmostEqual(d['span_V'],.6)
        self.assertAlmostEqual(d['predicted_inner_V'],1.4)
        self.assertAlmostEqual(d['predicted_current_A_per_um'],3.8e-12)
        self.assertIsNone(s.check_actual(d,1.4,d['parent']['outer_V']+.6,3.8e-12,p))

    def test_current_limit_dominates_large_resistor_gain(self):
        p=s.Policy(outer_max_V=1e6,accepted_inner_max_V=1.,current_fraction=.01)
        d=s.propose(records(),p,1e6,1e6)
        self.assertEqual(d['limiting'],'current_V')
        self.assertAlmostEqual(d['predicted_current_A_per_um']-3.4e-12,.01*3.4e-12)

    def test_negative_and_inconsistent_slopes_do_not_amplify(self):
        rs=records(slope=-1e-14)
        d=s.propose(rs,s.Policy(outer_max_V=100.),10.,100.)
        self.assertFalse(d['fit']['trusted'])
        self.assertLess(d['fit']['slope_A_per_um_V'],0)
        self.assertLessEqual(d['span_V'],.2)
        rs=records();rs[-1]['point']['current_A_per_um']*=10
        self.assertFalse(s.estimate(s.samples(rs))['trusted'])

    def test_duplicate_frame_records_and_small_spacing(self):
        rs=records(origin=5.)
        rs.append(rs[-1].copy())
        self.assertEqual(len(s.samples(rs)),3)
        self.assertAlmostEqual(s.samples(rs)[-1]['inner_V'],1.2)
        rs.append(records(xs=(1.3,),origin=15.)[0])
        self.assertEqual(len(s.samples(rs)),1)
        self.assertFalse(s.propose(rs,s.Policy(),.1,10.)['fit']['trusted'])

    def test_goal_and_independent_caps(self):
        p=s.Policy(outer_max_V=.5,accepted_inner_max_V=.02,trial_inner_max_V=.007)
        d=s.propose(records(),p,1.,.01)
        self.assertEqual(d['span_V'],.01)
        self.assertEqual(d['limiting'],'remaining_V')
        self.assertIn('internal',s.check_actual(d,1.23,d['parent']['outer_V']+.01,3.4e-12,p))
        self.assertIn('external',s.check_actual(d,1.21,d['parent']['outer_V']+.6,3.4e-12,p))
        self.assertIn('current',s.check_actual(d,1.21,d['parent']['outer_V']+.01,4.4e-12,p))

    def test_feedback_shrinks_for_error_cost_and_recovery(self):
        p=s.Policy();d=s.propose(records(slope=1e-16),p,.1,10.)
        kwargs=dict(proposal=d,inner=d['predicted_inner_V'],current=d['predicted_current_A_per_um'],
                    actual_span=.1,policy=p,cost={'iteration_rows':4,'elapsed_s':10.})
        self.assertAlmostEqual(s.feedback(**kwargs)['next_span_V'],.14)
        with self.assertRaises(ValueError):
            s.feedback(**dict(kwargs,actual_span=2.27e-7),minimum_actual_span=1e-5)
        self.assertAlmostEqual(s.feedback(**kwargs,minimum_actual_span=1e-5)['next_span_V'],.14)
        self.assertEqual(s.feedback(**kwargs,recovered=True)['next_span_V'],.05)
        kwargs['cost']={'iteration_rows':80,'elapsed_s':10.}
        self.assertEqual(s.feedback(**kwargs)['next_span_V'],.05)
        kwargs['inner']+=.1
        self.assertEqual(s.feedback(**kwargs)['next_span_V'],.05)

    def test_cost_includes_failed_and_correction_histories(self):
        with tempfile.TemporaryDirectory() as tmp:
            for folder in ('failed','corrector','accepted'):
                p=Path(tmp)/folder;p.mkdir()
                (p/'iv_newton_history.csv').write_text('iteration\n0\n1\n2\n')
            c=s.measure_cost(tmp,3.)
            self.assertEqual(c,dict(iteration_rows=6,history_files=3,elapsed_s=3.))

    def test_bad_policy_rejected(self):
        for kw in ({'outer_max_V':float('nan')},{'trial_inner_max_V':0},{'growth':3}):
            with self.assertRaises(ValueError):s.Policy(**kw)


if __name__=='__main__':unittest.main()
