"""Guard, isolation and acceptance tests for the linked D5 entry point."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import run_templates_ldmos_linked_d5 as linked


class LinkedD5Test(unittest.TestCase):
    def test_auger_density_profiles_change_only_the_explicit_source_model(self):
        folder=ROOT/'reference_tcad/templates_ldmos_sentaurus2022/profiles'
        for physics in ('d4','d5'):
            baseline=linked.read(folder/f'linked_{physics}_fermi_accurate_config.json')
            bundle=linked.read(folder/f'linked_{physics}_auger_density_inputs.json')
            profile=linked.read(ROOT/bundle['template'])
            enhancement=profile['solver'].pop('auger_density_dependence')
            self.assertTrue(enhancement['enabled'])
            self.assertEqual(enhancement['electron'],{'enhancement':3.46667,'reference_density_m3':1e18})
            self.assertEqual(enhancement['hole'],{'enhancement':8.25688,'reference_density_m3':1e18})
            self.assertEqual(profile['solver'],baseline['solver'])
            self.assertEqual(profile['contacts'],baseline['contacts'])
            self.assertEqual(profile['materials_file'],baseline['materials_file'])
            self.assertEqual(profile['discretization'],baseline['discretization'])
            self.assertEqual(profile['scaling']['mode'],'unit_scaling')

    def test_accurate_fermi_profiles_preserve_original_curve_gates_and_physics(self):
        for physics,original in [('d4','linked_d4'),('d5','linked_d5_auger_units')]:
            folder=ROOT/'reference_tcad/templates_ldmos_sentaurus2022/profiles'
            baseline=linked.read(folder/f'{original}_config.json')
            bundle=linked.read(folder/f'linked_{physics}_fermi_accurate_inputs.json')
            profile=linked.read(ROOT/bundle['template'])
            self.assertEqual(profile['solver'],baseline['solver'])
            self.assertEqual(profile['discretization'],baseline['discretization'])
            self.assertEqual(profile['contacts'],baseline['contacts'])
            self.assertEqual(linked.digest(ROOT/bundle['template']),bundle['files'][bundle['template']])
            material=profile['materials_file'].replace('@workspace/', '')
            self.assertEqual(linked.digest(ROOT/material),bundle['files'][material])
            self.assertNotEqual(profile['materials_file'],baseline['materials_file'])

    def test_d4_qualified_bundle_keeps_its_physics_and_frame_policy(self):
        bundle=linked.read(ROOT/'reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d4_inputs.json')
        profile=linked.read(ROOT/bundle['template'])
        self.assertEqual(bundle['schema'],'vela.templates_ldmos.linked_d4_inputs.v1')
        self.assertEqual(linked.digest(ROOT/bundle['template']),bundle['files'][bundle['template']])
        self.assertEqual(profile['solver']['mobility']['model'],'ialmob')
        self.assertEqual(profile['solver']['block_absolute_convergence'],linked.GATES)
        self.assertEqual(linked.frame_offset_for(bundle,4),28.)
        self.assertEqual(linked.frame_offset_for(bundle,8),4.)
        self.assertTrue(all('D4-classical' in p for p in bundle['references'].values()))

    def test_qualified_frame_policy_is_explicit_and_overridable(self):
        self.assertEqual(linked.frame_offset_for({},8),28.)
        bundle={'frame_offsets_V':{'4':28.,'8':16.}}
        self.assertEqual(linked.frame_offset_for(bundle,8),16.)
        self.assertEqual(linked.frame_offset_for(bundle,4),28.)
        self.assertEqual(linked.frame_offset_for(bundle,8,8.),8.)
        for value in (0.,float('nan'),float('inf'),-1.,41.):
            with self.assertRaises(ValueError):linked.frame_offset_for(bundle,8,value)

    def test_resume_requires_inactive_matching_exact_prefix(self):
        from copy import deepcopy
        plan=dict(runner_sha256='runner',gate_V=8,max_step_V=.2,
                  original_blocks=deepcopy(linked.GATES),execution_mode='subprocess')
        ledger=dict(status='stopped_after_failure',active_child=None,
                    exact_points=[{'bias_V':0.},{'bias_V':1.}],accepted_bias_V=1.,
                    frame_V=0.,runs=[{'status':'completed'}])
        check=lambda p,l:linked.validate_resume_contract(p,l,[0.,1.,2.],'runner',8,.2,16.,False)
        check(plan,ledger)
        for key,value in [('active_child','child'),('status','running'),('frame_V',28.),
                          ('accepted_bias_V',2.),('exact_points',[{'bias_V':0.},{'bias_V':.5}]),
                          ('runs',[{'status':'running'}])]:
            changed=deepcopy(ledger);changed[key]=value
            with self.assertRaises(ValueError,msg=key):check(plan,changed)
        for key,value in [('runner_sha256','other'),('gate_V',4),('max_step_V',.4),
                          ('execution_mode','dc_worker'),('original_blocks',{})]:
            changed=deepcopy(plan);changed[key]=value
            with self.assertRaises(ValueError,msg=key):check(changed,ledger)

    def test_d4_requires_explicit_profile_and_preserves_original_gates(self):
        from copy import deepcopy
        config=deepcopy(self.base)
        config['solver']['mobility'].update(model='ialmob',ialmob={'high_field':True,'reference_density_m3':1e18})
        original=linked.PHYSICS_PROFILE
        self.addCleanup(setattr,linked,'PHYSICS_PROFILE',original)
        linked.PHYSICS_PROFILE='D5'
        with self.assertRaises(ValueError):linked.assert_linked_physics_contract(config)
        linked.PHYSICS_PROFILE='D4'
        linked.assert_linked_physics_contract(config)
        config['solver']['mobility']['ialmob']['geometry_file']='@workspace/geometry.json'
        expanded=linked.prepare_config(config,ROOT/'out',ROOT/'seed',1.,.1,8)
        self.assertEqual(expanded['solver']['mobility']['ialmob']['geometry_file'],
                         str(linked.ROOT)+'/geometry.json')
        config['solver']['mobility']['ialmob']['reference_density_m3']=1e19
        with self.assertRaises(ValueError):linked.assert_linked_physics_contract(config)
        config['solver']['mobility']['ialmob']['reference_density_m3']=1e18
        config['solver']['block_absolute_convergence']['hole_residual_ceiling']*=10
        with self.assertRaises(ValueError):
            linked.prepare_config(config,ROOT/'out',ROOT/'seed',1.,.1,8)

    def test_default_bundle_uses_corrected_auger_units(self):
        bundle=linked.read(linked.DEFAULT_BUNDLE)
        profile=linked.read(ROOT/bundle['template'])
        self.assertEqual(linked.digest(ROOT/bundle['template']),bundle['files'][bundle['template']])
        # ConfigParser consumes the explicit coefficients in cm^6/s.
        self.assertEqual(profile['scaling']['mode'],'unit_scaling')
        self.assertEqual(profile['solver']['auger_cn_m6_per_s'],2.9e-31)
        self.assertEqual(profile['solver']['auger_cp_m6_per_s'],1.028e-31)

    def setUp(self):
        self.base = json.loads((ROOT / 'reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_config.json').read_text())

    def test_config_preserves_physics_and_physical_contacts(self):
        for gate in (4, 8):
            for frame, target in ((0., 1.), (28., 40.)):
                cfg = linked.prepare_config(self.base, ROOT/'build-release/isolated', ROOT/'seed.csv', target, .2, gate, frame)
                contacts = {c['name']: c['bias'] for c in cfg['contacts']}
                self.assertEqual(contacts['gate']-contacts['source'], gate)
                self.assertEqual(contacts['drain']-contacts['source'], target)
                self.assertEqual(cfg['solver']['mobility'], self.base['solver']['mobility'])
                self.assertEqual(cfg['solver']['block_absolute_convergence'], linked.GATES)
                self.assertEqual(cfg['sweep']['bias_points'], [target-frame])
                self.assertEqual(cfg['solver']['quasi_fermi_update_limit_V'], .2)
                self.assertNotIn('@', json.dumps(cfg))

    def test_reject_changed_acceptance_gate(self):
        self.base['solver']['block_absolute_convergence']['electron_residual_ceiling'] *= 10
        with self.assertRaises(ValueError):
            linked.prepare_config(self.base, ROOT/'out', ROOT/'seed', 1., .1, 4)

    def test_preflight_is_read_only_and_checks_dependencies(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            data = root/'input'; data.write_bytes(b'fixture')
            runner = root/'runner'; runner.write_bytes(b'not an executable')
            bundle = root/'bundle.json'; manifest = root/'manifest.json'
            bundle.write_text(json.dumps({'files': {'input': linked.digest(data)}}))
            manifest.write_text(json.dumps({'runner':str(runner), 'runner_sha256':linked.digest(runner), 'frozen_sources':{}}))
            before = set(root.iterdir())
            linked.preflight(bundle, root, manifest, 4)
            self.assertEqual(before, set(root.iterdir()))
            data.write_bytes(b'changed')
            with self.assertRaises(ValueError):
                linked.preflight(bundle, root, manifest, 4)

    def test_predictor_guard_protects_frame_clipped_step_and_retry(self):
        sweep = object.__new__(linked.Sweep)
        sweep.ledger = {'frame_V':0., 'rollbacks':[], 'transfers':[]}
        guard = lambda: sweep.prediction_guard('direct', .1, .3)
        self.assertEqual(guard(), 'no_history')
        last = dict(parent_V=0., bias_V=.1, accepted_step_V=.1, proposed_step_V=.1, frame_V=0.)
        sweep.ledger['transfers'] = [last]
        self.assertIsNone(guard())
        sweep.ledger['frame_V'] = 28.
        self.assertEqual(guard(), 'frame_changed')
        sweep.ledger['frame_V'] = 0.
        last['proposed_step_V'] = .2
        self.assertEqual(guard(), 'previous_step_clipped')
        last['proposed_step_V'] = .1
        self.assertEqual(sweep.prediction_guard('direct', .1, .31), 'extrapolation_ratio_exceeded')
        sweep.ledger['rollbacks'] = [{'parent_V': .1}]
        self.assertEqual(guard(), 'retry')
        self.assertEqual(sweep.prediction_guard('reclose', .1, .1), 'non_direct')

    def test_strict_blocks_reject_nan_and_over_limit(self):
        row = dict(final_psi_residual_norm=5e-8, final_electron_continuity_residual_norm=1e-11, final_hole_continuity_residual_norm=3e-10)
        self.assertTrue(linked.blocks_pass(row))
        for value in (float('nan'), float('inf'), 1.001e-11):
            row['final_electron_continuity_residual_norm'] = value
            self.assertFalse(linked.blocks_pass(row))



if __name__ == '__main__':
    unittest.main()
