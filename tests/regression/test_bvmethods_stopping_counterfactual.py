import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SCRIPT=Path(__file__).resolve().parents[2]/'scripts/diagnostics/bvmethods_stopping_counterfactual.py'
spec=importlib.util.spec_from_file_location('stopping_counterfactual',SCRIPT)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class StoppingReviewTests(unittest.TestCase):
    def test_OR_uses_applied_report_and_keeps_the_load_gate(self):
        row={'device_residual_inf':'1.01157e-8','load_residual_V':'1e-9','update_valid':'1','update_rms':'.25'}
        self.assertTrue(module.evaluate(row)['or_scalar_eligible'])
        self.assertFalse(module.evaluate(row)['and_scalar_eligible'])
        row['load_residual_V']='.01'
        self.assertFalse(module.evaluate(row)['or_scalar_eligible'])
        row['load_residual_V']='1e-9';row['update_valid']='0';row['candidate_update_rms']='.1'
        self.assertFalse(module.evaluate(row)['or_scalar_eligible'])

    def test_outgoing_rejection_does_not_erase_a_prior_applied_report(self):
        row={'device_residual':'1.01e-8','load_V':'0','update_valid':'1','rms':'.69','rejected':'1'}
        e=module.evaluate(row);self.assertTrue(e['or_scalar_eligible']);self.assertTrue(e['outgoing_candidate_rejected'])
        row['rms']='nan';self.assertFalse(module.evaluate(row)['or_scalar_eligible'])

    def test_residual_branch_can_pass_with_large_or_unavailable_update(self):
        row={'device_residual':'4.08e-9','load_V':'4.37e-9','update_valid':'1','rms':'1083'}
        result=module.evaluate(row)
        self.assertTrue(result['residual_branch'])
        self.assertFalse(result['update_branch'])
        self.assertTrue(result['or_scalar_eligible'])
        self.assertFalse(result['and_scalar_eligible'])
        row.update(update_valid='0',rms='nan')
        self.assertTrue(module.evaluate(row)['or_scalar_eligible'])
        row['load_V']='.01'
        self.assertFalse(module.evaluate(row)['or_scalar_eligible'])

    def test_configuration_changes_only_the_requested_combination(self):
        cfg={'sweep':{'external_circuit':{'solver':'coupled_newton','load_line_convergence_mode':'enforce','coupled_update_convergence':{'mode':'enforce_joint','digits':5}}}}
        old=json.dumps(cfg,sort_keys=True);changed=module.configure_mode(cfg)
        self.assertEqual(json.dumps(cfg,sort_keys=True),old)
        self.assertEqual(changed['sweep']['external_circuit']['coupled_update_convergence']['mode'],'enforce')
        changed['sweep']['external_circuit']['coupled_update_convergence']['mode']='enforce_joint'
        self.assertEqual(changed,cfg)

    def test_repeated_iteration_is_one_solve_but_a_reset_starts_a_new_solve(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'updates.csv'
            p.write_text('outer_V,step_V,retry,iteration,device_residual,load_V,update_valid,rms\n1,.01,0,0,1,0,0,nan\n1,.01,0,0,1,0,0,nan\n1,.01,0,1,1.01e-8,0,1,.5\n1,.01,0,0,1,0,0,nan\n1,.01,0,1,1e-10,0,1,.5\n')
            a=module.analyze(p);self.assertEqual(a['attempts'],2);self.assertEqual(a['OR_stops_earlier_attempts'],1)

if __name__=='__main__':unittest.main()
