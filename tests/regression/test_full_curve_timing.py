"""Full-curve schedule, original qualification and immutable resume contracts."""
import copy
import sys
from pathlib import Path
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from run_templates_ldmos_full_curve_timing import CONFIGS, schedule, require_joint_pass, validate_resume
import run_templates_ldmos_full_curve_timing as timing


class FullCurveTimingTests(unittest.TestCase):
    def test_all_backends_use_one_solver_thread(self):
        self.assertEqual({c['threads'] for c in CONFIGS.values()},{1})
        # Historical two-thread data cannot be resumed as single-thread repeats.
        old=dict(status='paused',configurations=copy.deepcopy(CONFIGS),schedule=schedule())
        old['configurations']['T1']['threads']=2
        with self.assertRaises(ValueError):validate_resume(old)

    def test_plan_counts_and_unique_cases(self):
        rows=schedule()
        self.assertEqual(len(rows),48)
        self.assertEqual(len({r['key'] for r in rows}),48)
        self.assertEqual(sum(r['points'] for r in rows),1350)
        self.assertEqual(sum(r['profile']=='D5' for r in rows),24)
        for r in rows:
            self.assertEqual(r['points'],8 if r['stage']=='P3short' else 31)

    def test_balanced_priority_order_and_explicit_control(self):
        rows=schedule()
        for repeat in range(3):
            groups={g:[r['config'] for r in rows if r['profile']=='D5' and r['round']==repeat and r['gate']==g and r['config'] in ('U0','U1','T1')] for g in (4,8)}
            self.assertEqual(groups[4],['U0','U1','T1'][repeat:]+['U0','U1','T1'][:repeat])
            self.assertEqual(groups[8],groups[4][::-1])
        self.assertEqual(CONFIGS['U0']|{'reuse':True},CONFIGS['U1'])

    def test_final_joint_gate_cannot_be_replaced_by_engineering(self):
        good={k:dict(status='pass') for k in ('engineering','final')}
        require_joint_pass(good)
        for level in good:
            bad=copy.deepcopy(good);bad[level]['status']='fail'
            with self.assertRaises(ValueError):require_joint_pass(bad)
        with self.assertRaises(ValueError):require_joint_pass({'engineering':{'status':'pass'}})

    def test_resume_rejects_active_failed_and_changed_contracts(self):
        base=dict(status='paused',configurations=copy.deepcopy(CONFIGS),schedule=schedule())
        validate_resume(base)
        for status in ('running','failed','pass'):
            with self.assertRaises(ValueError):validate_resume(base|{'status':status})
        with self.assertRaises(ValueError):validate_resume(base|{'active':'case'})
        bad=copy.deepcopy(base);bad['configurations']['T1']['threads']=4
        with self.assertRaises(ValueError):validate_resume(bad)
        bad=copy.deepcopy(base);bad['schedule'][0]['points']=8
        with self.assertRaises(ValueError):validate_resume(bad)

    def test_late_control_triggers_comparison_and_full_joint_gate(self):
        rows=[dict(r,status='pass',name=r['key']) for r in schedule()[:6]]
        report=dict(cases=rows,comparisons={},repeat_comparisons={},joint_qualifications={})
        with patch.object(timing,'compare',return_value={'equivalent':True}) as compare, \
             patch.object(timing,'read',return_value={'references':{'4':'r4.csv','8':'r8.csv'}}), \
             patch.object(timing,'analyze',return_value={k:{'status':'pass'} for k in ('engineering','final')}) as joint:
            timing.check_ready(Path('unused'),report)
            self.assertEqual(compare.call_count,4)
            self.assertEqual(joint.call_count,3)
            timing.check_ready(Path('unused'),report)
            self.assertEqual(compare.call_count,4)
            self.assertEqual(joint.call_count,3)

    def test_failed_cross_backend_state_blocks_qualification(self):
        rows=[dict(r,status='pass',name=r['key']) for r in schedule()[:2]]
        report=dict(cases=rows,comparisons={},repeat_comparisons={},joint_qualifications={})
        with patch.object(timing,'compare',return_value={'equivalent':False}), \
             patch.object(timing,'analyze') as joint:
            with self.assertRaises(ValueError):timing.check_ready(Path('unused'),report)
            joint.assert_not_called()


if __name__=='__main__':unittest.main()
