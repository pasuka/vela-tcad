"""Fresh equal-thread repeats cannot inherit historical two-thread qualification."""
import copy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import run_templates_ldmos_single_thread_repeats as study
import run_templates_ldmos_full_curve_timing as shared


class SingleThreadRepeatTests(unittest.TestCase):
    def test_qualification_precedes_five_backend_three_rounds(self):
        rows=study.schedule()
        self.assertEqual(len(rows),30)
        self.assertEqual(len({r['key'] for r in rows}),30)
        self.assertEqual([(r['config'],r['gate']) for r in rows[:2]],[('T1',4),('T1',8)])
        self.assertEqual({(c['threads'],c['reuse']) for c in study.CONFIGURATIONS.values()},{(1,True)})
        for config in study.IDS:
            for gate in (4,8):self.assertEqual([r['round'] for r in rows if r['config']==config and r['gate']==gate],[0,1,2])
        for repeat in (1,2):
            four=[r['config'] for r in rows if r['round']==repeat and r['gate']==4]
            eight=[r['config'] for r in rows if r['round']==repeat and r['gate']==8]
            self.assertEqual(four,eight[::-1])

    def test_resume_rejects_thread_mode_and_schedule_changes(self):
        report=dict(status='paused',configurations=copy.deepcopy(study.CONFIGURATIONS),schedule=study.schedule(),control_config='U1')
        study.validate_resume(report)
        for mutate in (lambda r:r['configurations']['T1'].update(threads=2),
                       lambda r:r['configurations']['U1'].update(reuse=False),
                       lambda r:r.update(control_config='U0'),lambda r:r.update(status='running')):
            bad=copy.deepcopy(report);mutate(bad)
            with self.assertRaises(ValueError):study.validate_resume(bad)

    def test_failing_single_thread_qualification_stops_before_repeats(self):
        row=dict(study.schedule()[0],status='pass',name='qualification')
        report=dict(cases=[row],qualification_controls={'4':'old'},qualification_comparisons={})
        with patch.object(study,'compare',return_value={'equivalent':False}):
            with self.assertRaises(ValueError):study.check_qualification(Path('unused'),report)

    def test_reuse_umfpack_is_explicit_state_control(self):
        rows=[dict(r,status='pass',name=r['key']) for r in study.schedule() if r['round']==0]
        report=dict(cases=rows,control_config='U1',comparisons={},repeat_comparisons={},joint_qualifications={})
        with patch.object(shared,'compare',return_value={'equivalent':True}) as compare, \
             patch.object(shared,'read',return_value={'references':{'4':'ref4','8':'ref8'}}), \
             patch.object(shared,'analyze',return_value={k:{'status':'pass'} for k in ('engineering','final')}) as joint:
            shared.check_ready(Path('unused'),report)
            self.assertEqual(compare.call_count,8)
            self.assertEqual(joint.call_count,5)
            self.assertTrue(all(v['control'].endswith('_U1') for v in report['comparisons'].values()))


if __name__=='__main__':unittest.main()
