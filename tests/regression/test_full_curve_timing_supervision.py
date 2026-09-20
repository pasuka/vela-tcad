"""Automatic advancement preserves user stops and shortlist review boundaries."""
import copy
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from supervise_templates_ldmos_full_curve_timing import review_reasons
from run_templates_ldmos_full_curve_timing import schedule


class SupervisionTests(unittest.TestCase):
    def report(self):
        return dict(status='paused',stop_reason='stage_boundary',joint_qualifications={},
            cases=[dict(x,status='pass',wall_seconds=100 if x['config']=='U1' else 120,
                        audit=dict(rollbacks=0)) for x in schedule() if x['stage']=='P1'])

    def test_qualified_fixed_shortlist_can_advance(self):
        self.assertEqual(review_reasons(self.report(),'P1'),[])

    def test_user_pause_failure_and_partial_stage_do_not_advance(self):
        for status,reason in [('paused','user_interrupt'),('failed','stage_boundary'),('running','stage_boundary')]:
            report=self.report();report.update(status=status,stop_reason=reason)
            self.assertTrue(review_reasons(report,'P1'))
        report=self.report();report['cases'].pop()
        self.assertTrue(review_reasons(report,'P1'))

    def test_faster_alternative_requires_review_instead_of_silent_exclusion(self):
        report=self.report();next(x for x in report['cases'] if x['config']=='M1')['wall_seconds']=90
        self.assertTrue(review_reasons(report,'P1'))
        report=self.report();report['joint_qualifications']['x']={'engineering':{'status':'pass'},'final':{'status':'fail'}}
        self.assertTrue(review_reasons(report,'P1'))


if __name__=='__main__':unittest.main()
