import math
import unittest
from scripts.slot_population_acceptance import accept, readback_audit


class AcceptanceTests(unittest.TestCase):
    def decision(self, r, u, **changes):
        guards=dict(closure=True,boundary=True,linear=True,mapping=True,positive=True)
        guards.update(changes)
        return accept(r,u,**guards)

    def test_update_branch_survives_rhs_readback_failure(self):
        for r,u in [(1.7011140352543726e-10,.5317051528707994),
                    (2.55809161842238e-10,.9221515641846297)]:
            self.assertTrue(self.decision(r,u))
            audit=readback_audit(r,r,False)
            self.assertTrue(audit['consistent']);self.assertFalse(audit['rhs_pass'])

    def test_rhs_branch_and_no_branch(self):
        self.assertTrue(self.decision(6.69e-11,4.09))
        self.assertFalse(self.decision(5.61e-10,45.36))
        self.assertFalse(self.decision(1.1e-10,1.))

    def test_guards_remain_mandatory(self):
        for guard in ['closure','boundary','linear','mapping','positive']:
            self.assertFalse(self.decision(1e-12,.1,**{guard:False}))

    def test_invalid_metrics_and_changed_state(self):
        for x in [math.nan,math.inf,-1.]:
            self.assertFalse(self.decision(x,.1))
            self.assertFalse(self.decision(1e-12,x))
            self.assertFalse(readback_audit(1e-10,x,True)['consistent'])
        self.assertFalse(readback_audit(1e-10,2e-10,True)['consistent'])

    def test_primary_only_reports_closure_and_boundary(self):
        for r,u in [(1e-10,20.),(1e-7,.9)]:
            self.assertTrue(self.decision(r,u,closure=False,boundary=False,policy="primary_only"))
        self.assertFalse(self.decision(1.01e-10,1.,policy="primary_only"))

    def test_primary_only_preserves_numerical_guards(self):
        for guard in ['linear','mapping','positive']:
            self.assertFalse(self.decision(1e-12,.1,policy="primary_only",**{guard:False}))
        self.assertFalse(self.decision(math.nan,.1,policy="primary_only"))
        self.assertFalse(readback_audit(1e-10,2e-10,False)['consistent'])
        with self.assertRaises(ValueError):self.decision(1e-12,.1,policy="unknown")


if __name__=='__main__':unittest.main()
