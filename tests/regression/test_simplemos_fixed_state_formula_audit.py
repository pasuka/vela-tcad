"""Guard units, conservative cancellation, and zero-current field masking."""
import sys
import math
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
import run_simplemos_fixed_state_formula_audit as audit


class FixedStateAuditTests(unittest.TestCase):
    def test_contact_flux_telescopes_and_has_physical_depth_units(self):
        mesh={'contacts':[{'name':'left','node_ids':[0]}, {'name':'right','node_ids':[2]}]}
        edges=[{'node0':0,'node1':1},{'node0':1,'node1':2}]
        flux=np.array([3.,3.]);zero=np.zeros(2)
        values=audit.sum_contacts(edges,mesh,flux,zero)
        self.assertAlmostEqual(values['right']/(audit.Q*1e-6),3.)
        self.assertEqual(values['left']+values['right'],0.)

    def test_zero_candidate_is_not_removed_from_reference_mask(self):
        result=audit.vector_metric([[1.,0.],[1.,0.]],[[1.,0.],[0.,0.]])
        self.assertEqual(result['nodes'],2)
        self.assertEqual(result['candidate_zero_nodes'],1)
        self.assertGreater(result['p95_dex'],250.)

    def test_bernoulli_equilibrium_identity(self):
        eta=np.array([-10.,-1e-7,0.,1e-7,10.],dtype=np.longdouble)
        np.testing.assert_allclose(audit.bernoulli(-eta),np.exp(eta)*audit.bernoulli(eta),rtol=2e-15)

    def test_density_flux_direction_and_equilibrium(self):
        eta=np.array([0.,1.],dtype=np.longdouble)
        got=audit.electron_density_flux(np.array([2.,1.]),np.array([1.,math.e]),eta,1.)
        self.assertEqual(got[0],1.)
        self.assertLess(abs(got[1]),1e-15)


if __name__=='__main__': unittest.main()
