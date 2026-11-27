"""Gauge invariants of an initial guess, not device qualification."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from slot_fermi_poly_predictor import shift_fields,monitor_bias_change


class FermiPredictorTest(unittest.TestCase):
    def test_advisory_warning_does_not_raise_or_suppress_repeated_steps(self):
        with self.assertNoLogs('slot_fermi_poly_predictor',level='WARNING'):
            self.assertFalse(monitor_bias_change(.5,.6)['warning_threshold_exceeded'])
        with self.assertLogs('slot_fermi_poly_predictor',level='WARNING') as logs:
            for parent,target in ((.5,.75),(.75,1.),(.75,.5)):
                result=monitor_bias_change(parent,target)
                self.assertTrue(result['warning_threshold_exceeded'] and result['warning_only'])
        self.assertEqual(len(logs.output),3)
        self.assertTrue(all('continuing unqualified prediction' in x for x in logs.output))
        with self.assertRaises(ValueError):monitor_bias_change(0.,float('nan'))

    def test_occupations_and_unselected_nodes_preserved(self):
        f = dict(psi=np.array([.42, .31, -.2]), phin=np.array([0., .1, .03]),
                 phip=np.array([0., .1, .04]), electrons_m3=np.array([3e26,2e26,1e12]),
                 holes_m3=np.array([1e3,2e3,1e20]))
        for carrier, qf in [('electron','phin'),('hole','phip')]:
            f[carrier+'_qf_reference_V']=f[qf].copy()
            f[carrier+'_qf_increment_V']=np.zeros(3)
        for delta in (.1,-.05):
            p=shift_fields(f,{0:delta,1:delta})
            for q in ('phin','phip'):
                np.testing.assert_allclose(p['psi']-p[q],f['psi']-f[q],atol=1e-16,rtol=0)
            for k in ('electrons_m3','holes_m3'):
                self.assertTrue(np.array_equal(f[k],p[k]))
            for k in f:
                self.assertEqual(f[k][2],p[k][2])
            self.assertEqual(f['psi'][0],.42)

    def test_unsupported_coordinates_and_nonfinite_rejected(self):
        with self.assertRaises(ValueError): shift_fields({'packed_psi':np.zeros(1)}, {})
        with self.assertRaises(ValueError): shift_fields({}, {0:float('nan')})


if __name__=='__main__': unittest.main()
