import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from slot_idvd_secant_predictor import extrapolate_fields

class SecantTests(unittest.TestCase):
    def fields(self,inc,density):
        f={k:np.array([0.1,0.2]) for k in ('psi','phin','phip')}
        for c in ('electron','hole'):
            f[c+'_qf_reference_V']=np.array([5.,5.])
            f[c+'_qf_increment_V']=np.array([inc,inc])
        f['electrons_m3']=np.array([density,0.]);f['holes_m3']=np.array([density,0.])
        return f
    def test_sub_ulp_qf_increment_and_positive_density(self):
        old=self.fields(1e-30,2.);parent=self.fields(2e-30,4.)
        out=extrapolate_fields(old,parent,1.)
        self.assertAlmostEqual(out['electron_qf_increment_V'][0]/1e-30,3.)
        self.assertAlmostEqual(out['electrons_m3'][0],8.)
        self.assertEqual(out['electrons_m3'][1],0.)
        self.assertTrue(np.array_equal(parent['electrons_m3'],[4.,0.]))
    def test_invalid_history_ratio_and_packed_rejected(self):
        f=self.fields(1e-30,2.)
        for ratio in (-1,0,3.1,float('nan')):
            with self.assertRaises(ValueError):extrapolate_fields(f,f,ratio)
        f['packed_psi']=np.zeros(2)
        with self.assertRaises(ValueError):extrapolate_fields(f,f,1.)

if __name__=='__main__':unittest.main()
