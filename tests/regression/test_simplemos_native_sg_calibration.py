"""Physical identities independent of the native reference data."""
import math
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import audit_simplemos_masetti_box_mobility_20260907 as sg

class NativeSgCalibrationTests(unittest.TestCase):
 def state(self,psi,fn,fp):
  return 1e10*math.exp((psi-fn)/sg.VT),1e10*math.exp((fp-psi)/sg.VT)

 def flux(self,psi0,psi1,fn0,fn1,fp0,fp1,mode):
  n0,p0=self.state(psi0,fn0,fp0);n1,p1=self.state(psi1,fn1,fp1)
  return tuple(float(x) for x in sg.sg_flux(psi0,psi1,n0,n1,p0,p1,fn0,fn1,fp0,fp1,mode))

 def test_equilibrium_has_zero_flux_across_built_in_potential(self):
  self.assertEqual(self.flux(-.13,.21,.07,.07,.07,.07,'qf'),(0.,0.))

 def test_orientation_reversal_changes_both_flux_signs(self):
  f=self.flux(-.03,.11,.02,.09,.015,.025,'qf')
  r=self.flux(.11,-.03,.09,.02,.025,.015,'qf')
  for x,y in zip(f,r):self.assertLess(abs(x+y)/max(abs(x),abs(y)),2e-14)

 def test_density_and_electrochemical_forms_agree(self):
  f=self.flux(-.03,.11,.02,.09,.015,.025,'qf')
  r=self.flux(-.03,.11,.02,.09,.015,.025,'density')
  for x,y in zip(f,r):self.assertLess(abs(x/y-1),2e-14)

 def test_potential_gauge_shift_leaves_currents_unchanged(self):
  f=self.flux(-.03,.11,.02,.09,.015,.025,'qf')
  r=self.flux(.27,.41,.32,.39,.315,.325,'qf')
  for x,y in zip(f,r):self.assertLess(abs(x/y-1),2e-14)

 def test_conventional_electron_and_hole_signs_differ(self):
  e,h=self.flux(0.,0.,0.,.01,0.,.01,'qf')
  self.assertGreater(e,0.);self.assertLess(h,0.)
  self.assertLess(-sg.Q*e,0.);self.assertLess(sg.Q*h,0.)

if __name__=='__main__':unittest.main()
