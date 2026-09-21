"""Physical identities and supported-source constraints for the material audit."""
import math
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from prepare_templates_ldmos_fermi_materials import material_values

SOURCE='''
Bandgap {
 Eg0 = 1.12
 alpha = 0.0004
 beta = 600
 Tpar = 300
 Chi0 = 4.05
 dEg0(OldSlotboom) = 0
}
eDOSMass {
 Formula = 1
 a = 0.19
 ml = 0.91
 mm = 0
}
hDOSMass {
 * Ignore this equation containing braces: mh = {a/b}
 Formula = 1
 a = 0.5
 b = 0
 c = 0
 d = 0
 e = 0
 f = 0
 g = 0
 h = 0
 i = 0
 mm = 0
}
'''
CONSTANTS=dict(kb=1.380649e-23,q=1.602176634e-19,h=6.62607015e-34,m0=9.1093837015e-31)

class FermiMaterialsTest(unittest.TestCase):
    def test_density_identity_and_reference_temperature(self):
        r=material_values(SOURCE,CONSTANTS)
        self.assertEqual(r['bandgap_eV'],1.12)
        vt=CONSTANTS['kb']/CONSTANTS['q']*300
        log_ni=math.log(r['intrinsic_carrier_density_cm3'])
        expected=.5*(math.log(r['conduction_band_density_of_states_cm3'])+math.log(r['valence_band_density_of_states_cm3']))-1.12/(2*vt)
        self.assertAlmostEqual(log_ni,expected,places=13)
        pref=2*(2*math.pi*CONSTANTS['m0']*CONSTANTS['kb']*300/CONSTANTS['h']**2)**1.5/1e6
        self.assertAlmostEqual(r['valence_band_density_of_states_cm3']/pref,.5,places=14)

    def test_reject_unmapped_physics(self):
        for source in [SOURCE.replace('Formula = 1','Formula = 2'),SOURCE.replace('dEg0(OldSlotboom) = 0','dEg0(OldSlotboom) = -0.01595')]:
            with self.assertRaises(ValueError):material_values(source,CONSTANTS)

if __name__=='__main__':unittest.main()
