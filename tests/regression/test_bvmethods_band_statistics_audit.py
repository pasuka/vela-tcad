from pathlib import Path
import sys,unittest
from decimal import Decimal as D,localcontext
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from bvmethods_band_statistics_audit import edge_driver,explicit_driver_flux


class BandStatisticsAuditTests(unittest.TestCase):
    def test_equilibrium_electrochemical_driver(self):
        for c,s in [('electron',-1),('hole',1)]:
            d=edge_driver(s*D('.025'),s*D('.05'),D(1),D(2),D('.025'),c)
            self.assertEqual(d['qf_drop_V'],0)
            self.assertEqual(explicit_driver_flux(D('1e18'),D('2e18'),d['qf_drop_V'],D('.025'),c),0)

    def test_band_reference_translation_leaves_edge_driver(self):
        a=edge_driver(D('.2'),D('.21'),D(2),D('2.1'),D('.025'),'electron')
        b=edge_driver(D('1000.2'),D('1000.21'),D(2),D('2.1'),D('.025'),'electron')
        self.assertEqual(a,b)

    def test_edge_reversal_and_tiny_nonzero_drive(self):
        with localcontext() as ctx:
            ctx.prec=80
            a=explicit_driver_flux(D('1e18'),D('2e18'),D('1e-20'),D('.025'),'electron')
            b=explicit_driver_flux(D('2e18'),D('1e18'),D('-1e-20'),D('.025'),'electron')
            self.assertNotEqual(a,0)
            self.assertLess(abs(a/-b-1),D('1e-55'))

    def test_invalid_physics_arguments(self):
        with self.assertRaises(ValueError):edge_driver(D(0),D(0),D(0),D(0),D(0),'electron')
        with self.assertRaises(ValueError):explicit_driver_flux(D(-1),D(1),D(0),D(1),'hole')
        with self.assertRaises(ValueError):edge_driver(D('NaN'),D(0),D(0),D(0),D(1),'electron')


if __name__=='__main__':unittest.main()
