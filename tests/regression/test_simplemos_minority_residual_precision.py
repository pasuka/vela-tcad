"""Physical invariants for the independent fixed-state precision reference."""
from decimal import Decimal as D,localcontext
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from analyze_simplemos_minority_residual_20260906 import em1,edge_flux,srh,exact


def edge():
    e={k:D(0) for k in ('node0','node1','epsi0','epsi1','ephi0','ephi1','hpsi0','hpsi1','hphi0','hphi1')}
    e.update(node1=D(1),ni0=D(1),ni1=D(1),ncoef=D(1),pcoef=D(1),scale=D(1),Vt=D(1));return e


class PrecisionTests(unittest.TestCase):
    def test_tiny_expm1_retains_signal(self):
        with localcontext() as c:
            c.prec=60;x=D('1e-80');self.assertEqual(x.exp()-1,0);self.assertEqual(em1(x),x)

    def test_equal_quasi_fermi_has_zero_flux_with_band_bending(self):
        with localcontext() as c:
            c.prec=100;e=edge();e.update(hpsi0=D(-20),hpsi1=D(10),epsi0=D(-20),epsi1=D(10))
            for carrier in ('electron','hole'):self.assertEqual(edge_flux(e,[],carrier,'kernel'),0)

    def test_hole_small_signal_has_expected_sign_and_magnitude(self):
        with localcontext() as c:
            c.prec=100;e=edge();e['hphi1']=D('1e-30');flux=edge_flux(e,[],'hole','kernel')
            self.assertLess(flux,0);self.assertLess(abs(flux/D('-1e-30')-1),D('1e-29'))

    def test_split_reference_recovers_drive_lost_by_reanchoring(self):
        with localcontext() as c:
            c.prec=100;e=edge();e['hphi0']=exact(.05);e['hphi1']=exact(.05)
            nodes=[dict(psi=D(0),href=D(0),hinc=exact(.05)),dict(psi=D(0),href=exact(.05),hinc=D('-1e-25'))]
            self.assertEqual(edge_flux(e,nodes,'hole','kernel'),0)
            self.assertGreater(edge_flux(e,nodes,'hole','partition'),0)

    def test_srh_equilibrium_and_linear_response(self):
        with localcontext() as c:
            c.prec=100;n=dict(ni=D(1),Vt=D(1),nsrh=D(1),p=D(1),dphi=D(0),taun=D(2),taup=D(3),volume=D(1),source_factor=D(1),scale=D(1))
            self.assertEqual(srh(n,'kernel'),0);n['dphi']=D('1e-30')
            self.assertLess(abs(srh(n,'kernel')/D('1e-31')-1),D('1e-29'))

    def test_60_and_100_digits_agree_on_nonlinear_flux(self):
        values=[]
        for precision in (60,100):
            with localcontext() as c:
                c.prec=precision;e=edge();e.update(hpsi0=D('.5'),hpsi1=D('.2'),hphi0=D('.1'),hphi1=D('.03'))
                values.append(edge_flux(e,[],'hole','kernel'))
        with localcontext() as c:
            c.prec=110;self.assertLess(abs(values[0]/values[1]-1),D('1e-55'))


if __name__=='__main__':unittest.main()
