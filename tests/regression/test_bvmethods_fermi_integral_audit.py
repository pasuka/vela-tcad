"""Run against an explicitly supplied private defining-integral executable."""
import os,subprocess,unittest
from decimal import Decimal as D,localcontext


class FermiIntegralAuditTests(unittest.TestCase):
    def evaluate(self, queries):
        p=subprocess.run([os.environ['VELA_FERMI_AUDIT_EXE']],input=queries,text=True,capture_output=True,timeout=60,check=True)
        return [[D(x) for x in line.split()] for line in p.stdout.splitlines()]

    def test_nondegenerate_series(self):
        with localcontext() as ctx:
            ctx.prec=70
            _,v,d,ev,ed=self.evaluate('f -40\n')[0];z=D(-40).exp()
            expected=z*(1-z/(D(2)*D(2).sqrt())+z*z/(D(3)*D(3).sqrt()))
            self.assertLess(abs(v/expected-1),D('1e-45'))
            self.assertLess(ev/abs(v),D('1e-32'))

    def test_density_and_compressibility_are_positive_monotone(self):
        rows=self.evaluate('f -2\nf 0\nf 8\n')
        self.assertTrue(all(row[1]>row[2]>0 for row in rows))
        self.assertTrue(all(a[1]<b[1] for a,b in zip(rows,rows[1:])))

    def test_analytic_integral_derivative_matches_direction(self):
        rows=self.evaluate('f 7.99999\nf 8\nf 8.00001\n')
        slope=(rows[2][1]-rows[0][1])/D('.00002')
        self.assertLess(abs(slope/rows[1][2]-1),D('1e-10'))

    def test_inverse_and_invalid_domain(self):
        with localcontext() as ctx:
            ctx.prec=70
            rows=self.evaluate('r 0.03\nr 18\n')
            self.assertLess(rows[0][0],rows[1][0])
            for row,target in zip(rows,[D('.03'),D(18)]):self.assertLess(abs(row[1]/target-1),D('1e-34'))
        p=subprocess.run([os.environ['VELA_FERMI_AUDIT_EXE']],input='r -1\n',text=True,capture_output=True)
        self.assertNotEqual(p.returncode,0)


if __name__=='__main__':unittest.main()
