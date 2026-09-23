import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from run_templates_ldmos_d5_targeted_controls import targets,secant_rows,set_target


class ControlsTest(unittest.TestCase):
    def test_exact_endpoint_and_positive_increments(self):
        t=targets(20,20+4/3,.2)
        self.assertEqual(len(t),7)
        self.assertEqual(t[-1],20+4/3)
        self.assertTrue(all(0<b-a<=.2+1e-12 for a,b in zip([20]+t,t)))

    def test_qf_cap_is_independent_and_reference_is_preserved(self):
        cfg={'contacts':[{'name':'drain','bias':0},{'name':'source','bias':-28}], 'sweep':{},'solver':{}}
        set_target(cfg,33,.8,.2,28,'seed.csv')
        self.assertEqual(cfg['sweep']['bias_points'],[5])
        self.assertEqual(cfg['sweep']['max_step'],.8)
        self.assertEqual(cfg['solver']['quasi_fermi_update_limit_V'],.2)
        self.assertEqual(cfg['contacts'][1]['bias'],-28)

    def test_secant_state_identity_and_guard(self):
        old=[dict(node_id='1',psi='1',phin='2',phip='3')]
        cur=[dict(node_id='1',psi='2',phin='3',phip='4')]
        out=secant_rows(cur,old,.5)[0]
        self.assertEqual(float(out['psi']),2.5)
        self.assertEqual(out['electron_qf_reference_V'],out['phin'])
        self.assertEqual(out['electron_qf_increment_V'],'0')
        self.assertEqual(cur[0]['psi'],'2')
        with self.assertRaises(ValueError):secant_rows(cur,old,3)
        with self.assertRaises(ValueError):secant_rows(cur,[dict(old[0],node_id='2')],1)


if __name__=='__main__':unittest.main()
