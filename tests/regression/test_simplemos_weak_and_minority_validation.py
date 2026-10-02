"""Numerical and contract properties for weak-column and minority-field validation."""
import csv
import json
import math
from pathlib import Path
import sys
import unittest
import numpy as np

REPO=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'scripts'))
from analyze_simplemos_weak_columns_20260906 import richardson
from analyze_simplemos_minority_invariance_20260906 import complete_qualified_group
import validate_simplemos_minority_invariance_20260906 as minority


class NumericalValidationTests(unittest.TestCase):
    def test_richardson_improves_smooth_nonlinear_reference(self):
        def central(h):return np.array([(math.exp(h)-math.exp(-h))/(2*h)])
        value,error=richardson(central(.04),central(.02),central(.01))
        self.assertLess(abs(value[0]-1),1e-8)
        self.assertLess(error,1e-7)

    def test_term_difference_recovers_charge_below_total_residual_precision(self):
        h=.001;p=1e-15;plus=p*math.exp(h);minus=p*math.exp(-h)
        whole=((1e5-plus)-(1e5-minus))/(2*h)
        component=-(plus-minus)/(2*h)
        self.assertEqual(whole,0.)
        self.assertAlmostEqual(component/(-p),1.,places=6)

    def test_symmetric_edge_difference_remains_conservative(self):
        derivative=np.longdouble('1e-18')
        assembled=np.array([derivative,-derivative],dtype=np.longdouble)
        self.assertEqual(sum(assembled),0.)

    def test_invariance_requires_all_four_qualified_initializations(self):
        group=[({'initialization':label,'strict_qualification_with_visibility':True},)
               for label in ('vela','native','hole_plus','hole_minus')]
        self.assertFalse(complete_qualified_group([]))
        self.assertFalse(complete_qualified_group(group[:2]))
        self.assertTrue(complete_qualified_group(group))
        group[-1][0]['strict_qualification_with_visibility']=False
        self.assertFalse(complete_qualified_group(group))

    def test_initialization_changes_preserve_contacts_and_density_relation(self):
        contract=minority.a.read(minority.OUT/'contract.json')
        self.assertEqual(len(contract['jobs']),32)
        for job in contract['jobs']:
            if job['initialization'] not in ('hole_plus','hole_minus'):continue
            path=Path(job['config']);cfg=minority.a.read(path)
            old=minority.a.read(Path(job['legacy_config']))
            self.assertEqual(cfg['contacts'],old['contacts'])
            geo=minority.d.matrix.spatial.m73.Geometry(job['device'])
            base=minority.d.ordered(path.parent.parent/'vela/initial.csv',geo.count)
            shifted=minority.d.ordered(path.parent/'initial.csv',geo.count)
            for i in geo.contact_nodes:self.assertEqual(base[i],shifted[i])
            offset=.05 if job['initialization']=='hole_plus' else -.05
            changed=[(x,y) for x,y in zip(base,shifted) if x!=y]
            self.assertTrue(changed)
            for x,y in changed:
                self.assertAlmostEqual(float(y['holes_m3'])/float(x['holes_m3']),math.exp(offset/minority.d.fixed.upstream.VT),places=12)
                self.assertEqual(x['electrons_m3'],y['electrons_m3'])
                self.assertEqual(x['psi'],y['psi'])

    def test_new_gate_is_explicit_and_does_not_relax_old_tolerance(self):
        for job in minority.a.read(minority.OUT/'contract.json')['jobs']:
            cfg=minority.a.read(Path(job['config']))['solver']['carrier_row_convergence']
            original=minority.a.read(Path(job['legacy_config']))['solver']['carrier_row_convergence']
            self.assertEqual(cfg['eps_row'],original['eps_row'])
            self.assertEqual(cfg['scale_floor'],0.)
            self.assertGreater(cfg['min_carrier_density_m3'],0.)
            self.assertEqual(original['min_carrier_density_m3'],0.)
            self.assertGreater(original['min_flux_scale'],cfg['min_flux_scale'])

    def test_each_minority_job_has_the_advertised_terminal_biases(self):
        for job in minority.a.read(minority.OUT/'contract.json')['jobs']:
            cfg=minority.a.read(Path(job['config']))
            biases={contact['name']:contact['bias'] for contact in cfg['contacts']}
            self.assertEqual(biases['gate'],job['vg'])
            self.assertEqual(biases['drain'],job['vd'])
            self.assertEqual(biases['source'],0.)
            self.assertEqual(biases['substrate'],0.)


if __name__=='__main__':unittest.main()
