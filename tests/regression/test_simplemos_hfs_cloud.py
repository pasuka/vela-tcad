"""Acceptance failure-path checks for the portable HFS experiment harness."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import simplemos_hfs_cloud_20260926 as h


class Acceptance(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.geo=dict(count=2,free_si=[1],contacts=[0])
        self.terms=[dict(node_id=i,**{c+'_'+k:v for c in ('electron','hole')
                    for k,v in dict(flux_abs_sum=1.,recombination=0.,impact=0.,residual=0.).items()}) for i in (0,1)]
        h.csvout(self.root/'all_row.csv',self.terms)
        h.csvout(self.root/'acceptance_edges.csv',[dict(node0=0,node1=1,electron_flux=0.,hole_flux=0.)])
        for name in ('config','all_row','acceptance_edges','functional'):
            h.write(self.root/(name+'.status.json'),dict(exit_code=0,converged=True,
                    contact_currents_A_per_um=dict(drain=1.,source=-1.),
                    current_A_per_um=1.,contact_current_extractor_A_per_um=1.))

    def check(self):return h.qualify(self.root,self.geo)

    def test_closed_state_passes_all_rows(self):
        r=self.check();self.assertTrue(r['qualified']);self.assertEqual(r['carrier_rows'],2)

    def test_small_terminal_error_cannot_hide_failed_row(self):
        self.terms[1]['hole_residual']=2e-6
        h.csvout(self.root/'all_row.csv',self.terms)
        r=self.check();self.assertFalse(r['qualified']);self.assertEqual(r['row_violations'],1)

    def test_nonfinite_row_rejected(self):
        self.terms[1]['hole_residual']='nan'
        h.csvout(self.root/'all_row.csv',self.terms)
        self.assertFalse(self.check()['qualified'])

    def test_active_source_requires_global_closure(self):
        self.terms[1]['electron_recombination']=1e-9
        h.csvout(self.root/'all_row.csv',self.terms)
        r=self.check();self.assertFalse(r['qualified']);self.assertTrue(r['closure']['electron']['active'])

    def test_inactive_source_keeps_frozen_floor(self):
        self.terms[1]['electron_recombination']=1e-12
        h.csvout(self.root/'all_row.csv',self.terms)
        self.assertTrue(self.check()['qualified'])

    def test_kcl_and_port_are_independent_gates(self):
        p=self.root/'config.status.json';s=h.read(p);s['contact_currents_A_per_um']['source']=-1+1e-7;h.write(p,s)
        self.assertFalse(self.check()['qualified'])
        s['contact_currents_A_per_um']['source']=-1.;h.write(p,s)
        p=self.root/'functional.status.json';s=h.read(p);s['current_A_per_um']=1+2e-8;h.write(p,s)
        self.assertFalse(self.check()['qualified'])

    def test_missing_or_duplicate_node_is_not_accepted(self):
        h.csvout(self.root/'all_row.csv',[self.terms[0],self.terms[0]])
        with self.assertRaises(AssertionError):self.check()

    def test_failed_probe_is_not_accepted(self):
        p=self.root/'acceptance_edges.status.json';s=h.read(p);s['exit_code']=1;h.write(p,s)
        self.assertFalse(self.check()['qualified'])

    def test_dual_initialization_rejects_field_difference_despite_equal_current(self):
        state=[dict(node_id=i,psi=0.,phin=0.,phip=0.,electrons_m3=1e20,holes_m3=1e5) for i in (0,1)]
        a=self.root/'a';b=self.root/'b';h.csvout(a/'state.csv',state)
        changed=copy.deepcopy(state);changed[1]['phip']=2e-6;h.csvout(b/'state.csv',changed)
        x=dict(dest=str(a),qualified=True,current_A_per_um=1.)
        y=dict(dest=str(b),qualified=True,current_A_per_um=1.)
        self.assertFalse(h.dual(x,y,self.geo)['qualified'])


if __name__=='__main__':unittest.main()
