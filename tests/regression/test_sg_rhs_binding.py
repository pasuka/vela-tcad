import importlib.util
from pathlib import Path
import unittest
spec=importlib.util.spec_from_file_location('sg_binding',Path(__file__).resolve().parents[2]/'scripts/audit_sg_rhs_binding.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class Binding(unittest.TestCase):
    def edges(self):return [{'cell':0,'nodes':[0,1],'product':3,'factor':1,'snapshot':'rhs1'},
                           {'cell':0,'nodes':[2,0],'product':1,'factor':1,'snapshot':'rhs1'}]
    def test_signed_conservation_and_sources(self):
        r=m.bind_row(0,{0:[0,1,2]},self.edges(),[{'value':-1,'snapshot':'rhs1'}],-1,'rhs1')
        self.assertEqual(r['transport'],2);self.assertEqual(r['error'],0)
    def test_reverse_edge_preserves_row(self):
        edges=self.edges();edges[0].update(nodes=[1,0],product=-3)
        self.assertEqual(m.bind_row(0,{0:[0,1,2]},edges,[],-2,'rhs1')['error'],0)
    def test_zero_edge_is_required(self):
        edges=self.edges();edges[0]['product']=0
        self.assertEqual(m.bind_row(0,{0:[0,1,2]},edges,[],1,'rhs1')['error'],0)
        with self.assertRaises(ValueError):m.bind_row(0,{0:[0,1,2]},edges[1:],[],1,'rhs1')
    def test_duplicate_and_foreign_support_rejected(self):
        for edges in [self.edges()+[self.edges()[0]],[dict(e,cell=1)for e in self.edges()]]:
            with self.assertRaises(ValueError):m.bind_row(0,{0:[0,1,2]},edges,[],0,'rhs1')
    def test_mixed_snapshot_rejected(self):
        for edges,sources in [([dict(e,snapshot='rhs2')for e in self.edges()],[]),(self.edges(),[{'value':1,'snapshot':'rhs2'}])]:
            with self.assertRaises(ValueError):m.bind_row(0,{0:[0,1,2]},edges,sources,0,'rhs1')
    def test_zero_budget_does_not_hide_residual(self):
        edges=[dict(e,product=0)for e in self.edges()]
        self.assertIsNone(m.bind_row(0,{0:[0,1,2]},edges,[],1,'rhs1')['error'])
    def test_nonfinite_and_wrong_mesh_rejected(self):
        for value in ['NaN','Infinity',True]:
            with self.assertRaises(ValueError):m.finite(value)
        with self.assertRaises(ValueError):m.bind_row(0,{0:[0,0,1]},self.edges(),[],0,'rhs1')

if __name__=='__main__':unittest.main()
