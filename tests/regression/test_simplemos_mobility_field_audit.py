import importlib.util
import math
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
def module(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/(name+'.py'))
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result

a=module('audit_simplemos_mobility_fields_20260930')
p=module('prepare_simplemos_mobility_native_20260930')
s=module('summarize_simplemos_mobility_fields_20260930')


class MobilityFieldAuditTests(unittest.TestCase):
    def test_summary_rejects_incomplete_and_duplicate_coverage(self):
        for records in ([],[{'case':'n17_vd_1','index':0}]*816):
            with self.assertRaises(ValueError):s.summarize(records)

    def test_hole_failure_remains_separate_from_electron_pass(self):
        records=[]
        for n in range(17,25):
            for vd in ('0p05','1'):
                for i in range(51):
                    carriers={}
                    for car in ('electron','hole'):
                        carriers[car]={metric:{'max':1e-8} for metric in s.METRICS}
                        carriers[car].update(same_export_cell_gate_1e_7_passed=car=='electron',
                            same_export_worst_cell={},native_node_display_hypotheses={m:{'max':.1} for m in ('equal','cell_area','box_vertex')})
                    records.append(dict(case=f'n{n}_vd_{vd}',index=i,vg=i*.05,
                                        native_load_identity={'potential':0.},carriers=carriers))
        result=s.summarize(records)
        self.assertEqual(result['carriers']['electron']['same_export_cell_gate_1e_7']['passed_points'],816)
        self.assertEqual(result['carriers']['hole']['same_export_cell_gate_1e_7']['failed_points'],816)

    def test_uniform_cell_mobility_preserved(self):
        parts={(0,1):[(0,.3),(1,.7)],(1,2):[(0,2.)]}
        values=a.edge_reconstruct(parts,{0:3.,1:3.})
        self.assertEqual(set(values),{(0,1),(1,2)})
        for value in values.values():self.assertLessEqual(abs(value-3.),math.ulp(3.))

    def test_geometric_average_not_endpoint_average(self):
        parts={(0,1):[(0,1.),(1,3.)],(1,2):[(1,0.)]}
        self.assertEqual(a.edge_reconstruct(parts,{0:2.,1:10.}),{(0,1):8.})

    def test_nonfinite_cannot_be_reported_as_small(self):
        for data in ([],[math.nan],[math.inf]):
            with self.assertRaises(ValueError):a.stats(data)

    def test_duplicate_scalar_indices_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            path=Path(name)/'scalar.csv';path.write_text('node_id,component0\n0,1\n0,2\n')
            with self.assertRaises(ValueError):a.scalar(path)

    def test_native_export_covers_exactly_816_without_solve(self):
        with tempfile.TemporaryDirectory() as name:
            out=Path(name)/'native';p.prepare(out)
            decks=list(out.glob('n*/native.cmd'));self.assertEqual(len(decks),8)
            for deck in decks:
                text=deck.read_text()
                self.assertEqual(text.count('Load(FilePrefix='),102)
                self.assertEqual(text.count('Plot(FilePrefix='),102)
                self.assertNotIn('Coupled',text);self.assertNotIn('Quasistationary',text)
                for vd in ('0p05','1'):
                    for i in range(51):
                        self.assertIn(f'{deck.parent.name}_vd_{vd}_state_{i:04d}',text)
            with self.assertRaises(ValueError):p.prepare(out)


if __name__=='__main__':unittest.main()
