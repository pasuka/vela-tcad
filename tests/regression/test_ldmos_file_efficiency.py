import csv
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest
from types import SimpleNamespace

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from ldmos_file_efficiency import CompletedFileCache,EfficiencyPolicy


class FileReuseTest(unittest.TestCase):
    def test_predictor_mutations_cannot_contaminate_cached_state(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'state.csv';p.write_text('node_id,psi\n1,2\n')
            def rows(p):
                with p.open() as f:return list(csv.DictReader(f))
            cache=CompletedFileCache(rows,lambda p:hashlib.sha256(p.read_bytes()).hexdigest())
            cache.rows(p);self.assertEqual(cache.stats['csv_bypass'],1)
            cache.seal(p.parent)
            first=cache.rows(p);first[0]['psi']='123';first.append({'psi':'999'})
            self.assertEqual(cache.rows(p),[{'node_id':'1','psi':'2'}])
            p.write_text('node_id,psi\n1,345\n')
            self.assertEqual(cache.rows(p)[0]['psi'],'345')
            self.assertEqual(cache.stats['csv_invalidated'],1)
            cache.disable();p.write_text('node_id,psi\n1,7\n')
            self.assertEqual(cache.rows(p)[0]['psi'],'7')
            self.assertEqual(len(cache.csv),0)

    def test_digest_replacement_invalidation_and_cold_verification(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a';p.write_bytes(b'a')
            c=CompletedFileCache(None,lambda p:hashlib.sha256(p.read_bytes()).hexdigest())
            h=c.digest(p);self.assertEqual(c.digest(p),h)
            other=p.with_suffix('.tmp');other.write_bytes(b'b');other.replace(p)
            self.assertNotEqual(c.digest(p),h)
            c.disable();self.assertEqual(c.digest(p),hashlib.sha256(b'b').hexdigest())
            self.assertEqual(c.stats['digest_hit'],1)

    def test_csv_capacity_is_bounded(self):
        with tempfile.TemporaryDirectory() as d:
            c=CompletedFileCache(lambda p:[{'value':p.read_text()}],None,entries=2);c.seal(d)
            for i in range(5):
                p=Path(d)/str(i);p.write_text(str(i));c.rows(p)
            self.assertEqual(len(c.csv),2)

    def test_write_during_read_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'state.csv';p.write_text('old')
            def rows(p):
                p.write_text('new state written during read');return [{'value':'old'}]
            c=CompletedFileCache(rows,None);c.seal(d)
            with self.assertRaises(RuntimeError):c.rows(p)


class EfficiencyTest(unittest.TestCase):
    def test_cheaper_per_volt_can_grow_despite_more_single_point_work(self):
        p=EfficiencyPolicy();p.next(.1,.1,2.,1.,False,0)
        step,detail=p.next(.2,.2,3.,1.,False,0)
        self.assertGreater(step,.2);self.assertEqual(detail['cost_per_V'],15.)

    def test_cost_per_volt_regression_reverses_size_change(self):
        for h,c in ((.2,6.),(.05,2.)):
            p=EfficiencyPolicy();p.next(.1,.1,2.,1.,False,0)
            step,detail=p.next(h,h,c,1.,False,0)
            self.assertEqual(step,.1);self.assertIn('cost_regression',detail['reason'])

    def test_fragments_frames_and_recovery_do_not_train_cost_trend(self):
        p=EfficiencyPolicy();p.next(.1,.1,2.,1.,False,0)
        self.assertEqual(p.next(.2,.02,1.,1.,False,0)[0],.2)
        self.assertIsNone(p.previous)
        p.next(.1,.1,2.,1.,False,0)
        self.assertGreater(p.next(.2,.2,6.,1.,False,28)[0],.2)
        self.assertEqual(p.next(.2,.2,6.,1.,True,28)[0],.1)
        self.assertIsNone(p.previous)

    def test_uniform_cost_unit_change_preserves_decision(self):
        a=EfficiencyPolicy();b=EfficiencyPolicy()
        for h,c in ((.1,2.),(.2,3.),(.27,6.)):
            self.assertEqual(a.next(h,h,c,1,False,0)[0],b.next(h,h,10*c,1,False,0)[0])

    def test_failed_outer_attempt_is_charged_to_accepted_progress(self):
        from run_templates_ldmos_d5_cost_controls import install_candidate
        class Sweep:
            def __init__(self):self.ledger={'runs':[],'frame_V':0}
            def save(self):pass
            def child(self,*args):
                self.ledger['runs'].append({})
                return dict(cpu_seconds={'total':2.},curve=[]),Path('unused')
        b=SimpleNamespace(Sweep=Sweep,MAXIMUM=.8,read=lambda p:{'counters':{}},rows=lambda p:[])
        s=install_candidate(b,efficiency=True)()
        # First target attempt is not accepted; feedback only occurs after retry.
        s.child(None,0,.2,.2,'direct');s.child(None,0,.1,.1,'direct')
        s.feedback('fixed',.1,.1,0)
        row=s.ledger['cost_feedback'][-1]
        self.assertAlmostEqual(row['efficiency']['estimated_total_cost_seconds'],4.76)
        self.assertEqual(len(row['attempts']),2)
        self.assertEqual(s.pending,[])

if __name__=='__main__':unittest.main()
