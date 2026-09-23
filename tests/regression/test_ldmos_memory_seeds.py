import hashlib
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from ldmos_memory_seeds import MemorySeeds,replay_prediction
import dd_state_binary as binary


class MemorySeedTest(unittest.TestCase):
    def test_bounded_storage_spills_and_failure_flush_preserves_exact_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            a,b=Path(tmp)/'a.vds',Path(tmp)/'b.vds';store=MemorySeeds(4)
            store.add(a,b'123');store.add(b,b'456')
            self.assertEqual(a.read_bytes(),b'123');self.assertFalse(b.exists())
            self.assertEqual(store.audit()['memory_seeds'],1)
            store.flush();self.assertEqual(b.read_bytes(),b'456')
            self.assertEqual(store.audit()['memory_seeds'],0)

    def test_immutable_seed_and_cold_spill_audit(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'seed.vds';store=MemorySeeds(2)
            store.add(path,b'abc')
            with self.assertRaises(ValueError):store.add(path,b'def')
            path.write_bytes(b'bad')
            with self.assertRaises(ValueError):store.audit()

    def test_replay_uses_persisted_parents_and_checks_hash(self):
        row=dict(node_id='0',psi=.2,phin=.1,phip=-.1,electrons_m3=1e12,holes_m3=1e10,
                 electron_qf_increment_V=0.,hole_qf_increment_V=0.,electron_qf_reference_V=.1,hole_qf_reference_V=-.1)
        old=dict(row,psi=.1,phin=.05,phip=-.05,electron_qf_reference_V=.05,hole_qf_reference_V=-.05)
        # Expected physical secant doubles the previous change at ratio 1.
        expected=dict(row,psi=float(.2+(.2-.1)),phin=float(.1+(.1-.05)),phip=float(-.1+(-.1+.05)))
        expected.update(electron_qf_reference_V=expected['phin'],hole_qf_reference_V=expected['phip'])
        metadata=dict(mode='outer_secant_guarded_v1',ratio=1.,current_state='new',previous_state='old',
                      current_sha256='new',previous_sha256='old',predicted_sha256=hashlib.sha256(binary.encode([expected])).hexdigest())
        out=replay_prediction(metadata,lambda p:[dict(row if str(p)=='new' else old)],lambda p:str(p))
        self.assertEqual(binary.decode(out)[0]['psi'],expected['psi'])
        with self.assertRaises(ValueError):replay_prediction(metadata,lambda p:[row],lambda p:'changed')


if __name__=='__main__':unittest.main()
