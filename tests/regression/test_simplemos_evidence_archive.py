"""Integrity and safe restoration of external historical evidence."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import warnings
import zipfile

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('archive', ROOT / 'scripts/simplemos_evidence_archive.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


class EvidenceArchiveTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.names = [a.PREFIX + 'sample/summary.json', a.PREFIX + 'sample/nodes.csv']
        self.contents = [b'{"passed":false}\n', b'node,value\n0,1\n']
        self.archive = self.root / 'archive.zip'
        with zipfile.ZipFile(self.archive, 'w') as z:
            for name, data in zip(self.names, self.contents):
                z.writestr(name, data)
        self.data = {'schema': 'vela.simplemos.local-evidence.v1', 'file_count': 2,
                     'archive': {}, 'files': [
                         {'path': n, 'bytes': len(d), 'sha256': hashlib.sha256(d).hexdigest(),
                          'storage': 'git' if i == 0 else 'archive'}
                         for i, (n, d) in enumerate(zip(self.names, self.contents))]}
        self.reseal_archive()

    def reseal_archive(self):
        self.data['archive'] = {'bytes': self.archive.stat().st_size, 'sha256': a.sha(self.archive)}

    def test_restore_exact_bytes_and_identical_noop(self):
        dst = self.root / 'restore'
        result = a.restore(self.data, self.archive, dst)
        self.assertEqual(result['created'], 1)
        self.assertFalse((dst / self.names[0]).exists())
        self.assertEqual((dst / self.names[1]).read_bytes(), self.contents[1])
        stamp = (dst / self.names[1]).stat().st_mtime_ns
        self.assertEqual(a.restore(self.data, self.archive, dst)['existing_identical'], 1)
        self.assertEqual((dst / self.names[1]).stat().st_mtime_ns, stamp)

    def test_conflict_is_rejected_before_any_write(self):
        dst = self.root / 'restore'
        path = dst / self.names[1]
        path.parent.mkdir(parents=True)
        path.write_bytes(b'user work')
        with self.assertRaisesRegex(ValueError, 'overwrite'):
            a.restore(self.data, self.archive, dst, self.names)
        self.assertFalse((dst / self.names[0]).exists())
        self.assertEqual(path.read_bytes(), b'user work')

    def test_bad_outer_hash_rejected(self):
        with self.archive.open('ab') as f:
            f.write(b'changed')
        with self.assertRaisesRegex(ValueError, 'byte identity'):
            a.verify_archive(self.data, self.archive)

    def test_bad_member_hash_even_after_outer_reseal(self):
        with zipfile.ZipFile(self.archive, 'w') as z:
            z.writestr(self.names[0], self.contents[0])
            z.writestr(self.names[1], b'node,value\n0,2\n')
        self.reseal_archive()
        with self.assertRaisesRegex(ValueError, 'member hash'):
            a.restore(self.data, self.archive, self.root / 'restore')
        self.assertFalse((self.root / 'restore').exists())

    def test_duplicate_archive_member_rejected(self):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            with zipfile.ZipFile(self.archive, 'a') as z:
                z.writestr(self.names[1], self.contents[1])
        self.reseal_archive()
        with self.assertRaisesRegex(ValueError, 'coverage'):
            a.verify_archive(self.data, self.archive)

    def test_traversal_and_unknown_selection_rejected(self):
        for name in (a.PREFIX + '../escape', a.PREFIX + 'x/../../escape',
                     'D:/outside', a.PREFIX + 'x\\escape', a.PREFIX + 'x:stream'):
            with self.assertRaises(ValueError):
                a.safe_path(self.root, name)
        with self.assertRaisesRegex(ValueError, 'selection'):
            a.restore(self.data, self.archive, self.root, [a.PREFIX + 'unknown.csv'])

    def test_duplicate_manifest_casefold_rejected(self):
        self.data['files'][1]['path'] = self.names[0].upper().replace(a.PREFIX.upper(), a.PREFIX)
        path = self.root / 'manifest.json'
        path.write_text(json.dumps(self.data))
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            a.load_manifest(path)

    def test_missing_archived_copy_is_optional_but_retained_is_required(self):
        dst = self.root / 'checkout'
        with self.assertRaisesRegex(ValueError, 'Missing'):
            a.verify_files(self.data, dst)
        a.restore(self.data, self.archive, dst, [self.names[0]])
        result = a.verify_files(self.data, dst)
        self.assertEqual(result, {'retained_verified': 1, 'local_archive_copies_verified': 0})

    def test_shipped_manifest_and_precise_ignore_scope(self):
        data = a.load_manifest(a.MANIFEST)
        result = a.verify_files(data, ROOT)
        self.assertEqual(data['file_count'], 1946)
        self.assertEqual(result['retained_verified'], 1424)
        expected = {'/' + x['path'][len(a.PREFIX):]
                    for x in data['files'] if x['storage'] == 'archive'}
        lines = (ROOT / a.PREFIX / '.gitignore').read_text().splitlines()
        actual = {line for line in lines if line and not line.startswith('#')}
        self.assertEqual(actual, expected)
        self.assertEqual(len(actual), 522)

    def test_remote_receipt_matches_frozen_manifest(self):
        data = a.load_manifest(a.MANIFEST)
        receipt = json.loads((a.MANIFEST.parent / 'remote_verification.json').read_text())
        self.assertTrue(receipt['passed'])
        self.assertEqual(receipt['manifest_sha256'], a.sha(a.MANIFEST))
        self.assertEqual(receipt['archive_sha256'], data['archive']['sha256'])
        self.assertEqual(receipt['archive_bytes'], data['archive']['bytes'])
        self.assertEqual(receipt['archive_path'], data['archive']['remote_path'])
        self.assertEqual(receipt['members_verified'], data['file_count'])
        self.assertEqual(receipt['uncompressed_bytes'], sum(x['bytes'] for x in data['files']))


if __name__ == '__main__':
    unittest.main()
