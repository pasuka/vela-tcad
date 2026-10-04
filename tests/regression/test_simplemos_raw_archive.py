import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('raw', ROOT / 'scripts/simplemos_raw_archive.py')
raw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(raw)


class RawArchiveTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        # tempfile cleanup also needs the Windows extended-length root.
        self.tmp.name = str(raw.safe(self.root, 'child').parent)
        self.src = self.root / 'src'
        self.src.mkdir()
        self.content = b'\x89HDF\r\n\x1a\nbyte preservation\x00'
        (self.src / 'a.h5').write_bytes(self.content)
        s = (self.src / 'a.h5').stat()
        self.inventory = {'files': [{'path': 'a.h5', 'bytes': s.st_size, 'mtime_ns': s.st_mtime_ns}]}
        self.archive = self.root / 'raw.zip'

    def test_pack_verify_restore_and_no_overwrite(self):
        raw.pack(self.src, self.inventory, self.archive)
        r = raw.verify(self.archive)
        self.assertEqual(r['files_verified'], 1)
        dest = self.root / 'dest'
        self.assertEqual(raw.verify(self.archive, r['archive_sha256'], dest)['restored'], 1)
        self.assertEqual((dest / 'a.h5').read_bytes(), self.content)
        self.assertEqual(raw.verify(self.archive, r['archive_sha256'], dest)['restored'], 0)
        (dest / 'a.h5').write_bytes(b'preserve user file')
        with self.assertRaisesRegex(ValueError, 'Refusing'):
            raw.verify(self.archive, r['archive_sha256'], dest)
        self.assertEqual((dest / 'a.h5').read_bytes(), b'preserve user file')

    def test_changed_inventory_fails(self):
        (self.src / 'a.h5').write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'changed'):
            raw.pack(self.src, self.inventory, self.archive)

    def test_long_nested_path_roundtrip(self):
        name = '/'.join(['long_directory_segment'] * 15 + ['state.h5'])
        path = raw.safe(self.src, name)
        path.parent.mkdir(parents=True)
        path.write_bytes(self.content)
        s = path.stat()
        self.inventory['files'] = [{'path': name, 'bytes': s.st_size, 'mtime_ns': s.st_mtime_ns}]
        raw.pack(self.src, self.inventory, self.archive)
        result = raw.verify(self.archive, destination=self.root / 'dest')
        self.assertEqual(result['restored'], 1)
        self.assertEqual(raw.safe(self.root / 'dest', name).read_bytes(), self.content)

    def test_corruption_and_missing_member_fail_before_restore(self):
        data = raw.pack(self.src, self.inventory, self.archive)
        with self.assertRaisesRegex(ValueError, 'Archive SHA256'):
            raw.verify(self.archive, '0' * 64)
        with zipfile.ZipFile(self.archive, 'w') as z:
            z.writestr('manifest.json', json.dumps(data))
            z.writestr('files/a.h5', b'X' * len(self.content))
        with self.assertRaisesRegex(ValueError, 'Member SHA256'):
            raw.verify(self.archive, destination=self.root / 'dest')
        self.assertFalse((self.root / 'dest').exists())
        with zipfile.ZipFile(self.archive, 'w') as z:
            z.writestr('manifest.json', json.dumps(data))
        with self.assertRaisesRegex(ValueError, 'coverage'):
            raw.verify(self.archive)

    def test_unsafe_and_case_duplicate_paths_fail(self):
        for name in ('../escape', '/tmp/escape', 'a/../../x', '.git/config', 'x:stream', 'x\\y'):
            with self.assertRaises(ValueError):
                raw.safe(self.src, name)
        self.inventory['files'].append(dict(self.inventory['files'][0], path='A.H5'))
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            raw.pack(self.src, self.inventory, self.archive)


if __name__ == '__main__':
    unittest.main()
