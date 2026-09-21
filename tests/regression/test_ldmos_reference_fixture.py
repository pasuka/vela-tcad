"""Portable neutral inputs must not silently depend on ignored staging files."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from run_templates_ldmos_reference import verify_fixture, FIXTURE
from run_templates_ldmos_stage4_d5 import read_points, assert_d5_contract
from run_templates_ldmos_linked_d5 import read, digest


class ReferenceFixtureTests(unittest.TestCase):
    def test_all_runtime_inputs_are_tracked_portable_and_hash_qualified(self):
        bundle=verify_fixture()
        self.assertEqual(len(bundle['files']),9)
        template=read(ROOT/bundle['template'])
        assert_d5_contract(template)
        text=json.dumps(template)
        self.assertNotIn('reference_staging',text)
        for gate in ('4','8'):
            points=read_points(ROOT/bundle['references'][gate])
            self.assertEqual(len(points),31)
            self.assertEqual((points[0],points[-1]),(0.,40.))
            self.assertIn(bundle['seeds'][gate],bundle['files'])
        provenance=read(FIXTURE.parent/'provenance.json')
        for item in provenance['files']:
            self.assertEqual(digest(FIXTURE.parent/item['path']),item['sha256'])

    def test_external_or_modified_dependency_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fixture=root/'case';fixture.mkdir()
            p=root/'external.csv';p.write_text('outside')
            b=fixture/'inputs.json';b.write_text(json.dumps(dict(files={'external.csv':digest(p)})))
            with self.assertRaises(ValueError):verify_fixture(root,b)
            p=fixture/'state.csv';p.write_text('new')
            b.write_text(json.dumps(dict(files={'case/state.csv':'not-the-hash'})))
            with self.assertRaises(ValueError):verify_fixture(root,b)


if __name__=='__main__':unittest.main()
