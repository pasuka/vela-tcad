import json
import re
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from prepare_wp14_remaining_references import materialize_bvmethods, materialize_gate, SOURCE, VARIANTS


class RemainingSourcesTest(unittest.TestCase):
    def test_native_precision_repair_preserves_physics_contacts_goals_and_stop(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            original = output / 'original'
            repaired = output / 'repaired'
            materialize_bvmethods(original)
            materialize_bvmethods(repaired, numerical_precision=True)
            for name in ('d_doping_enormal', 'external_resistor'):
                old = (original / (name + '_sdevice.cmd')).read_text()
                new = (repaired / (name + '_sdevice.cmd')).read_text()
                for section in ('Physics', 'Electrode', 'Goal', 'BreakCriteria'):
                    pattern = rf'\b{section}\s*\{{(.*?)\n?\}}'
                    self.assertEqual(re.findall(pattern, old, re.S), re.findall(pattern, new, re.S))
                if name == 'd_doping_enormal':
                    self.assertIn('Digits=8', new)
                else:
                    self.assertEqual(old, new)
            external = (repaired / 'external_resistor_sdevice.cmd').read_text()
            self.assertNotIn('ExtendedPrecision', external)
            self.assertIn('Resistor=1e7', external)
            self.assertIn('Voltage=100000.0', external)
            self.assertIn('AbsVal=1.443e-3', external)

    def test_bvmethods_preserves_exact_sources_and_declared_ablation_branches(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            manifest = materialize_bvmethods(output)
            for name in ('models.par', 'full_physics_constant_srh.par',
                         'aba_poisson_sdevice.cmd', 'external_resistor_sdevice.cmd'):
                self.assertEqual((SOURCE / name).read_bytes(), (output / name).read_bytes())
            for name, (srh, normal) in VARIANTS.items():
                text = (output / (name + '_sdevice.cmd')).read_text()
                self.assertNotIn('SRH_PARAMETER_FILE', text)
                self.assertNotIn('SRH_MODEL', text)
                self.assertNotIn('ENORMAL_MODEL', text)
                self.assertEqual(srh, 'SRH(DopingDep)' in text)
                self.assertEqual(normal, 'HighFieldsaturation(GradQuasiFermi) Enormal' in text)
                self.assertIn('Voltage=6.0', text)
                self.assertIn('Current=1e-4', text)
                self.assertIn('Fermi', text)
                self.assertIn('Band2Band(E2)', text)
                self.assertIn('Avalanche(Eparallel)', text)
            self.assertEqual('bvmethods_sentaurus2018/v1/mesh/n1_msh.tdr', manifest['mesh_origin'])

    def test_gate_uses_one_geometry_and_current_library_for_native_reference_controls(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            manifest = materialize_gate(output)
            self.assertEqual(1e17, manifest['donors_cm3'])
            self.assertEqual(0.01, manifest['geometry_um']['oxide_thickness'])
            self.assertEqual((ROOT / 'share/vela/materials/models.par').read_bytes(),
                             (output / 'models.par').read_bytes())
            self.assertEqual(4.611003610676079, manifest['reference_intrinsic_work_function_eV'])
            expected = manifest['variants']['barrier_zero']['expected_flatband_V']
            self.assertAlmostEqual(0.55, manifest['variants']['barrier_positive']['expected_flatband_V'] - expected)
            self.assertAlmostEqual(-0.3, manifest['variants']['barrier_negative']['expected_flatband_V'] - expected)
            self.assertAlmostEqual(4.633743995291464, manifest['sentaurus_gate_reference_eV'])
            self.assertAlmostEqual(-0.02274038461538463,
                manifest['variants']['workfunction']['expected_flatband_V'] - expected)
            self.assertEqual(expected, manifest['variants']['workfunction_zero']['expected_flatband_V'])
            for name, entry in manifest['variants'].items():
                text = (output / (name + '_sdevice.cmd')).read_text()
                self.assertIn('Grid="gate_msh.tdr"', text)
                self.assertIn(entry['attribute'], text)
                self.assertNotIn('Fermi', text.split('Plot')[0])
                self.assertEqual(41, entry['samples'])
            fixture = json.loads((output / 'vela_fixture.json').read_text())['mesh']
            self.assertEqual({'Si','SiO2'}, {r['material'] for r in fixture['regions']})


if __name__ == '__main__':
    unittest.main()
