"""Prevent signed Poisson override from silently bypassing SLOT AverageBox."""
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]


class AverageBoxColdConfig(unittest.TestCase):
    def setUp(self):
        self.config = json.loads((ROOT / 'configs/regression/slot_ldmos_idvd_average_box_aligned_control.json').read_text())

    def test_dielectric_and_charge_use_box_geometry(self):
        c = self.config
        self.assertEqual(c['mesh_geometry']['cell_box_policy'], 'average_intersection')
        self.assertEqual(c['mesh_geometry']['poisson_permittivity_policy'], 'cell_material')
        self.assertEqual(c['discretization']['poisson_permittivity_coupling'], 'edge_mean')
        self.assertEqual(c['discretization']['poisson_charge_volume_policy'], 'material_local_box')

    def test_new_cold_lineage(self):
        c = self.config
        self.assertTrue(all(x['bias'] == 0. for x in c['contacts']))
        self.assertNotIn('initial_state_file', c['sweep'])
        self.assertNotIn('state_file', c)
        self.assertEqual(c['sweep']['bias_points'], [0.])
        self.assertEqual(c['sweep']['initialization']['mode'], 'poisson_block')

    def test_full_qualification_preserved(self):
        s = self.config['solver']
        self.assertEqual(s['abstol'], 1e-9)
        self.assertEqual(s['reltol'], 1e-8)
        self.assertEqual(s['carrier_statistics'], 'fermi_dirac')
        self.assertEqual(s['carrier_row_convergence']['mode'], 'enforce')
        self.assertEqual(s['carrier_row_convergence']['eps_row'], .001)
        self.assertEqual(s['global_continuity_closure']['mode'], 'enforce')
        self.assertEqual(s['global_continuity_closure']['tolerance'], .01)


if __name__ == '__main__':
    unittest.main()
