"""IDL material bundle, cold lineage and acceptance contract (no native data required)."""
import copy
import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('idl', ROOT/'scripts/run_slot_ldmos_idl.py')
idl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(idl)


class IDLContract(unittest.TestCase):
    def test_material_reference_and_intrinsic_density_are_coherent(self):
        import math
        material = json.loads((ROOT/'configs/regression/slot_ldmos_idl_bennett_materials.json').read_text())
        p = next(x for x in material['materials'] if x['name'] == 'PolySilicon')
        vt = .025851995266484913
        ni = math.sqrt(p['Nc_m3']*p['Nv_m3'])*math.exp(-p['bandgap_eV']/(2*vt))
        self.assertLess(abs(ni/p['ni']-1), 3e-5)
        vacuum = p['electron_affinity_eV'] + p['bandgap_eV']/2 + vt/2*math.log(p['Nc_m3']/p['Nv_m3'])
        self.assertLess(abs((4.633703995135664-vacuum)-p['intrinsic_reference_potential_V']), 1e-10)
        self.assertNotIn('ni_m3', p)

    def test_cold_configuration_and_original_gates(self):
        c = json.loads((ROOT/'configs/regression/slot_ldmos_idl_bennett_control.json').read_text())
        self.assertNotIn('initial_state_file', c['sweep'])
        self.assertEqual(c['sweep']['bias_points'], [0.])
        self.assertTrue(all(x['bias']==0 for x in c['contacts']))
        self.assertEqual({x['name'] for x in c['contacts']}, {'gate','drain','source','substrate'})
        self.assertEqual(c['solver']['bandgap_narrowing']['model'], 'bennett_wilson')
        self.assertEqual(c['solver']['impact_ionization']['model'], 'none')
        self.assertEqual(c['solver']['recombination_materials'], ['Si'])
        self.assertEqual(c['solver']['carrier_row_convergence']['eps_row'], .001)
        self.assertEqual(c['solver']['global_continuity_closure']['tolerance'], .01)
        self.assertEqual(c['solver']['abstol'], 1e-9)

    def test_successful_process_is_not_sufficient_for_acceptance(self):
        row = dict(converged='1',solver_method='newton',carrier_row_violations='0',
            global_continuity_closure_satisfied='1',qf_bounds_violations='0',
            carrier_row_max_ratio='.0009',current_total_A_per_um='1e-8',bias_V='2.6')
        self.assertTrue(idl.qualify(row,2.6))
        for key,value in [('converged','0'), ('carrier_row_violations','1'),
                          ('global_continuity_closure_satisfied','0'),
                          ('qf_bounds_violations','1'), ('carrier_row_max_ratio','.00101'),
                          ('carrier_row_max_ratio','nan'), ('current_total_A_per_um','nan'),
                          ('solver_method','probe'), ('bias_V','2.5')]:
            rejected = copy.deepcopy(row); rejected[key] = value
            self.assertFalse(idl.qualify(rejected,2.6), key)


if __name__ == '__main__':
    unittest.main()
