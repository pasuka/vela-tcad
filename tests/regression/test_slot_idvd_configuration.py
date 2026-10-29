"""IdVd physical package and qualification boundaries, without native data."""
import json
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]

class IdVdContract(unittest.TestCase):
    def setUp(self):
        self.c=json.loads((ROOT/'configs/regression/slot_ldmos_idvd_fermi_control.json').read_text())

    def test_si_and_poly_keep_distinct_physics(self):
        solver=self.c['solver']
        self.assertEqual(solver['carrier_statistics'],'fermi_dirac')
        self.assertTrue(solver['bandgap_narrowing']['fermi_statistics_correction'])
        self.assertEqual(solver['mobility']['constant_materials'],['PolySilicon'])
        self.assertEqual(solver['recombination_materials'],['Si'])
        self.assertEqual(solver['impact_ionization']['current_approximation'],'phumob_element_density_grad_qf')
        self.assertEqual(solver['impact_ionization']['source_jacobian'],'local_ad')
        self.assertNotIn('ialmob',solver['mobility'])

    def test_new_cold_lineage_preserves_acceptance(self):
        self.assertNotIn('initial_state_file',self.c['sweep'])
        self.assertTrue(all(c['bias']==0. for c in self.c['contacts']))
        s=self.c['solver']
        self.assertEqual(s['carrier_row_convergence']['mode'],'enforce')
        self.assertEqual(s['carrier_row_convergence']['eps_row'],.001)
        self.assertEqual(s['global_continuity_closure']['mode'],'enforce')
        self.assertEqual(s['global_continuity_closure']['tolerance'],.01)
        self.assertEqual(s['abstol'],1e-9)
        self.assertEqual(s['reltol'],1e-8)

    def test_sg_source_is_opt_in_and_preserves_the_aligned_system(self):
        base=json.loads((ROOT/'configs/regression/slot_ldmos_idvd_average_box_aligned_control.json').read_text())
        new=json.loads((ROOT/'configs/regression/slot_ldmos_idvd_sg_source_control.json').read_text())
        old_solver=base['solver'];new_solver=new['solver']
        self.assertEqual(old_solver['impact_ionization']['current_approximation'],'phumob_element_density_grad_qf')
        self.assertEqual(new_solver['impact_ionization']['current_approximation'],'phumob_element_sg_best_vertex')
        new_solver['impact_ionization']['current_approximation']=old_solver['impact_ionization']['current_approximation']
        # Paths may differ, but final physical/qualification parameters may not.
        def normalize(v):
            if isinstance(v,str):return v.replace('build/s10sourcealign/cold_template','build/s09alignment/cold')
            if isinstance(v,dict):return {k:normalize(x) for k,x in v.items()}
            if isinstance(v,list):return [normalize(x) for x in v]
            return v
        self.assertEqual(normalize(new_solver),old_solver)
        for key in ['mesh_geometry','discretization','materials_file','mesh_file','node_doping_file','contacts']:
            self.assertEqual(new[key],base[key])
        self.assertNotIn('initial_state_file',new['sweep'])

    def test_material_file_is_independent_and_not_pre_narrowed(self):
        path=ROOT/self.c['materials_file']
        self.assertNotEqual(path.name,'slot_ldmos_idl_bennett_materials.json')
        materials={m['name']:m for m in json.loads(path.read_text())['materials']}
        si,poly=materials['Si'],materials['PolySilicon']
        self.assertGreater(si['bandgap_eV'],1.12)
        self.assertGreater(poly['bandgap_eV'],1.12)
        self.assertEqual(poly['mun'],1417.)
        self.assertEqual(poly['mup'],470.5)
        self.assertNotEqual(poly['intrinsic_reference_potential_V'],0.)

if __name__=='__main__': unittest.main()
