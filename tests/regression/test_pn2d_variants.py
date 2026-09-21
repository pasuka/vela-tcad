"""Prospective variant isolation, exact voltage targets, and input semantics."""
import importlib.util
import math
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import pn2d_variants as v
import compare_pn2d_variants as compare
from analyze_pn2d_sections import nodal_section


class Variants(unittest.TestCase):
    def setUp(self):
        self.spec = v.read_json(v.SPEC)

    def test_reconstructed_section_linear_field_units_and_shared_edges(self):
        # Divergence-free J=(2+3y, 4-2x) must have an x-independent x cut.
        nodes=[{'id':i,'x':x,'y':y} for i,(x,y) in enumerate([(0,0),(1,0),(2,0),(0,.5),(1,.5),(2,.5)])]
        mesh={'nodes':nodes,'triangles':[{'node_ids':t} for t in [[0,1,4],[0,4,3],[1,2,5],[1,5,4]]]}
        field=[(2+3*n['y'],4-2*n['x']) for n in nodes]
        for x in [0,.25,1,1.75,2]:
            self.assertAlmostEqual(nodal_section(mesh,field,'x',x),1.375e-8,delta=1e-22)
        for y in [0,.2,.5]:
            self.assertAlmostEqual(nodal_section(mesh,field,'y',y),4e-8,delta=1e-22)
        self.assertEqual(nodal_section(mesh,field,'x',3),0)
        self.assertAlmostEqual(nodal_section(mesh,field,'x',1,2),2.75e-8,delta=1e-22)

    def test_current_recovery_area_weighting_and_regional_error(self):
        from diagnose_pn2d_current_reconstruction import metrics
        # Same vector on unequal control volumes must have zero error and
        # preserve the physical transverse/longitudinal ratio.
        exact=metrics([(3.,4.),(3.,4.)],[(3.,4.),(3.,4.)],[.1,.9],[True,True])
        self.assertEqual(exact['relative_l2'],0.)
        self.assertAlmostEqual(exact['Jy_over_Jx_l2'],4/3)
        # A localized unit error on 10% of area has RMS sqrt(0.1), while
        # excluding that region removes the error without changing other data.
        whole=metrics([(1.,0.),(1.,0.)],[(2.,0.),(1.,0.)],[.1,.9],[True,True])
        outside=metrics([(1.,0.),(1.,0.)],[(2.,0.),(1.,0.)],[.1,.9],[False,True])
        self.assertAlmostEqual(whole['relative_l2'],math.sqrt(.1))
        self.assertAlmostEqual(whole['squared_error_integral'],.1)
        self.assertEqual(outside['relative_l2'],0.)
        self.assertAlmostEqual(outside['area_um2'],.9)

    def test_eight_unique_one_factor_experiments(self):
        rows = self.spec['experiments']
        self.assertEqual([r['id'] for r in rows], ['P0','P1','P2','P3','D1','D2','G1','G2'])
        self.assertTrue(all(len(r['overrides']) == (0 if r['id']=='P0' else 1) for r in rows))
        self.assertEqual(self.spec['aliases'], {'D0':'P0','G0':'P0'})

    def test_fixed_exact_goals_and_independent_branches(self):
        cfg=self.spec['baseline']
        for branch, count in [('forward',41),('reverse',11)]:
            deck=v.vela_config(self.spec,cfg,branch)
            points=deck['sweep']['bias_points']
            self.assertEqual(len(points),count)
            self.assertEqual(points[0],0)
            self.assertNotIn('initial_state_file',deck['sweep'])
            source=v.sdevice_text(self.spec,cfg,branch)
            for bias in points[1:]:
                self.assertIn(f'Voltage={bias:.12g} ',source)
            self.assertNotIn('Intervals=',source)
            self.assertNotIn('Load(',source)

    def test_contact_segments_preserve_geometry_and_doping(self):
        for fraction, expected in [(1.0,5),(0.5,3),(0.25,1)]:
            cfg={**self.spec['baseline'],'anode_fraction':fraction}
            source=v.sde_text(self.spec,cfg,self.spec['meshes']['M0'])
            self.assertEqual(source.count('(sdegeo:insert-vertex'),4)
            self.assertEqual(source.count('(sdegeo:define-2d-contact'),expected+1)
            self.assertIn('(define XJ 1.0)',source)
            self.assertIn('"BoronActiveConcentration"\n  1e+17',source)

    def test_doping_reaches_sde_not_region_average(self):
        source=v.sde_text(self.spec,{**self.spec['baseline'],'NA_cm3':1e18},self.spec['meshes']['M0'])
        self.assertIn('"BoronActiveConcentration"\n  1e+18',source)
        self.assertIn('"PhosphorusActiveConcentration"\n  1e+17',source)

    def test_bgn_includes_gap_reference_convention(self):
        on=v.material(True)['materials'][0]
        off=v.material(False)['materials'][0]
        self.assertAlmostEqual(on['ni']/off['ni'],math.exp(0.01595/(2*8.617333262145e-5*300)))
        self.assertEqual(on['mun'],1417)
        self.assertEqual(on['mup'],470.5)

    def test_generation_is_reproducible_and_paths_portable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            a=v.generate(root)
            b=v.generate(root)
            self.assertEqual(a,b)
            self.assertEqual(len(a),1+8*(len(self.spec['meshes'])-1))
            self.assertTrue(all('\\' not in e['directory'] for e in a))
            self.assertEqual((root/'P0/M0/pn2d_sde.cmd').read_bytes(),(root/'P2/M0/pn2d_sde.cmd').read_bytes())

    def test_manifest_cannot_silently_change_unsupported_physics(self):
        self.spec['temperature_K']=350
        with self.assertRaises(ValueError): v.validate_spec(self.spec)

    def test_exact_points_cannot_be_interpolated(self):
        data=[{'bias_V':0.0},{'bias_V':0.04}]
        with self.assertRaises(ValueError): compare.exact_row(data,0.02)
        self.assertEqual(compare.exact_row(data,0.04)['bias_V'],0.04)

    def test_partial_curve_is_not_comparison_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'live.csv'
            path.write_text('bias_V,current\n0.02\n')
            with self.assertRaises(ValueError): compare.rows(path)

    def test_cached_comparison_invalidates_changed_imported_doping(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);v.generate(root)
            before=compare.compare(root)['results'][1]['evidence_fingerprint']
            path=root/'P0/M0/inputs/doping.csv';path.parent.mkdir(exist_ok=True)
            path.write_text('node_id,donors_cm3,acceptors_cm3\n0,0,1e18\n')
            after=compare.compare(root)['results'][1]['evidence_fingerprint']
            self.assertNotEqual(before,after)

    def test_current_floor_sign_and_nonfinite(self):
        gate=v.read_json(v.CONTRACT)['current']
        self.assertTrue(compare.compare_current(0,1e-25,gate)['pass'])
        self.assertFalse(compare.compare_current(1e-10,-1e-10,gate)['pass'])
        self.assertFalse(compare.compare_current(1e-18,1e-17,gate)['pass'])
        with self.assertRaises(ValueError): compare.compare_current(1,float('nan'),gate)

    def test_node_matching_rejects_duplicate_and_missing_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'field.csv'
            path.write_text('node_id,component0\n0,1\n0,2\n')
            with self.assertRaises(ValueError): compare.node_field(path,[0,1])

    def test_mesh_and_stage_gates_fail_closed(self):
        self.assertEqual(compare.mesh_gate('P0',{},v.read_json(v.CONTRACT))['status'],'not_run')
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);v.generate(root)
            v.write_json(root/'results.json',{'stage_A_accepted':False})
            with self.assertRaises(ValueError): v.selected(root,'D1','M0')

    def test_stage_gate_does_not_trust_an_unbacked_success_flag(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);v.generate(root)
            v.write_json(root/'results.json',{'stage_A_accepted':True,'contract_sha256':v.sha(v.CONTRACT),'results':[]})
            with self.assertRaises(ValueError): v.selected(root,'D1','M2')

    def test_vtk_matching_rejects_equal_count_reordered_nodes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'state.vtk'
            path.write_text('POINTS 2 double\n0 0 0\n1 0 0\n')
            compare.validate_vtk_coordinates(path,[(0,0),(1,0)],1e-10)
            with self.assertRaises(ValueError): compare.validate_vtk_coordinates(path,[(1,0),(0,0)],1e-10)

    def test_stage_c_requires_completed_stage_b(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);v.generate(root)
            v.write_json(root/'results.json',{'stage_A_accepted':True,'contract_sha256':v.sha(v.CONTRACT),
                         'results':[],'refined_mesh':{'D1':{'status':'not_run'},'D2':{'status':'not_run'}}})
            with self.assertRaisesRegex(ValueError,'Stage B'):
                v.selected(root,'G1','M2')
            self.assertEqual(len(v.selected(root,'G1','M2',allow_modeling=True)),1)

    def test_refinement_preserves_thresholds_and_coarse_inputs(self):
        continuation=v.read_json(v.SPEC.parent/'mesh_refinement.json')
        self.assertEqual(continuation['frozen_contract_sha256'],v.sha(v.CONTRACT))
        self.assertEqual(continuation['original_pair'],v.read_json(v.CONTRACT)['mesh']['required_pair'])
        self.assertEqual(continuation['additional_pair'],['M2','M3'])
        self.assertIn(['M1','M2'],continuation['previous_additional_pairs'])
        a=v.sde_text(self.spec,self.spec['baseline'],self.spec['meshes']['M1'])
        b=v.sde_text(self.spec,self.spec['baseline'],self.spec['meshes']['M2'])
        self.assertEqual(a.replace('"Junction.Mesh"\n  0.005 0.01\n  0.0025 0.0025','"Junction.Mesh"\n  0.0025 0.01\n  0.00125 0.0025'),b)

    def test_m3_uses_solved_goals_and_only_required_spatial_outputs(self):
        for branch in ['forward','reverse']:
            deck=v.sdevice_text(self.spec,self.spec['baseline'],branch,self.spec['meshes']['M3'])
            for bias in v.biases(self.spec['branches'][branch])[1:]:
                self.assertIn(f'Voltage={bias:.12g} ',deck)
            self.assertIn('InitialStep=1 ',deck)
            self.assertEqual(deck.count('Plot(FilePrefix='),4 if branch=='forward' else 2)
            self.assertIn('Digits=8',deck)

    def test_native_precision_control_preserves_physics_and_stage_a(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);v.generate(root)
            for identifier in ['P0','P1','P2','P3']:
                self.assertNotIn('ExtendedPrecision',(root/identifier/'M2/forward_sdevice.cmd').read_text())
            for identifier in ['D1','D2','G1','G2']:
                for branch,count in [('forward',40),('reverse',10)]:
                    text=(root/identifier/'M2'/f'{branch}_sdevice.cmd').read_text()
                    self.assertIn('ExtendedPrecision(80)',text)
                    self.assertIn('RHSMin=1e-20',text)
                    self.assertIn('Digits=8',text)
                    self.assertEqual(text.count('InitialStep=1 '),count)
                    self.assertIn('Mobility(DopingDependence)',text)
                    self.assertIn('Recombination(SRH)',text)
                    self.assertIn('EffectiveIntrinsicDensity(OldSlotboom)',text)
                    self.assertNotIn('Load(',text)

    def test_conditioning_does_not_relax_convergence(self):
        deck=v.vela_config(self.spec,self.spec['baseline'],'forward')
        self.assertEqual(deck['solver']['continuity_row_scaling']['flux_fraction'],1.0)
        self.assertEqual(deck['solver']['global_continuity_closure']['tolerance'],0.01)
        self.assertEqual(deck['solver']['carrier_row_convergence']['eps_row'],0.001)

    def test_contact_regions_follow_actual_electrode_endpoints(self):
        mesh={'nodes':[{'id':i,'x':0.0,'y':y} for i,y in enumerate([0,.125,.1875,.25,.3125,.375,.5])],
              'contacts':[{'name':'Anode','node_ids':[1,2,3,4,5]}]}
        mask=compare.spatial_regions(mesh,v.read_json(v.CONTRACT))['contact_edges']
        self.assertEqual(mask,[False,True,False,False,False,True,False])
        mesh['contacts'][0]['node_ids']=[2,3,4]
        self.assertEqual(compare.spatial_regions(mesh,v.read_json(v.CONTRACT))['contact_edges'],
                         [False,False,True,False,True,False,False])

    def test_endpoint_refinement_changes_only_endpoint_size_lines(self):
        coarse=v.sde_text(self.spec,self.spec['baseline'],self.spec['meshes']['M3']).splitlines()
        fine=v.sde_text(self.spec,self.spec['baseline'],self.spec['meshes']['E1']).splitlines()
        changed=[(a,b) for a,b in zip(coarse,fine,strict=True) if a!=b]
        self.assertEqual(len(changed),4)
        self.assertTrue(all('define-refinement-size "EdgeMesh' in a and 'define-refinement-size "EdgeMesh' in b for a,b in changed))
        contract=v.read_json(v.SPEC.parent/'mesh_refinement.json')['contact_edge_control']
        self.assertEqual(contract['ids'],['P0','G1','G2'])
        self.assertTrue(contract['declared_before_G_results'])

    def test_conservation_requires_numeric_closure_not_only_success_flag(self):
        gate=v.read_json(v.CONTRACT)
        state={'global_continuity_closure_satisfied':'true',
               'global_electron_integrated_source':'1','global_hole_integrated_source':'1',
               'global_electron_continuity_closure_ratio':'0.001',
               'global_hole_continuity_closure_ratio':'0.002'}
        self.assertTrue(compare.continuity_satisfied(state,gate))
        for value in ['0.011','nan','-0.001']:
            state['global_hole_continuity_closure_ratio']=value
            self.assertFalse(compare.continuity_satisfied(state,gate))

    def test_zero_source_checks_through_current_balance(self):
        gate=v.read_json(v.CONTRACT)
        state={'global_continuity_closure_satisfied':'true',
               'global_electron_integrated_source':'0','global_hole_integrated_source':'0',
               'global_electron_continuity_closure_ratio':'1','global_hole_continuity_closure_ratio':'1'}
        ports={name:{'current_electron_A_per_um':value,'current_hole_A_per_um':-value}
               for name,value in [('Anode',1e-10),('Cathode',-1e-10+1e-22)]}
        self.assertTrue(compare.continuity_satisfied(state,gate,ports))
        ports['Cathode']['current_electron_A_per_um']=-.98e-10
        self.assertFalse(compare.continuity_satisfied(state,gate,ports))


if __name__=='__main__': unittest.main()
