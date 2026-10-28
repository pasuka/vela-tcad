import unittest
from scripts.vela_sdevice.deck_emit.bvmethods import apply_bvmethods_profile
from tests.regression.test_sdevice_deck_emit import qualification, material_result, emit
from vela_sdevice.cmd_parser import parse_cmd_text
from vela_sdevice import schema_validator as sv
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


class BVMethodsCmdProfileTests(unittest.TestCase):
    def test_joint_density_stop_is_valid_through_real_sidecar_schema(self):
        schema = sv.load_schema(REPO / "schemas/vela.sdevice_sidecar.v1.schema.json")
        sidecar = {"schema": "vela.sdevice_sidecar.v1", "solver": {
            "state_variables": "potential_density", "population_derivative_recovery": "automatic"},
            "sweep": {"external_circuit": {"coupled_update_convergence": {"mode": "enforce_joint"}}}}
        sv.validate(schema, sidecar)
        sidecar["sweep"]["external_circuit"]["coupled_update_convergence"]["mode"] = "invalid"
        with self.assertRaises(sv.ValidationError):
            sv.validate(schema, sidecar)

    def test_population_mode_is_explicit_and_emits_independent_protection(self):
        d, ir, stage, _ = self.inputs()
        sidecar = {"solver": {"state_variables": "potential_density",
            "population_update": {"protection": "local_nonpositive", "local_retained_fraction": .5}}}
        sidecar["solver"]["population_derivative_recovery"] = "automatic"
        metadata = apply_bvmethods_profile(d, ir, stage, sidecar)
        self.assertEqual(metadata["coordinates"], "potential_density")
        self.assertEqual(d["solver"]["population_derivative_recovery"], "automatic")
        self.assertEqual(d["solver"]["state_variables"], "potential_density")
        self.assertEqual(d["solver"]["population_update"]["local_retained_fraction"], .5)
        d["solver"]["mobility"]["doping_concentration_basis"] = "total_impurity"
        schema = sv.load_schema(REPO / "configs/schema/vela-simulation.schema.json")
        sv.validate({**schema, "required": []}, {"solver": d["solver"]})
        d["solver"]["population_update"]["local_retained_fraction"] = 0.
        with self.assertRaises(sv.ValidationError):
            sv.validate({**schema, "required": []}, {"solver": d["solver"]})

    def test_population_mode_rejects_unknown_coordinates_and_qf_caps(self):
        for values in ({"state_variables": "bad"},
                       {"population_derivative_recovery": "automatic"},
                       {"state_variables": "potential_density", "method": "gummel"},
                       {"state_variables": "potential_density", "quasi_fermi_update_limit_V": .1}):
            d, ir, stage, _ = self.inputs()
            with self.assertRaises(ValueError):
                apply_bvmethods_profile(d, ir, stage, {"solver": values})

    def test_population_only_controls_and_gummel_are_rejected_before_execution(self):
        schema = sv.load_schema(REPO / "configs/schema/vela-simulation.schema.json")
        schema = {**schema, "required": []}
        for solver in ({"population_update": {"protection": "local_nonpositive"}},
                       {"population_derivative_recovery": "automatic"}):
            with self.assertRaises(sv.ValidationError):
                sv.validate(schema, {"solver": solver})
        d, ir, stage, _ = self.inputs()
        apply_bvmethods_profile(d, ir, stage, {"solver": {"state_variables": "potential_density"}})
        d["solver"]["mobility"]["doping_concentration_basis"] = "total_impurity"
        d["solver"]["method"] = "gummel"
        with self.assertRaises(sv.ValidationError):
            sv.validate(schema, {"solver": d["solver"]})

    def inputs(self):
        deck = {'solver': {'mobility': {'model': 'masetti_field_lombardi',
            'high_field_driving_force': 'quasi_fermi_gradient', 'electron_cr_cm3': 9.68e16},
            'impact_ionization': {'model':'van_overstraeten','driving_force':'quasi_fermi_gradient',
                                  'electron_a_low_cm_inv':703000.},
            'carrier_statistics': {'model':'fermi_dirac'}, 'recombination':['srh'],
            'bandgap_narrowing': {'model':'old_slotboom'}, 'temperature_K':300.,
            'method':'newton','max_iter':40}, 'mesh_geometry':{},
            'sweep': {'start':0.,'stop':100000.,'step':2500.,'initial_step':.01,
                'min_step':1e-5,'max_step':2500.,'growth_factor':1.41,'shrink_factor':.5,
                'vtk_fields':['ElectrostaticPotential','eDensity','hDensity'],
                'external_circuit':{'resistance_ohm_um':1e7},
                'boundary_control':{'adaptive_outer_continuation':True}}}
        return deck, {'math':{'aval_derivatives':True}}, {'required_t':[], 'sweep':{'mode':'voltage','step_V':2500.}}, {}

    def test_profile_retains_coefficients_and_cmd_step_limits(self):
        d,ir,stage,sidecar=self.inputs();apply_bvmethods_profile(d,ir,stage,sidecar)
        self.assertEqual(d['solver']['mobility']['electron_cr_cm3'],9.68e16)
        self.assertEqual(d['solver']['impact_ionization']['electron_a_low_cm_inv'],703000.)
        self.assertEqual(d['solver']['carrier_statistics']['model'],'fermi_dirac_joyce_dixon')
        self.assertEqual(d['sweep']['external_circuit']['resistance_ohm_um'],1e7)
        self.assertEqual(d['sweep']['external_circuit']['coupled_max_outer_step_V'],2500.)
        self.assertEqual(d['sweep']['external_circuit']['coupled_initial_outer_step_V'],.01)
        self.assertEqual(d['sweep']['step'],100000.)
        self.assertEqual(stage['sweep']['step_V'],100000.)
        self.assertNotIn('boundary_control',d['sweep'])
        self.assertEqual(d['solver']['quasi_fermi_update_limit_V'],0.)
        self.assertEqual(d['solver']['quasi_fermi_update_limit_minority_V'],0.)

    def test_implicit_legacy_cap_is_removed_but_explicit_control_is_preserved(self):
        d,ir,stage,s=self.inputs();d['solver']['quasi_fermi_update_limit_V']=.1
        apply_bvmethods_profile(d,ir,stage,s)
        self.assertEqual(d['solver']['quasi_fermi_update_limit_V'],0.)
        d,ir,stage,s=self.inputs()
        s={'solver':{'quasi_fermi_update_limit_V':.1,'quasi_fermi_update_limit_minority_V':.05}}
        apply_bvmethods_profile(d,ir,stage,s)
        self.assertEqual(d['solver']['quasi_fermi_update_limit_V'],.1)
        self.assertEqual(d['solver']['quasi_fermi_update_limit_minority_V'],.05)

    def test_explicit_output_grid_is_retained(self):
        d,ir,stage,s=self.inputs();stage['required_t']=[.5]
        apply_bvmethods_profile(d,ir,stage,s)
        self.assertEqual(d['sweep']['step'],2500.)

    def test_sentaurus_undamped_default_and_local_override(self):
        d,ir,stage,s=self.inputs();apply_bvmethods_profile(d,ir,stage,s)
        self.assertEqual(d['sweep']['external_circuit']['coupled_line_search_mode'],'full_step')
        for option,value in [('line_search_damping',.1),('not_damped',0)]:
            d,ir,stage,s=self.inputs();ir['math'][option]=value
            with self.assertRaisesRegex(ValueError,'undamped Sentaurus'):
                apply_bvmethods_profile(d,ir,stage,s)
            # A local Coupled setting overrides the global Math setting.
            d,ir,stage,s=self.inputs();ir['math'][option]=value
            apply_bvmethods_profile(d,ir,stage,s,
                source_stage={'line_search_damping':1.,'not_damped':100})
            self.assertEqual(d['sweep']['external_circuit']['coupled_line_search_mode'],'full_step')

    def test_damping_values_survive_math_and_local_coupled_parsing(self):
        ir=parse_cmd_text('''Math {LineSearchDamping=.1 NotDamped=5}
            Solve {Coupled(LineSearchDamping=1 NotDamped=100 Iterations=20){Poisson Electron Hole}
            Quasistationary(Goal{Name="drain" Voltage=1}){
              Coupled(LineSearchDamping=1 NotDamped=200){Poisson Electron Hole}}}''',path_name='d.cmd')
        self.assertEqual(ir['math']['line_search_damping'],.1)
        self.assertEqual(ir['math']['not_damped'],5)
        self.assertEqual(ir['stages'][0]['not_damped'],100)
        self.assertEqual(ir['stages'][1]['not_damped'],200)
        self.assertEqual(ir['stages'][1]['line_search_damping'],1.)

    def test_profile_does_not_invent_a_series_resistor_for_equilibrium(self):
        d,ir,stage,s=self.inputs();del d['sweep']['external_circuit']
        apply_bvmethods_profile(d,ir,stage,s)
        self.assertNotIn('external_circuit',d['sweep'])

    def test_incompatible_physics_or_outputs_are_rejected(self):
        for change in ['boltzmann','btbt','no_av_derivatives','current_plot','auger','coordinate']:
            with self.subTest(change=change):
                d,ir,stage,s=self.inputs()
                if change=='boltzmann':d['solver']['carrier_statistics']['model']='boltzmann'
                elif change=='btbt':d['solver']['band_to_band']={'model':'e2'}
                elif change=='no_av_derivatives':ir['math']['aval_derivatives']=False
                elif change=='current_plot':d['sweep']['vtk_fields']=['eCurrentDensity/Vector']
                elif change=='coordinate':s={'solver':{'quasi_fermi_reference':'contact_majority'}}
                else:d['solver']['recombination'].append('auger')
                with self.assertRaises(ValueError):apply_bvmethods_profile(d,ir,stage,s)

    def test_cmd_stage_chain_and_sidecar_validate_together(self):
        text = '''File {Grid="device.tdr" Parameter="models.par"}
          Electrode { {Name="drain" Voltage=0 Resistor=1e7}
            {Name="source" Voltage=0} {Name="gate" Voltage=0 Barrier=-.55}
            {Name="substrate" Voltage=0} }
          Physics { EffectiveIntrinsicDensity(OldSlotboom) Fermi
            Mobility(DopingDep HighFieldSaturation(GradQuasiFermi) Enormal)
            Recombination(SRH(DopingDep) Avalanche(GradQuasiFermi)) }
          Plot { Potential eDensity hDensity SpaceCharge Doping }
          Math {Extrapolate AvalDerivatives Iterations=20
            BreakCriteria {Current(Contact="drain" AbsVal=1.3e-3)} }
          Solve { Coupled(Iterations=100) {Poisson}
            Coupled {Poisson Electron Hole}
            Quasistationary(InitialStep=1e-7 MinStep=1e-10 MaxStep=.025 Increment=1.41
              Goal {Name="drain" Voltage=100000}) {Coupled {Poisson Electron Hole}} }'''
        ir = parse_cmd_text(text, path_name="device.cmd")
        self.assertTrue(ir['verdict']['accepted'], ir['errors'])
        sidecar = {'schema':'vela.sdevice_sidecar.v1',
            'bvmethods_operator_profile':'bvmethods_tri3_qf',
            'solver': {'max_iter':40, 'carrier_row_convergence':{'mode':'report'},
                       'global_continuity_closure':{'mode':'report'}},
            'sweep':{'external_circuit':{'solver':'coupled_newton',
                'coupled_line_search_mode':'report_only',
                'coupled_update_convergence':{'mode':'off','diagnostic_csv':'updates.csv'},
                'coupled_equation_tolerance':1e-8, 'coupled_max_step_retries':20,
                'coupled_continuation':{'enabled':True,'retry_from_accepted':True}}}}
        sv.validate(sv.load_schema(REPO/'schemas/vela.sdevice_sidecar.v1.schema.json'), sidecar)
        result = emit(ir, mesh_qualification=qualification(ir, ('Si','SiO2')),
                      materials_report=material_result(ir, ('Si','SiO2')).report, sidecar=sidecar)
        self.assertEqual(len(result.decks), 2)
        zero, scan = result.decks
        self.assertNotIn('external_circuit', zero['sweep'])
        self.assertEqual(zero['sweep']['initialization']['poisson_max_iter'],100)
        self.assertEqual(scan['sweep']['initial_state_file'],zero['sweep']['write_state_file'])
        self.assertEqual(scan['sweep']['external_circuit']['resistance_ohm_um'],1e7)
        self.assertEqual(scan['sweep']['external_circuit']['coupled_line_search_mode'],'report_only')
        self.assertEqual(scan['sweep']['external_circuit']['coupled_update_convergence'],
                         {'mode':'off','diagnostic_csv':'updates.csv'})
        self.assertEqual(scan['sweep']['external_circuit']['coupled_max_outer_step_V'],2500)
        self.assertEqual(scan['solver']['carrier_row_convergence'],{'mode':'report'})
        for deck in result.decks:
            self.assertEqual(deck['solver']['coupled_operator'],'bvmethods_tri3_qf')
            self.assertEqual(deck['solver']['quasi_fermi_update_limit_V'],0.)
            self.assertEqual(deck['solver']['quasi_fermi_update_limit_minority_V'],0.)
            sv.validate(sv.load_schema(REPO/'configs/schema/vela-simulation.schema.json'),deck)
        sv.validate(sv.load_schema(REPO/'schemas/vela.stage_plan.v1.schema.json'),result.stage_plan)

    def test_sidecar_cannot_override_cmd_outer_adaptation(self):
        schema = sv.load_schema(REPO/'schemas/vela.sdevice_sidecar.v1.schema.json')
        for key in ('resistance_ohm_um','coupled_max_outer_step_V','coupled_initial_outer_step_V'):
            self.assertTrue(sv.collect_errors(schema, {'schema':'vela.sdevice_sidecar.v1',
                'sweep':{'external_circuit':{key:1.}}}))


if __name__=='__main__':unittest.main()
