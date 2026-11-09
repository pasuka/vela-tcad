import unittest
from scripts.vela_sdevice.deck_emit.bvmethods import apply_bvmethods_profile
from tests.regression.test_sdevice_deck_emit import qualification, material_result, emit
from vela_sdevice.cmd_parser import parse_cmd_text
from vela_sdevice import schema_validator as sv
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


class BVMethodsCmdProfileTests(unittest.TestCase):
    def test_approximate_stopping_reads_math_update_references_and_preserves_override(self):
        d, ir, stage, _ = self.inputs()
        ir["math"].update(digits=6, err_ref_electron_cm3=2e10, err_ref_hole_cm3=3e10)
        d["sweep"]["external_circuit"]["coupled_update_convergence"]={"digits": 4}
        apply_bvmethods_profile(d, ir, stage, {})
        metric=d["sweep"]["external_circuit"]["coupled_update_convergence"]
        self.assertEqual(metric, {"mode": "report", "digits": 4,
            "electron_reference_cm3": 2e10, "hole_reference_cm3": 3e10})
        d, ir, stage, _ = self.inputs()
        ir["math"]["rel_err_control"]=False
        with self.assertRaisesRegex(ValueError, "RelErrControl"):
            apply_bvmethods_profile(d, ir, stage, {})

    def test_load_report_and_or_stopping_are_independent_schema_controls(self):
        for mode, load in [("enforce", "enforce"), ("enforce_joint", "report"), ("enforce", "report")]:
            controls = {"solver": "coupled_newton",
                "load_line_convergence_mode": load, "coupled_update_convergence": {"mode": mode}}
            for name, value in [("schemas/vela.sdevice_sidecar.v1.schema.json",
                                 {"schema": "vela.sdevice_sidecar.v1", "sweep": {"external_circuit": controls}}),
                                ("configs/schema/vela-simulation.schema.json", {"sweep": {"external_circuit": controls}})]:
                schema = sv.load_schema(REPO / name)["properties"]["sweep"]["properties"]["external_circuit"]
                schema = {**schema, "required": []}
                sv.validate(schema, controls)
                controls["load_line_convergence_mode"] = "bad"
                with self.assertRaises(sv.ValidationError):sv.validate(schema, controls)
                controls["load_line_convergence_mode"] = load

    def test_joint_density_stop_is_valid_through_real_sidecar_schema(self):
        schema = sv.load_schema(REPO / "schemas/vela.sdevice_sidecar.v1.schema.json")
        sidecar = {"schema": "vela.sdevice_sidecar.v1", "solver": {
            "state_variables": "potential_density", "population_derivative_recovery": "automatic"},
            "sweep": {"external_circuit": {"coupled_update_convergence": {"mode": "enforce_joint"}}}}
        sv.validate(schema, sidecar)
        sidecar["sweep"]["external_circuit"]["coupled_update_convergence"]["mode"] = "invalid"
        with self.assertRaises(sv.ValidationError):
            sv.validate(schema, sidecar)

    def test_direct_population_sg_is_explicit_and_validated(self):
        d, ir, stage, _ = self.inputs()
        apply_bvmethods_profile(d, ir, stage, {"solver": {
            "state_variables": "potential_density", "population_sg_current": "band_population"}})
        self.assertEqual(d["solver"]["population_sg_current"], "band_population")
        d["solver"]["mobility"]["doping_concentration_basis"] = "total_impurity"
        schema = sv.load_schema(REPO / "configs/schema/vela-simulation.schema.json")
        schema = {**schema, "required": []}
        sv.validate(schema, {"solver": d["solver"]})
        for value in ["bad", "band_population"]:
            invalid = {"population_sg_current": value}
            with self.assertRaises(sv.ValidationError):sv.validate(schema, {"solver": invalid})
            d, ir, stage, _ = self.inputs()
            with self.assertRaises(ValueError):apply_bvmethods_profile(d, ir, stage, {"solver": {"state_variables": "potential_quasi_fermi", **invalid}})

    def test_statistical_bootstrap_schema_requires_compatible_population_floor(self):
        import copy
        schema = sv.load_schema(REPO / "configs/schema/vela-simulation.schema.json")
        schema = {**schema, "required": []}
        d, ir, stage, _ = self.inputs()
        apply_bvmethods_profile(d, ir, stage, {"solver": {"state_variables": "potential_density"}})
        solver = d["solver"]
        solver["mobility"]["doping_concentration_basis"] = "total_impurity"
        solver["population_sg_bootstrap_linearization"] = "floor_adjacent_edges"
        solver["population_update"] = {"minimum_density_cm3": 1e-100}
        sv.validate(schema, {"solver": solver})
        forwarded, ir, stage, _ = self.inputs()
        apply_bvmethods_profile(forwarded, ir, stage, {"solver": {
            "state_variables": "potential_density", "population_sg_bootstrap_linearization": "floor_adjacent_edges",
            "population_update": {"minimum_density_cm3": 1e-100}}})
        self.assertEqual(forwarded["solver"]["population_sg_bootstrap_linearization"], "floor_adjacent_edges")
        defaults, ir, stage, _ = self.inputs()
        apply_bvmethods_profile(defaults, ir, stage, {"solver": {
            "population_sg_bootstrap_linearization": "floor_adjacent_edges"}})
        self.assertEqual(defaults["solver"]["population_update"]["minimum_density_cm3"], 1e-100)
        for control in [{"population_sg_bootstrap_linearization": "bad"}]:
            bad, ir, stage, _ = self.inputs()
            with self.assertRaises(ValueError):apply_bvmethods_profile(bad, ir, stage, {"solver": control})
        for key, value in [("population_update", {"minimum_density_cm3": 0}),
                           ("population_sg_current", "band_population"),
                           ("population_sg_bootstrap_linearization", "bad")]:
            bad = copy.deepcopy(solver); bad[key] = value
            with self.assertRaises(sv.ValidationError):sv.validate(schema, {"solver": bad})

    def test_population_floor_unit_alias_does_not_emit_two_floor_fields(self):
        d, ir, stage, _ = self.inputs()
        apply_bvmethods_profile(d, ir, stage, {"solver": {
            "population_update": {"minimum_density_m3": 1e-100}}})
        self.assertNotIn("minimum_density_cm3", d["solver"]["population_update"])
        self.assertEqual(d["solver"]["population_update"]["minimum_density_m3"], 1e-100)

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
                       {"state_variables": "potential_quasi_fermi", "population_derivative_recovery": "automatic"},
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

    def test_population_deck_removes_inherited_qf_controls_but_preserves_physics(self):
        d, ir, stage, _ = self.inputs()
        d["solver"].update(quasi_fermi_update_limit_V=.1,
            quasi_fermi_update_limit_mode="uniform_trust_region",
            quasi_fermi_trust_region_growth_factor=2.)
        apply_bvmethods_profile(d, ir, stage, {"solver": {"state_variables": "potential_density"}})
        for key in ("quasi_fermi_update_limit_V", "quasi_fermi_update_limit_minority_V",
                    "quasi_fermi_reference", "quasi_fermi_update_limit_mode", "quasi_fermi_trust_region_growth_factor"):
            self.assertNotIn(key, d["solver"])
        self.assertEqual(d["solver"]["mobility"]["high_field_driving_force"], "quasi_fermi_gradient")
        self.assertEqual(d["solver"]["impact_ionization"]["driving_force"], "quasi_fermi_gradient")

    def test_population_qf_policies_fail_in_cmd_and_json_before_execution(self):
        schema = {**sv.load_schema(REPO / "configs/schema/vela-simulation.schema.json"), "required": []}
        for option in ({"quasi_fermi_recenter_on_stall": True},
                       {"quasi_fermi_trust_region_growth_factor": 2.},
                       {"contact_majority_qf_branch_drop_limit_V": 1e-9},
                       {"stable_merit_comparison": True}):
            d, ir, stage, _ = self.inputs()
            with self.assertRaises(ValueError):
                apply_bvmethods_profile(d, ir, stage, {"solver": {"state_variables": "potential_density", **option}})
            with self.assertRaises(sv.ValidationError):
                sv.validate(schema, {"solver": {"state_variables": "potential_density",
                    "coupled_operator": "bvmethods_tri3_qf", "carrier_statistics": {"model": "fermi_dirac_joyce_dixon"}, **option}})

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
        self.assertEqual(d['sweep']['external_circuit']['coupled_step_control'],'quasistationary')
        self.assertNotIn('coupled_continuation',d['sweep']['external_circuit'])
        self.assertNotIn('max_inner_voltage_step_V',d['sweep']['external_circuit'])
        self.assertEqual(d['sweep']['step'],100000.)
        self.assertEqual(stage['sweep']['step_V'],100000.)
        self.assertNotIn('boundary_control',d['sweep'])
        self.assertNotIn('quasi_fermi_update_limit_V',d['solver'])
        self.assertNotIn('quasi_fermi_update_limit_minority_V',d['solver'])

    def test_implicit_legacy_cap_is_removed_but_explicit_control_is_preserved(self):
        d,ir,stage,s=self.inputs();d['solver']['quasi_fermi_update_limit_V']=.1
        apply_bvmethods_profile(d,ir,stage,s)
        self.assertNotIn('quasi_fermi_update_limit_V',d['solver'])
        d,ir,stage,s=self.inputs()
        s={'solver':{'state_variables':'potential_quasi_fermi','quasi_fermi_update_limit_V':.1,'quasi_fermi_update_limit_minority_V':.05}}
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

    def test_default_population_profile_and_explicit_qf_rollback(self):
        d, ir, stage, _ = self.inputs()
        apply_bvmethods_profile(d, ir, stage, {})
        self.assertEqual(d["solver"]["state_variables"], "potential_density")
        self.assertEqual(d["solver"]["primary_variable_coordinates"], "unit")
        self.assertEqual(d["solver"]["population_derivative_recovery"], "after_protection")
        self.assertEqual(d["solver"]["population_stopping_rule"], "native_approx")
        legacy, ir, stage, _ = self.inputs()
        apply_bvmethods_profile(legacy, ir, stage, {"solver": {"state_variables": "potential_quasi_fermi"}})
        self.assertNotIn("population_stopping_rule", legacy["solver"])
        self.assertEqual(legacy["solver"]["quasi_fermi_update_limit_V"], 0.)

    def test_post_protection_policy_survives_profile_and_schema(self):
        d, ir, stage, _ = self.inputs()
        apply_bvmethods_profile(d, ir, stage, {"solver": {
            "state_variables": "potential_density",
            "population_derivative_recovery": "after_protection"}})
        self.assertEqual(d["solver"]["population_derivative_recovery"], "after_protection")
        d["solver"]["mobility"]["doping_concentration_basis"] = "total_impurity"
        schema = sv.load_schema(REPO / "configs/schema/vela-simulation.schema.json")
        sv.validate({**schema, "required": []}, {"solver": d["solver"]})
        with self.assertRaises(sv.ValidationError):
            sv.validate({**schema, "required": []}, {"solver": {
                "state_variables": "potential_quasi_fermi",
                "population_derivative_recovery": "after_protection"}})

    def test_explicit_basic_state_plot_projection_is_schema_bounded(self):
        schema = sv.load_schema(REPO / "schemas/vela.sdevice_sidecar.v1.schema.json")
        valid = {"schema": "vela.sdevice_sidecar.v1", "sweep": {
            "vtk_fields": ["ElectrostaticPotential", "eDensity", "hDensity"]}}
        sv.validate(schema, valid)
        for bad in ([], ["UnqualifiedField"], ["eDensity", "eDensity"]):
            with self.subTest(fields=bad), self.assertRaises(sv.ValidationError):
                sv.validate(schema, {**valid, "sweep": {"vtk_fields": bad}})
        with self.assertRaises(sv.ValidationError):
            sv.validate(schema, {**valid, "sweep": {"vtk_fields": ["eDensity"], "stop": 1}})

    def test_implicit_population_state_has_the_same_schema_contract(self):
        d, ir, stage, _ = self.inputs()
        apply_bvmethods_profile(d, ir, stage, {})
        solver = d["solver"]
        solver["mobility"]["doping_concentration_basis"] = "total_impurity"
        del solver["state_variables"]
        solver["population_sg_bootstrap_linearization"] = "floor_adjacent_edges"
        schema = {**sv.load_schema(REPO / "configs/schema/vela-simulation.schema.json"),
                  "required": []}
        sv.validate(schema, {"solver": solver})
        for changed in ({"state_variables": "potential_quasi_fermi"},
                        {"coupled_operator": "legacy_edge"},
                        {"quasi_fermi_trust_region_growth_factor": 2.}):
            with self.subTest(changed=changed), self.assertRaises(sv.ValidationError):
                sv.validate(schema, {"solver": {**solver, **changed}})

    def test_resistor_budget_uses_sweep_coupled_not_seed_budget(self):
        for local, math, explicit, expected in [(20, 40, None, 20),
                                                (None, 20, None, 20),
                                                (20, 40, 12, 12)]:
            with self.subTest(local=local, math=math, explicit=explicit):
                d, ir, stage, sidecar = self.inputs()
                d["solver"]["max_iter"] = 100
                ir["math"]["iterations"] = math
                if explicit is not None:
                    d["sweep"]["external_circuit"]["max_iterations"] = explicit
                apply_bvmethods_profile(d, ir, stage, sidecar,
                    source_stage={"iterations": local, "not_damped": 100})
                self.assertEqual(d["sweep"]["external_circuit"]["max_iterations"], expected)

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
                'load_line_convergence_mode':'report',
                'coupled_update_convergence':{'mode':'off','diagnostic_csv':'updates.csv'},
                'coupled_equation_tolerance':1e-8, 'coupled_max_step_retries':20,
                'coupled_step_control':'legacy_inner_guard',
                'coupled_continuation':{'enabled':True,'retry_from_accepted':True}}}}
        sv.validate(sv.load_schema(REPO/'schemas/vela.sdevice_sidecar.v1.schema.json'), sidecar)
        result = emit(ir, mesh_qualification=qualification(ir, ('Si','SiO2')),
                      materials_report=material_result(ir, ('Si','SiO2')).report, sidecar=sidecar)
        self.assertEqual(len(result.decks), 2)
        # A deliberate output projection must be visible and leave the solve
        # contract intact. Unsupported physical diagnostics still reject by default.
        import copy
        wide_ir = parse_cmd_text(text.replace("Plot { Potential", "Plot { eCurrent/Vector Potential"),
                                 path_name="device.cmd")
        with self.assertRaises(ValueError):
            emit(wide_ir, mesh_qualification=qualification(wide_ir, ('Si', 'SiO2')),
                 materials_report=material_result(wide_ir, ('Si', 'SiO2')).report,
                 sidecar=sidecar)
        projected_sidecar = copy.deepcopy(sidecar)
        projected_sidecar['sweep']['vtk_fields'] = list(result.decks[0]['sweep']['vtk_fields'])
        projected = emit(wide_ir, mesh_qualification=qualification(wide_ir, ('Si', 'SiO2')),
                         materials_report=material_result(wide_ir, ('Si', 'SiO2')).report,
                         sidecar=projected_sidecar)
        for original, narrowed in zip(result.decks, projected.decks):
            self.assertEqual(original['solver'], narrowed['solver'])
            self.assertEqual(original['sweep'], narrowed['sweep'])
        self.assertTrue(any(d.get('message_id') == 'W207' and
                            d.get('vela_key') == 'sweep.vtk_fields' for d in projected.decisions))
        zero, scan = result.decks
        self.assertNotIn('external_circuit', zero['sweep'])
        self.assertEqual(zero['sweep']['initialization']['poisson_max_iter'],100)
        self.assertEqual(scan['sweep']['initial_state_file'],zero['sweep']['write_state_file'])
        self.assertEqual(scan['sweep']['external_circuit']['resistance_ohm_um'],1e7)
        self.assertEqual(scan['sweep']['external_circuit']['max_iterations'],20)
        for name, cmd, override, expected in (
            ("local", text.replace("{Coupled {Poisson Electron Hole}}",
                                   "{Coupled(Iterations=12) { Poisson Electron Hole}}"), None, 12),
            ("default", text.replace(" Iterations=20", ""), None, 50),
            ("override", text, 17, 17),
        ):
            with self.subTest(budget=name):
                case = parse_cmd_text(cmd, path_name="device.cmd")
                self.assertTrue(case["verdict"]["accepted"], case["errors"])
                circuit = dict(sidecar["sweep"]["external_circuit"])
                options = {**sidecar, "sweep": {"external_circuit": circuit}}
                if override is not None:
                    circuit["max_iterations"] = override
                scan_deck = emit(case, mesh_qualification=qualification(case, ("Si", "SiO2")),
                                 materials_report=material_result(case, ("Si", "SiO2")).report,
                                 sidecar=options).decks[1]
                self.assertEqual(scan_deck["sweep"]["external_circuit"]["max_iterations"], expected)
        self.assertEqual(scan['sweep']['external_circuit']['coupled_line_search_mode'],'report_only')
        self.assertEqual(scan['sweep']['external_circuit']['load_line_convergence_mode'],'report')
        self.assertEqual(scan['sweep']['external_circuit']['coupled_update_convergence'],
                         {'mode':'off','diagnostic_csv':'updates.csv','digits':5})
        self.assertEqual(scan['sweep']['external_circuit']['coupled_max_outer_step_V'],2500)
        self.assertEqual(scan['solver']['carrier_row_convergence'],{'mode':'report'})
        for deck in result.decks:
            self.assertEqual(deck['solver']['coupled_operator'],'bvmethods_tri3_qf')
            self.assertNotIn('quasi_fermi_update_limit_V',deck['solver'])
            self.assertNotIn('quasi_fermi_update_limit_minority_V',deck['solver'])
            sv.validate(sv.load_schema(REPO/'configs/schema/vela-simulation.schema.json'),deck)
        sv.validate(sv.load_schema(REPO/'schemas/vela.stage_plan.v1.schema.json'),result.stage_plan)

    def test_native_external_policy_rejects_custom_guards(self):
        for guard in ({'max_inner_voltage_step_V':.5},
                      {'coupled_continuation':{'enabled':True}},
                      {'coupled_continuation':{'retry_from_accepted':True}}):
            d,ir,stage,_=self.inputs()
            d['sweep']['external_circuit'].update(guard)
            with self.assertRaises(ValueError):apply_bvmethods_profile(d,ir,stage,{})
            d,ir,stage,_=self.inputs()
            d['sweep']['external_circuit'].update(guard)
            d['sweep']['external_circuit']['coupled_step_control']='legacy_inner_guard'
            apply_bvmethods_profile(d,ir,stage,{})
            self.assertEqual(d['sweep']['external_circuit']['coupled_step_control'],'legacy_inner_guard')

    def test_sidecar_cannot_override_cmd_outer_adaptation(self):
        schema = sv.load_schema(REPO/'schemas/vela.sdevice_sidecar.v1.schema.json')
        for key in ('resistance_ohm_um','coupled_max_outer_step_V','coupled_initial_outer_step_V'):
            self.assertTrue(sv.collect_errors(schema, {'schema':'vela.sdevice_sidecar.v1',
                'sweep':{'external_circuit':{key:1.}}}))


    def test_current_reference_snapshot_has_verified_scoped_members(self):
        import hashlib, json
        root = REPO / 'reference_tcad/bvmethods_sentaurus2018'
        pointer = json.loads((root / 'current_validation.json').read_text(encoding='utf-8'))
        self.assertFalse(pointer['native_store_registration_updated'])
        self.assertEqual(json.loads((root / 'artifact_manifest.json').read_text(encoding='utf-8'))['version'], 'v4')
        manifest_path = root / pointer['manifest']
        self.assertEqual(hashlib.sha256(manifest_path.read_bytes()).hexdigest(), pointer['manifest_sha256'])
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        self.assertEqual(manifest['kind'], 'local_text_qualification_not_native_store_registration')
        for name, receipt in manifest['files'].items():
            with self.subTest(member=name):
                path = manifest_path.parent / name
                self.assertTrue(path.resolve().is_relative_to(manifest_path.parent.resolve()))
                self.assertNotIn(path.suffix, ('.tdr', '.plt', '.h5', '.log', '.bin'))
                self.assertEqual(len(path.read_bytes()), receipt['bytes'])
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), receipt['sha256'])
        self.assertEqual(manifest['files']['source/pp5_des.cmd']['sha256'],
                         '47ed4a2330f20bcb190535175aa5dc6fb180cf49315093659c79d8e450c2a518')

    def test_current_reference_configs_preserve_the_two_stage_contract(self):
        import json
        root = REPO / 'reference_tcad/bvmethods_sentaurus2018/u2022_12_20261109'
        sidecar = json.loads((root / 'source/pp5_des.vela.json').read_text(encoding='utf-8'))
        sv.validate(sv.load_schema(REPO / 'schemas/vela.sdevice_sidecar.v1.schema.json'), sidecar)
        schema = sv.load_schema(REPO / 'configs/schema/vela-simulation.schema.json')
        decks = [json.loads((root / ('vela/simulation_stage%02d.json' % i)).read_text(encoding='utf-8')) for i in (0, 1)]
        for deck in decks:
            sv.validate(schema, deck)
            self.assertEqual(deck['physical_constants_profile'], 'sentaurus_u2022_12')
            self.assertEqual(deck['solver']['primary_variable_coordinates'], 'unit')
            self.assertEqual(deck['solver']['state_variables'], 'potential_density')
            self.assertEqual(deck['solver']['population_derivative_recovery'], 'after_protection')
            self.assertEqual(deck['solver']['population_stopping_rule'], 'native_approx')
            self.assertEqual(deck['solver']['band_to_band']['model'], 'none')
        self.assertEqual(decks[0]['sweep']['initialization']['poisson_max_iter'], 100)
        self.assertEqual(decks[1]['sweep']['initial_state_file'], decks[0]['sweep']['write_state_file'])
        circuit = decks[1]['sweep']['external_circuit']
        self.assertEqual(circuit['max_iterations'], 20)
        self.assertEqual(circuit['resistance_ohm_um'], 1e7)
        self.assertEqual(decks[1]['sweep']['stop'], 100000.)
        materials = json.loads((root / 'vela/materials.json').read_text(encoding='utf-8'))
        self.assertEqual(materials['physical_constants_profile'], 'sentaurus_u2022_12')
        self.assertEqual({x['name'] for x in materials['materials']}, {'Si', 'SiO2', 'Nitride'})

    def test_current_reference_curves_keep_native_and_candidate_roles(self):
        import csv, json, math
        root = REPO / 'reference_tcad/bvmethods_sentaurus2018/u2022_12_20261109'
        summary = json.loads((root / 'comparison.json').read_text(encoding='utf-8'))
        for name, case in summary['candidates'].items():
            self.assertEqual(case['role'], 'Vela_candidate_not_native_reference')
            self.assertIsNone(case['attempts']); self.assertIsNone(case['failed'])
            with (root / case['curve']).open(encoding='utf-8', newline='') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), case['accepted'])
            self.assertEqual(sum(int(x['newton_iterations']) for x in rows), case['successful_newton_iterations'])
            volts = [float(x['outer_voltage_V']) for x in rows]
            self.assertEqual(volts, sorted(set(volts)))
            self.assertTrue(all(math.isfinite(float(v)) for x in rows for v in x.values()))
            self.assertGreaterEqual(abs(float(rows[-1]['drain_total_current_A_per_um'])), .0013)
            self.assertLessEqual(abs(case['BV_inner_V'] / summary['contract']['existing_BV_acceptance_anchor_V'] - 1), .005)
            self.assertEqual(case['F_pass'], len(rows)); self.assertEqual(case['E_pass'], len(rows))
        native = summary['native']
        self.assertEqual(native['role'], 'Sentaurus_native_reference')
        self.assertEqual(native['version'], 'U-2022.12')
        self.assertEqual(native['raw_accepted_csv_rows'], 47)
        self.assertEqual(native['zero_initialization_records'], 3)
        with (root / native['curve']).open(encoding='utf-8', newline='') as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 45)
        self.assertEqual(float(rows[0]['outer_voltage_V']), 0.)
        self.assertEqual(float(rows[-1]['outer_voltage_V']), 100.)
        self.assertFalse(native['full_BV_qualified'])
        self.assertFalse(summary['qualifications']['strict_native_stopping_equivalence'])
        self.assertFalse(summary['qualifications']['speed_noninferiority'])
        self.assertEqual(summary['candidates']['cmd_full_bv']['carrier_report_failed'], 70)
        self.assertEqual(summary['candidates']['cmd_full_bv']['continuity_report_failed'], 13)

    def test_reference_freezer_rejects_overwrite_and_unqualified_evidence(self):
        import tempfile
        from scripts.freeze_bvmethods_u2022_validation import freeze, VERSION
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); existing = root / VERSION; existing.mkdir()
            sentinel = existing / 'retained'; sentinel.write_bytes(b'original')
            with self.assertRaisesRegex(ValueError, 'already exists'):
                freeze(root / 'absent', root / 'absent', root / 'absent', existing)
            self.assertEqual(sentinel.read_bytes(), b'original')
            bad = root / 'bad.zip'; bad.write_bytes(b'unqualified')
            fresh = root / 'fresh' / VERSION
            with self.assertRaisesRegex(ValueError, 'unqualified'):
                freeze(bad, root / 'absent', root / 'absent', fresh)
            self.assertFalse(fresh.exists())


if __name__=='__main__':unittest.main()
