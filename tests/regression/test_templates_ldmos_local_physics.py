"""Physical and data-contract checks for native local-state diagnostics."""
import json,math,os,subprocess,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from audit_templates_ldmos_thermal_physics import auger_parameters,auger_at_native_density

class AugerIsolationTest(unittest.TestCase):
    def test_units_balance_and_generation_sign(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'model.par'
            path.write_text('Auger * coefficients:\n{\n A = 1e-31, 2e-31\n B = 0, 0\n C = 0, 0\n H = 0, 0\n N0 = 1e18, 2e18\n}\n')
            p=auger_parameters(path)
            self.assertEqual(p['N0'],[1e24,2e24])
            self.assertAlmostEqual(p['A'][0]/1e-43,1.)
            self.assertEqual(auger_at_native_density(2e23,3e22,450.,0.,0.,p),0.)
            for split in (-.01,.01):
                value=auger_at_native_density(2e23,3e22,450.,0.,split,p,with_generation=True)
                expected=(1e-43*2e23+2e-43*3e22)*2e23*3e22*(1-math.exp(-split/(1.380649e-23/1.602176634e-19*450.)))
                self.assertAlmostEqual(value/expected,1.,places=12)
                if split<0:self.assertEqual(auger_at_native_density(2e23,3e22,450.,0.,split,p),0.)
            with self.assertRaises(ValueError):auger_at_native_density(0.,1.,300.,0.,0.,p)

class NativeBundleTest(unittest.TestCase):
    def test_bundle_preserves_physics_and_uses_load_plot_only_with_lf(self):
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory);native=base/'native';control=base/'control';out=base/'out'
            (native/'bundle').mkdir(parents=True);(native/'raw').mkdir();control.mkdir()
            prefix='Electrode { hRecVelocity=1.93e6 }\nPhysics { Fermi }\nPlot { SRH Band2Band }\nFile { Output="vg8_closed.log" }\n'
            (control/'vg8_closed.cmd').write_text(prefix+'Solve { Coupled { Poisson Electron Hole Temperature } }\n')
            for name in ('n1_fps.tdr','sdevice.par','Siliconc100.par'):(native/'bundle'/name).write_bytes(b'audited input')
            for gate in (4,8):
                for index in (0,1,10,30):(native/'raw'/f'field_vg{gate}_{index:04d}_des.tdr').write_bytes(b'saved state')
            subprocess.run([sys.executable,str(ROOT/'scripts/prepare_templates_ldmos_local_physics.py'),'--native-full',str(native),'--control-bundle',str(control),'--output',str(out)],check=True)
            self.assertNotIn(b'\r',(out/'run.sh').read_bytes())
            self.assertEqual(len(json.loads((out/'manifest.json').read_text())['cases']),8)
            for path in out.glob('*.cmd'):
                text=path.read_text();self.assertNotIn('Coupled',text)
                self.assertIn('hRecVelocity=1.93e6',text);self.assertIn('SRH Auger Band2Band',text)
                self.assertIn('Load(FilePrefix=',text);self.assertNotIn(b'\r',path.read_bytes())
            single=base/'single'
            subprocess.run([sys.executable,str(ROOT/'scripts/prepare_templates_ldmos_local_physics.py'),'--native-full',str(native),'--control-bundle',str(control),'--output',str(single),'--single-process'],check=True)
            batch=(single/'local_batch.cmd').read_text()
            self.assertNotIn('Coupled',batch);self.assertEqual(batch.count('Load(FilePrefix='),8)
            self.assertEqual(batch.count('Plot(FilePrefix='),8)
            self.assertEqual((single/'run.sh').read_text().count('/bin/sdevice '),1)

class RepresentativeSummaryTest(unittest.TestCase):
    def test_reference_alignment_uses_voltage_after_duplicate_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);native=root/'native';(native/'normalized').mkdir(parents=True)
            points=[];cases=[]
            for gate,tag in ((4,1),(8,2)):
                lines=['bias_V,current_total_A_per_um']
                if gate==8:lines.append('0,1e-30')
                lines += [f'{i*40/30},{gate*i*40/30*1e-6}' for i in range(31)]
                (native/'normalized'/f'IdVd_Vg{tag}_n4_des_drain_curve.csv').write_text('\n'.join(lines))
                for index in (0,1,10,30):
                    bias=index*40/30
                    points.append(dict(gate=gate,index=index,bias_V=bias,current_A_per_m=gate*bias,updates=0,wall_seconds=1.,thermal=dict(gates=dict(peak_temperature_rise=dict(error_K=0.),temperature_rise_field=dict(rms_error_K=0.)))))
                    cases.append(dict(name=f'local_vg{gate}_{index:02d}',gate_V=gate,index=index,bias_V=bias,fixed_native=dict(auger_m3_per_s=dict(max_abs=0.)),reclosed_A={},native_state_max_difference={}))
            (root/'points.json').write_text(json.dumps(dict(status='pass',points=points)))
            (root/'local.json').write_text(json.dumps(dict(status='completed',cases=cases)))
            subprocess.run([sys.executable,str(ROOT/'scripts/analyze_templates_ldmos_local_physics.py'),'--points',str(root/'points.json'),'--local',str(root/'local.json'),'--legacy-local',str(root/'local.json'),'--native',str(native),'--output',str(root/'summary.json')],check=True,capture_output=True,text=True)
            result=json.loads((root/'summary.json').read_text())
            self.assertFalse(result['full_curve_qualification'])
            for gate in ('4','8'):self.assertLess(result['per_gate'][gate]['max_selected_current_error_percent'],1e-12)

@unittest.skipUnless(os.environ.get('VELA_SILICON_THERMAL_PROBE'),'Matching build probe required')
class ReferencedProbeTest(unittest.TestCase):
    def test_saved_sub_ulp_splitting_preserves_recombination(self):
        with tempfile.TemporaryDirectory() as directory:
            source=Path(directory)/'input.json'
            state=dict(id=0,potential_V=40.65,electron_qf_V=0.,hole_qf_V=0.,temperature_K=515.,donors_m3=1e23,acceptors_m3=0.,electron_qf_reference_V=40.,hole_qf_reference_V=40.)
            states=[state,dict(state,id=1,hole_qf_V=1e-18)]
            source.write_text(json.dumps(dict(states_SI=states)))
            result=subprocess.run([os.environ['VELA_SILICON_THERMAL_PROBE'],str(source)],check=True,capture_output=True,text=True)
            rows=json.loads(result.stdout)['results']
            self.assertEqual(rows[0]['auger_m3_per_s']['value'],0.)
            self.assertGreater(rows[1]['auger_m3_per_s']['value'],0.)
            self.assertGreater(rows[1]['srh_m3_per_s']['value'],0.)
            state.update(reference_conduction_band_eV=rows[0]['conduction_band_eV']['value'],reference_valence_band_eV=rows[0]['valence_band_eV']['value'])
            source.write_text(json.dumps(dict(states_SI=[state])))
            aligned=json.loads(subprocess.run([os.environ['VELA_SILICON_THERMAL_PROBE'],str(source)],check=True,capture_output=True,text=True).stdout)['results'][0]
            for carrier in ('electrons','holes'):
                self.assertAlmostEqual(aligned[f'{carrier}_at_native_band_m3']/aligned[f'{carrier}_m3']['value'],1.,places=11)
            state['hole_qf_V']=-1e-18
            for enabled in (False,True):
                source.write_text(json.dumps(dict(states_SI=[state],auger_with_generation=enabled)))
                diagnostic=json.loads(subprocess.run([os.environ['VELA_SILICON_THERMAL_PROBE'],str(source)],check=True,capture_output=True,text=True).stdout)
                self.assertEqual(diagnostic['auger_with_generation'],enabled)
                rate=diagnostic['results'][0]['auger_m3_per_s']['value']
                if enabled:self.assertLess(rate,0.)
                else:self.assertEqual(rate,0.)

if __name__=='__main__':unittest.main()
