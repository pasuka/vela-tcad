"""Validate explicit LDMOS electrothermal input preparation and provenance."""
import csv,json,os,subprocess,sys,tempfile,unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'scripts'))
SCRIPT=ROOT/'scripts/prepare_templates_ldmos_electrothermal.py'
class PreparationTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.base=Path(self.temp.name)
        self.export=self.base/'export';(self.export/'fields').mkdir(parents=True)
        nodes=[dict(id=i,x=x,y=y) for i,(x,y) in enumerate(((0,0),(2,0),(.5,.5),(0,1)))]
        mesh=dict(nodes=nodes,triangles=[dict(id=0,region_id=0,node_ids=[0,1,2]),dict(id=1,region_id=1,node_ids=[0,2,3])],contacts=[dict(name='source',node_ids=[0]),dict(name='drain',node_ids=[1]),dict(name='gate',node_ids=[3]),dict(name='th_lat',node_ids=[2])])
        self.save('mesh.json',mesh);self.save('geometry.json',{})
        self.save('thermal.json',dict(mesh_file=str(self.base/'mesh.json'),coordinate_to_metres=1e-6,region_conductivity=[],thermodes=[]))
        self.save('profile.json',dict(contacts=[dict(name='gate',flatband_voltage=-.56)],solver=dict(mobility=dict(electron_saturation_velocity_m_s=1.07e7,hole_saturation_velocity_m_s=8.37e6,ialmob=dict(geometry_file=str(self.base/'geometry.json'))),carrier_row_convergence=dict(eps_row=1e-8),continuity_row_scaling={},block_absolute_convergence={})))
        (self.base/'couples.csv').write_text('node0,node1,couple_m\n')
        for name,value in [('ElectrostaticPotential',.5),('eQuasiFermiPotential',0.),('hQuasiFermiPotential',0.),('LatticeTemperature',300.),('DonorConcentration',1e17),('AcceptorConcentration',0.)]:
            with (self.export/'fields'/f'{name}_region0.csv').open('w',newline='') as f:
                w=csv.writer(f);w.writerow(['node_id','component0']);w.writerows((i,value) for i in range(3 if name.endswith('Concentration') else 4))
    def save(self,name,data):
        (self.base/name).write_text(json.dumps(data),encoding='utf-8')
    def run_prepare(self,*extra):
        return subprocess.run([sys.executable,'-X','utf8',str(SCRIPT),'--thermal-input',str(self.base/'thermal.json'),'--export',str(self.export),'--profile',str(self.base/'profile.json'),'--couples',str(self.base/'couples.csv'),'--output',str(self.base/'out'),'--gate','8','--drain','40','--freeze-temperature','300',*extra],cwd=ROOT,capture_output=True,text=True)
    def test_units_gate_sign_and_qualified_obtuse_geometry(self):
        run=self.run_prepare();self.assertEqual(run.returncode,0,run.stderr)
        output=json.loads((self.base/'out/input.json').read_text());manifest=json.loads((self.base/'out/manifest.json').read_text())
        self.assertAlmostEqual(sum(output['silicon_area_m2']),.5e-12,delta=1e-26)
        edge=next(e for e in output['edge_geometry'] if e['nodes']==[0,1])
        # Qualified negative-cotangent fallback is the barycentric dual
        # length divided by edge length, here area/(3*length**2).
        self.assertAlmostEqual(edge['poisson_F_per_m'],11.7*8.8541878128e-12*.5/(3*4),delta=1e-25)
        gate=next(b for b in output['boundaries'] if b['kind']=='psi');self.assertEqual(gate['value'],8.56)
        self.assertEqual(len([b for b in output['boundaries'] if b['kind']=='temperature']),4)
        self.assertEqual(output['mobility_SI']['electron_saturation_velocity_m_s'],1.07e5)
        self.assertAlmostEqual(output['recombination_area_m2'][0],.25e-12,delta=1e-26)
        self.assertEqual(output['recombination_area_m2'][3],0.)
        self.assertEqual(output['donors_m3'],[1e23,1e23,1e23,0.])
        self.assertIn(str(self.base/'geometry.json'),manifest['sources_sha256'])

    def test_restart_retains_sub_ulp_quasi_fermi_increments(self):
        state=self.base/'state.csv'
        with state.open('w',newline='') as f:
            names=['node_id','psi','phin','phip','electron_qf_increment_V','hole_qf_increment_V','electron_qf_reference_V','hole_qf_reference_V']
            w=csv.writer(f);w.writerow(names)
            for i in range(4):w.writerow([i,-3.5,-4.,-4.,1e-20,-2e-20,-4.,-4.])
        run=self.run_prepare('--isothermal-state',str(state));self.assertEqual(run.returncode,0,run.stderr)
        output=json.loads((self.base/'out/input.json').read_text())
        self.assertEqual(output['state_interleaved'][1],0.)
        self.assertEqual(output['electron_qf_reference_V'][0],0.)
        self.assertEqual(output['referenced_state_interleaved'][1],1e-20)
        self.assertEqual(output['referenced_state_interleaved'][2],-2e-20)
    def test_missing_semiconductor_doping_is_rejected(self):
        (self.export/'fields/DonorConcentration_region0.csv').write_text('node_id,component0\n0,1e17\n1,1e17\n')
        run=self.run_prepare();self.assertNotEqual(run.returncode,0)
        self.assertIn('Incomplete semiconductor field DonorConcentration',run.stderr)
        self.assertFalse((self.base/'out/input.json').exists())

class CurveGateTest(unittest.TestCase):
    def result(self):
        return dict(diagnostic_stop='diagnostic_scaled_residual',electrical_block_gates=[dict(satisfied=True)]*3,carrier_row_gate=dict(satisfied=True),
            contacts=[dict(contact=name,total_outflow_A_per_m=value) for name,value in [('drain',1.),('source',-1.),('gate',0.),('substrate',0.)]],
            temperature_K=[310.,312.],lattice_source_W_per_m=2.,boundary_heat_W_per_m=2.)
    def test_electrical_and_thermal_gates_are_both_required(self):
        from scripts.run_templates_ldmos_electrothermal_curve import state_gate
        r=self.result();self.assertTrue(state_gate(r,2.)['pass_gate'])
        r['carrier_row_gate']['satisfied']=False;self.assertFalse(state_gate(r,2.)['pass_gate'])
        r=self.result();r['boundary_heat_W_per_m']=1.9;self.assertFalse(state_gate(r,2.)['pass_gate'])
        r=self.result();r['contacts'][1]['total_outflow_A_per_m']=-.99;self.assertFalse(state_gate(r,2.)['pass_gate'])
        r=self.result();r['diagnostic_stop']='diagnostic_iteration_limit';self.assertFalse(state_gate(r,2.)['pass_gate'])
    def test_zero_power_requires_exact_constant_temperature_evidence(self):
        from scripts.run_templates_ldmos_electrothermal_curve import state_gate
        r=self.result()
        for c in r['contacts']:c['total_outflow_A_per_m']=0.
        r.update(temperature_K=[300.,300.],lattice_source_W_per_m=0.,boundary_heat_W_per_m=0.)
        self.assertTrue(state_gate(r,0.)['pass_gate']);self.assertIsNone(state_gate(r,0.)['heat_balance_relative'])
        r['temperature_K'][0]=301.;self.assertFalse(state_gate(r,0.)['pass_gate'])
@unittest.skipUnless(os.environ.get('VELA_ELECTROTHERMAL_PROBE'),'Probe supplied by matching CTest build')
class ProbeInputTest(unittest.TestCase):
    def test_truncated_state_is_rejected_before_reference_access(self):
        with tempfile.TemporaryDirectory() as name:
            base=Path(name);mesh=base/'mesh.json';source=base/'input.json';target=base/'output.json'
            mesh.write_text(json.dumps(dict(nodes=[dict(id=i,x=x,y=y) for i,(x,y) in enumerate(((0,0),(1,0),(0,1)))],triangles=[dict(id=0,region_id=0,node_ids=[0,1,2])],regions=[dict(id=0,name='silicon',material='Silicon',cell_ids=[0])],contacts=[])),encoding='utf-8')
            cfg=dict(mesh_file=str(mesh),coordinate_to_metres=1.,region_conductivity=[dict(region_id=0,model='constant',value_W_per_m_K=1.)],thermodes=[],silicon_area_m2=[1/6]*3,fixed_charge_C_per_m=[0.]*3,edge_geometry=[dict(nodes=[a,b],poisson_F_per_m=1e-10,transport_weight=0.) for a,b in ((0,1),(1,2),(0,2))],donors_m3=[0.]*3,acceptors_m3=[0.]*3,mobility_SI=dict(model='constant'),boundaries=[],state_interleaved=[0.,0.,0.,300.]*3)
            source.write_text(json.dumps(cfg),encoding='utf-8')
            run=lambda:subprocess.run([os.environ['VELA_ELECTROTHERMAL_PROBE'],str(source),str(target)],capture_output=True,text=True)
            valid=run();self.assertEqual(valid.returncode,0,valid.stderr)
            target=base/'invalid_output.json';cfg['state_interleaved'].pop();source.write_text(json.dumps(cfg),encoding='utf-8')
            invalid=run();self.assertNotEqual(invalid.returncode,0);self.assertIn('Invalid electrothermal state size',invalid.stderr);self.assertFalse(target.exists())

class FullD0GateTest(unittest.TestCase):
    def setUp(self):
        import hashlib
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.base=Path(self.temp.name)
        self.native=self.base/'native';(self.native/'normalized').mkdir(parents=True);self.curves={};fields=[]
        self.contract=json.loads((ROOT/'reference_tcad/templates_ldmos_sentaurus2022/thermal/d0_thermal_acceptance.json').read_text(encoding='utf-8'))
        for gate,number in ((4,1),(8,2)):
            directory=self.base/f'vg{gate}';directory.mkdir();self.curves[gate]=directory;points=[];runs=[];curve=['bias_V,current_total_A_per_um'];balance=['point_index,bias_V,contact,current_total_A_per_um']
            for i in range(31):
                bias=i*40/30;current=gate*bias*1e-6;temperature=[300.+bias,300.+bias*.5]
                result=CurveGateTest().result();result.update(temperature_K=temperature,nodal_area_m2=[1.,2.],lattice_source_W_per_m=bias,boundary_heat_W_per_m=bias)
                for c in result['contacts']:c['total_outflow_A_per_m']=current*1e6*(1 if c['contact']=='drain' else -1 if c['contact']=='source' else 0)
                case=directory/f'point_{i}';case.mkdir();rp=case/'output.json';rp.write_text(json.dumps(result));points.append(dict(bias_V=bias,result=str(rp)));runs.append(dict(bias_V=bias,directory=str(case),gate=dict(pass_gate=True)))
                tp=self.native/f'{gate}_{i}.json';tp.write_text(json.dumps(dict(node_id=[0,1],temperature_K=temperature)))
                fields.append(dict(gate_V=gate,point_index=i,bias_V=bias,temperature_file=str(tp),temperature_sha256=hashlib.sha256(tp.read_bytes()).hexdigest()))
                curve.append(f'{bias},{current}')
                for c in result['contacts']:balance.append(f"{i},{bias},{c['contact']},{c['total_outflow_A_per_m']*1e-6}")
            (directory/'ledger.json').write_text(json.dumps(dict(status='complete',exact_points=points,runs=runs)))
            (directory/'curve.csv').write_text('\n'.join(curve));(directory/'terminal_balance.csv').write_text('\n'.join(balance))
            (self.native/'normalized'/f'IdVd_Vg{number}_n4_des_drain_curve.csv').write_text('\n'.join(curve))
        (self.native/'manifest.json').write_text(json.dumps(dict(fields=fields)))
    def score(self):
        from analyze_templates_ldmos_d0 import analyze
        return analyze(self.native,self.curves,self.contract)
    def test_interior_temperature_failure_cannot_hide_behind_endpoint(self):
        self.assertEqual(self.score()['status'],'pass')
        path=self.curves[4]/'point_10/output.json';data=json.loads(path.read_text());data['temperature_K'][0]+=20.;path.write_text(json.dumps(data))
        self.assertEqual(self.score()['status'],'fail')
    def test_incomplete_curve_is_rejected(self):
        path=self.curves[8]/'ledger.json';data=json.loads(path.read_text());data['exact_points'].pop();path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError,'31-point'):self.score()
if __name__=='__main__':unittest.main()
