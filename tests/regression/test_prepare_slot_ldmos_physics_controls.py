"""Isolation, native parameter selection and unit safeguards for SLOT controls."""
import csv
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from prepare_slot_ldmos_physics_controls import prepare,prepare_poly,native_values

PAR='''Material = "Oxide" {
Chi0=0.9
Eg0=9
}
Material = "Silicon" {
Chi0=4.0727
Eg0=1.12416
Auger * parameters
{
A=6.7e-32,7.2e-32
B=2.45e-31,4.5e-33
C=-2.2e-32,2.63e-32
H=3.46667,8.25688
N0=1e18,1e18
}
}
'''

class PhysicsControls(unittest.TestCase):
    def setup_files(self,folder):
        native=folder/'native.par';native.write_text(PAR)
        state=folder/'seed.csv';state.write_text('frozen seed')
        mats=folder/'materials.json';mats.write_text(json.dumps({'materials':[{'name':'Si','ni':1.0749877484352535e10},{'name':'PolySilicon','eps_r':11.7}]}))
        cfg={'scaling':{'mode':'unit_scaling'},'solver':{'auger_with_generation':False,
             'auger_cn_m6_per_s':2.9e-31,'auger_cp_m6_per_s':1.028e-31,
             'carrier_row_convergence':{'eps_row':.001}},
             'contacts':[{'name':'gate','type':'metal_gate','bias':0}],
             'materials_file':str(mats),'sweep':{'bias_points':[2.],'initial_state_file':str(state)}}
        base=folder/'base.json';base.write_text(json.dumps(cfg))
        return base,state,native,cfg

    def test_independent_controls_preserve_ni_contacts_and_acceptance(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);base,state,native,original=self.setup_files(folder)
            frozen={p:p.read_bytes() for p in (base,state,native,folder/'materials.json')}
            result=prepare(base,state,native,folder/'controls')
            for name,path in result['cases'].items():
                cfg=json.loads(Path(path).read_text())
                self.assertEqual(cfg['contacts'],original['contacts'])
                self.assertEqual(cfg['solver']['carrier_row_convergence'],original['solver']['carrier_row_convergence'])
                self.assertFalse(cfg['solver']['auger_with_generation'])
                self.assertEqual('auger_density_dependence' in cfg['solver'],name=='auger_density')
                mats=json.loads(Path(cfg['materials_file']).read_text())['materials']
                self.assertEqual(next(r['ni'] for r in mats if r['name']=='Si'),1.0749877484352535e10)
                self.assertEqual(next(r for r in mats if r['name']=='PolySilicon'),{'name':'PolySilicon','eps_r':11.7})
            density=json.loads(Path(result['cases']['auger_density']).read_text())['solver']['auger_density_dependence']
            self.assertEqual(density['electron'],{'enhancement':3.46667,'reference_density_m3':1e18})
            self.assertEqual(density['hole']['enhancement'],8.25688)
            for p,content in frozen.items():self.assertEqual(p.read_bytes(),content)
            with self.assertRaises(ValueError):prepare(base,state,native,folder/'controls')

    def test_rejects_ambiguous_units_and_temperature(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);base,state,native,cfg=self.setup_files(folder)
            cfg['scaling']['mode']='legacy_si';base.write_text(json.dumps(cfg))
            with self.assertRaises(ValueError):prepare(base,state,native,folder/'wrong_units')
            cfg['scaling']['mode']='unit_scaling';cfg['solver']['temperature_K']=400;base.write_text(json.dumps(cfg))
            with self.assertRaises(ValueError):prepare(base,state,native,folder/'wrong_temperature')

    def test_reads_material_scope_and_rejects_incomplete_auger(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'native.par';p.write_text(PAR)
            result=native_values(p)
            self.assertEqual(result['Si']['electron_affinity_eV'],4.0727)
            self.assertEqual(result['SiO2']['electron_affinity_eV'],.9)
            p.write_text(PAR.replace('N0=1e18,1e18',''))
            with self.assertRaises(ValueError):native_values(p)

    def test_poly_restoration_preserves_silicon_and_requires_native_units(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);base,state,native,cfg=self.setup_files(folder)
            def write_csv(path,rows):
                with path.open('w',newline='') as stream:
                    writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
            mesh={'nodes':[{'id':i,'x':float(i),'y':0.} for i in range(6)],
                  'triangles':[{'id':0,'region_id':0,'node_ids':[0,1,2]},
                               {'id':1,'region_id':3,'node_ids':[3,4,5]}],
                  'regions':[{'id':0,'material':'Si'},{'id':3,'material':'PolySilicon'}]}
            mesh_path=folder/'mesh.json';mesh_path.write_text(json.dumps(mesh))
            doping=[{'node_id':i,'donors_cm3':1e17 if i<3 else 0.,'acceptors_cm3':0.} for i in range(6)]
            doping_path=folder/'doping.csv';write_csv(doping_path,doping)
            states=[{'node_id':i,'psi':.1,'phin':0.,'phip':0.,'electrons_m3':1e20,'holes_m3':1e12} for i in range(6)]
            write_csv(state,states)
            cfg.update(mesh_file=str(mesh_path),node_doping_file=str(doping_path))
            cfg['solver']['mobility']={'model':'ialmob','ialmob':{}}
            cfg['contacts'] += [{'name':'SLOT','type':'metal_gate','flatband_voltage':0,'bias':0},
                                {'name':'drain','type':'ohmic','bias':2}]
            base.write_text(json.dumps(cfg))
            export=folder/'native';fields=export/'fields';fields.mkdir(parents=True)
            write_csv(export/'nodes.csv',[{'id':r['id'],'x_um':r['x'],'y_um':r['y']} for r in mesh['nodes']])
            write_csv(export/'elements.csv',[{'id':r['id'],**{'node'+str(k):v for k,v in enumerate(r['node_ids'])}} for r in mesh['triangles']])
            vt=.025851995266484913;ni=1.4e10;metadata=[]
            for name,unit,values in [
                ('eDensity','cm^-3',[ni*math.exp((.01*i+.2)/(2*vt)) for i in range(3)]),
                ('hDensity','cm^-3',[ni*math.exp((.01*i-.2)/(2*vt)) for i in range(3)]),
                ('eQuasiFermiPotential','V',[0.]*3),('hQuasiFermiPotential','V',[0.]*3),
                ('DonorConcentration','cm^-3',[1e18,2e18,3e18]),('AcceptorConcentration','cm^-3',[1e15]*3),
                ('BandgapNarrowing','eV',[.0,.01,.02]),('BandGap','eV',[1.1082]*3),
                ('ElectronAffinity','eV',[4.0727+.005*i for i in range(3)]),
                ('eMobility','cm^2*V^-1*s^-1',[1417.]*3),('hMobility','cm^2*V^-1*s^-1',[470.5]*3)]:
                write_csv(fields/(name+'_region3.csv'),[{'node_id':i+3,'component0':v} for i,v in enumerate(values)])
                metadata.append({'name':name,'region':3,'unit':unit})
            meta=export/'metadata.json';meta.write_text(json.dumps({'fields':metadata}))
            frozen={p:p.read_bytes() for p in (base,state,doping_path,folder/'materials.json')}
            result=prepare_poly(base,state,fields,folder/'poly')
            self.assertAlmostEqual(result['ni_cm3']/ni,1.,places=12)
            self.assertAlmostEqual(result['native_thermal_voltage_V'],vt,places=12)
            for name,path in result['cases'].items():
                control=json.loads(Path(path).read_text())
                self.assertEqual(control['solver']['recombination_materials'],['Si'])
                self.assertTrue(all(c['type']=='ohmic' for c in control['contacts']))
                self.assertEqual(control['solver']['carrier_row_convergence'],cfg['solver']['carrier_row_convergence'])
                with Path(control['node_doping_file']).open(newline='') as stream:
                    rows=list(csv.DictReader(stream))
                self.assertEqual(float(rows[0]['donors_cm3']),1e17)
                self.assertEqual(float(rows[3]['donors_cm3']),1e18 if name=='poly_doped' else 0.)
                with Path(control['sweep']['initial_state_file']).open(newline='') as stream:
                    seed=list(csv.DictReader(stream))
                self.assertEqual(float(seed[0]['psi']),.1)
                self.assertEqual(float(seed[0]['electrons_m3']),1e20)
            for p,content in frozen.items():self.assertEqual(p.read_bytes(),content)
            self.assertIn(str(fields/'eDensity_region3.csv'),result['sources'])
            metadata[0]['unit']='m^-3';meta.write_text(json.dumps({'fields':metadata}))
            with self.assertRaisesRegex(ValueError,'unit'):
                prepare_poly(base,state,fields,folder/'bad_units')

if __name__=='__main__':unittest.main()
