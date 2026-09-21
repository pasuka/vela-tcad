"""Independent band identity and unit checks for SLOT reference preparation."""
import csv
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from prepare_slot_ldmos_energy_reference import derive_reference
from run_slot_ldmos_noavalanche_continuation import render_cold_equilibrium


class BandReference(unittest.TestCase):
    def test_recovers_common_origin_and_material_difference(self):
        with tempfile.TemporaryDirectory() as tmp:
            export=Path(tmp);fields=export/'fields';fields.mkdir();metadata=[]
            vt=.025851995266484913;origin=4.6337;expected={}
            for region,chi,eg,nc,nv in [(0,4.0727,1.12416,2.85665e19,3.10463e19),
                                       (3,4.07274038,1.10820923,2.85825e19,3.10463e19)]:
                ni=math.sqrt(nc*nv)*math.exp(-eg/(2*vt))
                offset=origin-chi-eg/2-vt/2*math.log(nc/nv)
                expected[region]=(ni,offset)
                values={k:[] for k in ['eDensity','hDensity','ElectrostaticPotential',
                    'eQuasiFermiPotential','hQuasiFermiPotential','ConductionBandEnergy',
                    'ValenceBandEnergy','ElectronAffinity','BandGap','BandgapNarrowing']}
                for i in range(3):
                    psi=.15+i*.1;bgn=i*.02;fn=.001*i;fp=.004*i
                    ec=origin-psi-chi-bgn/2;ev=ec-eg+bgn
                    row={'eDensity':nc*math.exp((-fn-ec)/vt),'hDensity':nv*math.exp((ev+fp)/vt),
                        'ElectrostaticPotential':psi,'eQuasiFermiPotential':fn,'hQuasiFermiPotential':fp,
                        'ConductionBandEnergy':ec,'ValenceBandEnergy':ev,'ElectronAffinity':chi+bgn/2,
                        'BandGap':eg,'BandgapNarrowing':bgn}
                    for k,v in row.items():values[k].append({'node_id':i,'component0':v})
                for name,rows in values.items():
                    unit='cm^-3' if name.endswith('Density') else ('V' if name.endswith('Potential') else 'eV')
                    metadata.append({'name':name,'region':region,'unit':unit})
                    with (fields/f'{name}_region{region}.csv').open('w',newline='') as f:
                        writer=csv.DictWriter(f,fieldnames=['node_id','component0']);writer.writeheader();writer.writerows(rows)
            meta=export/'metadata.json';meta.write_text(json.dumps({'fields':metadata}))
            result=derive_reference(export,vt)
            for material,region in [('Si',0),('PolySilicon',3)]:
                self.assertAlmostEqual(result[material]['vacuum_reference_eV'],origin,places=12)
                self.assertAlmostEqual(result[material]['ni_cm3']/expected[region][0],1.,places=12)
                self.assertAlmostEqual(result[material]['intrinsic_reference_V'],expected[region][1],places=12)
            self.assertAlmostEqual(result['PolySilicon']['relative_to_Si_V'],expected[3][1]-expected[0][1],places=12)
            metadata[0]['unit']='m^-3';meta.write_text(json.dumps({'fields':metadata}))
            with self.assertRaisesRegex(ValueError,'unit'):derive_reference(export,vt)

    def test_rejects_invalid_thermal_voltage(self):
        for vt in [0.,-1.,float('nan'),float('inf')]:
            with self.assertRaises(ValueError):derive_reference(Path('not_read'),vt)

    def test_cold_equilibrium_has_no_imported_solution_and_preserves_physics_and_gates(self):
        base={'state_file':'foreign.csv','state_fields_dir':'foreign',
            'contacts':[{'name':'drain','bias':0.}],
            'solver':{'abstol':1e-9,'recombination_volume_policy':'material_local_mixed_voronoi',
                      'carrier_row_convergence':{'eps_row':.001}},
            'sweep':{'contact':'drain','initial_state_file':'foreign.csv'}}
        cold,_,_=render_cold_equilibrium(base,Path('output'))
        self.assertNotIn('initial_state_file',cold['sweep'])
        self.assertNotIn('state_file',cold);self.assertNotIn('state_fields_dir',cold)
        self.assertEqual(cold['sweep']['bias_points'],[0.])
        self.assertEqual(cold['solver']['abstol'],base['solver']['abstol'])
        self.assertEqual(cold['solver']['carrier_row_convergence'],base['solver']['carrier_row_convergence'])
        self.assertEqual(cold['solver']['recombination_volume_policy'],base['solver']['recombination_volume_policy'])
        self.assertEqual(base['sweep']['initial_state_file'],'foreign.csv')
        base['contacts'][0]['bias']=1.
        with self.assertRaisesRegex(ValueError,'zero terminal'):render_cold_equilibrium(base,Path('output'))


if __name__=='__main__':unittest.main()
