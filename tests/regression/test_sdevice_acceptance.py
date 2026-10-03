"""WP10 reference admission, sign, coverage and same-node field gates."""
from pathlib import Path
import base64
import copy
import json
import math
import struct
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'scripts'))
import sdevice_acceptance as a


class AdmissionTests(unittest.TestCase):
    def test_all_six_checks_required_and_partial_is_not_pass(self):
        report = {'status':'pass','mode':'online','compared_with_source':True,
                  'checks':[{'number':i,'status':'pass','items':[{'status':'pass'}]} for i in range(1,7)]}
        self.assertTrue(a.qualified(report))
        for change in ({'mode':'offline'}, {'partial_references':['early stop']}, {'compared_with_source':False},
                       {'checks':report['checks'][:5]}):
            self.assertFalse(a.qualified(dict(report,**change)))
        for status in ('skip','warn','fail','partial'):
            broken=copy.deepcopy(report);broken['checks'][2]['items'][0]['status']=status
            self.assertFalse(a.qualified(broken))

    def test_only_parameter_declaration_may_be_added(self):
        original='File { Grid="x.tdr" } Physics { Recombination(SRH) }'
        adapted=original.replace('File {','File { Parameter="models.par"')
        self.assertIn('inserted',a.augmented_source(original,adapted))
        for altered in (adapted.replace('SRH','Auger'),adapted+' Solve {Coupled{Poisson}}',adapted.replace('models.par','other.par')):
            with self.assertRaises(ValueError):a.augmented_source(original,altered)

    def test_missing_reference_environment_is_explicit_skip(self):
        import os
        from unittest.mock import patch
        with patch.dict(os.environ,{},clear=True):
            self.assertEqual(a.main(['--case','pn2d_0v','--output-dir','unused','--launcher','unused','--tdr-importer','unused']),77)


class CurveTests(unittest.TestCase):
    def setUp(self):
        self.contract=a.read_json(a.CONTRACT)
        self.case=self.contract['cases']['pn2d_iv']
        self.ref=[{'bias_V':b,'I':i} for b,i in [(0.,0.),(.1,1e-15),(.2,1e-13),(.3,1e-11)]]
        self.cand=[dict(r,currents=[r['I'],-r['I']]) for r in self.ref]

    def test_same_contact_three_point_sign_and_inconsistent_sign(self):
        self.assertEqual(a.sign_evidence(self.ref,self.cand)['factor'],1.)
        inverse=[dict(r,I=-r['I']) for r in self.cand]
        self.assertEqual(a.sign_evidence(self.ref,inverse)['factor'],-1.)
        inverse[-1]['I']*=-1
        with self.assertRaisesRegex(ValueError,'inconsistent'):a.sign_evidence(self.ref,inverse)
        with self.assertRaises(ValueError):a.sign_evidence(self.ref,self.cand[:2])

    def test_frozen_oom_sign_and_conservation_are_independent(self):
        self.assertEqual(a.curve_metrics(self.ref,self.cand,self.case,self.contract)['status'],'pass')
        for change in (lambda r:dict(r,I=10*r['I']),lambda r:dict(r,I=-r['I']),lambda r:dict(r,currents=[r['I'],0.])):
            self.assertEqual(a.curve_metrics(self.ref,[change(r) for r in self.cand],self.case,self.contract)['status'],'fail')
        with self.assertRaisesRegex(ValueError,'cover'):a.curve_metrics(self.ref,self.cand[:-1],self.case,self.contract)

    def test_nonfinite_measurements_fail_closed(self):
        for v in ('nan','inf','-inf'):
            with self.assertRaises(ValueError):a.finite(v)

    def test_opposite_sign_fails_even_without_shared_biases(self):
        cand=[{'bias_V':b,'I':-10**(20*b-17),'currents':[1e-12,-1e-12]} for b in (0.,.19,.31)]
        self.assertFalse(a.curve_metrics(self.ref,cand,self.case,self.contract)['same_sign'])


class FieldTests(unittest.TestCase):
    def test_inline_binary_and_malformed_shape(self):
        raw=struct.pack('<I2d',16,2.,3.)
        element=ET.fromstring('<DataArray type="Float64" format="binary">'+base64.b64encode(raw).decode()+'</DataArray>')
        self.assertEqual(a.array_values(element),[[2.],[3.]])
        element.set('NumberOfComponents','3')
        with self.assertRaises(ValueError):a.array_values(element)
        element=ET.fromstring('<DataArray>nan</DataArray>')
        with self.assertRaises(ValueError):a.array_values(element)

    def test_native_multistage_time_resets_and_accumulation_is_rejected(self):
        from vela_sdevice.plt import build_plt_document
        with tempfile.TemporaryDirectory() as temp:
            base=Path(temp);(base/'plt').mkdir()
            entry={'outer_voltage_V':0.,'inner_voltage_V':0.,'e_current_A_per_um':1.,'h_current_A_per_um':1.,'total_current_A_per_um':2.}
            paths=[base/'plt'/f'bjt_m0_{stage}_bjt_m0_des.plt' for stage in ('base','collector')]
            rows=[dict(time=t,base=entry,collector=entry) for t in (0.,.5,1.)]
            for p in paths:p.write_text(build_plt_document(['base','collector'],rows,sign_factor=1.))
            self.assertEqual(a.time_semantics(base)['semantics'],'reset_per_quasistationary_stage')
            paths[1].write_text(build_plt_document(['base','collector'],[dict(r,time=r['time']+1.) for r in rows],sign_factor=1.))
            with self.assertRaises(ValueError):a.time_semantics(base)

    def test_same_node_field_gate_detects_identity_units_and_values(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);(root/'fields').mkdir()
            (root/'nodes.csv').write_text('id,x_um,y_um\n0,1,0\n1,1,1\n2,2,0\n')
            mesh={'nodes':[{'id':0,'x':1.,'y':0.},{'id':1,'x':1.,'y':1.},{'id':2,'x':2.,'y':0.}],
                  'triangles':[{'node_ids':[0,1,2]}]}
            names=('ElectrostaticPotential','eDensity','hDensity','ElectricField')
            units=('V','cm^-3','cm^-3','V/cm')
            fields=[]
            arrays=[];metadata=[]
            for name,unit in zip(names,units):
                vector=name=='ElectricField';components=2 if vector else 1
                value='3,4' if vector else '1000000'
                (root/'fields'/f'{name}.csv').write_text('node_id,component0'+(',component1' if vector else '')+'\n'+''.join(f'{i},{value}\n' for i in range(3)))
                fields.append({'name':name,'region':0,'mapping_status':'complete','support_kind':'node','components':components,'unit':unit,'csv_file':name+'.csv'})
                arrays.append(f'<DataArray Name="{name}" NumberOfComponents="{components}">'+('3 4 '*3 if vector else '1000000 '*3)+'</DataArray>')
                metadata.append(f'<DataArray type="String" Name="Unit:{name}">'+ ' '.join(str(ord(c)) for c in unit)+' 0</DataArray>')
            (root/'field_manifest.json').write_text(json.dumps({'fields':fields}))
            text='<VTKFile><UnstructuredGrid><FieldData><DataArray type="Float64" Name="Bias_V">0</DataArray><DataArray type="Float64" Name="Converged">1</DataArray>'+''.join(metadata)+'</FieldData><Piece><Points><DataArray NumberOfComponents="3">1 0 0 1 1 0 2 0 0</DataArray></Points><PointData>'+''.join(arrays)+'</PointData></Piece></UnstructuredGrid></VTKFile>'
            vtu=root/'candidate.vtu';vtu.write_text(text)
            limits=a.read_json(a.CONTRACT)['spatial']
            self.assertEqual(a.field_metrics(root,vtu,mesh,limits)['status'],'pass')
            # Native TDR exports both magnitude and vector under the same name.
            fields.append(dict(fields[-1],components=1));(root/'field_manifest.json').write_text(json.dumps({'fields':fields}))
            self.assertEqual(a.field_metrics(root,vtu,mesh,limits)['status'],'pass')
            vtu.write_text(text.replace('1 0 0 1 1 0','1.01 0 0 1 1 0'))
            with self.assertRaisesRegex(ValueError,'coordinates'):a.field_metrics(root,vtu,mesh,limits)
            vtu.write_text(text.replace('1000000 1000000 1000000','100000000 100000000 100000000'))
            self.assertEqual(a.field_metrics(root,vtu,mesh,limits)['status'],'fail')
            vtu.write_text(text);fields[0]['unit']='mV';(root/'field_manifest.json').write_text(json.dumps({'fields':fields}))
            with self.assertRaisesRegex(ValueError,'units'):a.field_metrics(root,vtu,mesh,limits)


if __name__=='__main__':unittest.main()
