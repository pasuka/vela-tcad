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
from unittest.mock import patch
from types import SimpleNamespace
import shutil

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(REPO/'scripts'))
import sdevice_acceptance as a


def field_fixture(root, bias=0.0):
    """Three-node native export and identical accepted vector VTU, in physical units."""
    root.mkdir(parents=True, exist_ok=True)
    (root/'fields').mkdir()
    (root/'nodes.csv').write_text('id,x_um,y_um\n0,1,0\n1,1,1\n2,2,0\n')
    mesh = {'nodes':[{'id':0,'x':1.,'y':0.},{'id':1,'x':1.,'y':1.},{'id':2,'x':2.,'y':0.}],
            'triangles':[{'node_ids':[0,1,2]}]}
    fields, arrays, metadata = [], [], []
    for name, unit, values in [('ElectrostaticPotential','V',[0.]), ('eDensity','cm^-3',[1e6]),
                              ('hDensity','cm^-3',[1e6]), ('ElectricField','V/cm',[3.,4.])]:
        (root/'fields'/f'{name}.csv').write_text('node_id,'+','.join(f'component{i}' for i in range(len(values)))+'\n'
            + ''.join(f'{i},'+','.join(map(str,values))+'\n' for i in range(3)))
        fields.append({'name':name,'region':0,'mapping_status':'complete','support_kind':'node',
                       'components':len(values),'unit':unit,'csv_file':name+'.csv'})
        arrays.append(f'<DataArray Name="{name}" NumberOfComponents="{len(values)}">'+
                      ' '.join(map(str,values*3))+'</DataArray>')
        metadata.append(f'<DataArray type="String" Name="Unit:{name}">'+
                        ' '.join(str(ord(c)) for c in unit)+' 0</DataArray>')
    (root/'field_manifest.json').write_text(json.dumps({'fields':fields}))
    text=f'<VTKFile><UnstructuredGrid><FieldData><DataArray type="Float64" Name="Bias_V">{bias}</DataArray><DataArray type="Float64" Name="Converged">1</DataArray>'+''.join(metadata)+\
        '</FieldData><Piece><Points><DataArray NumberOfComponents="3">1 0 0 1 1 0 2 0 0</DataArray></Points><PointData>'+\
        ''.join(arrays)+'</PointData></Piece></UnstructuredGrid></VTKFile>'
    vtu = root/'stage01_candidate.vtu'; vtu.write_text(text)
    return mesh, vtu


class AvalancheQualificationTests(unittest.TestCase):
    def setUp(self):
        self.contract = a.read_json(a.CONTRACT)
        self.q = self.contract['biased_spatial_qualification']['pn2d_bv']
        deck = {'solver':{'temperature_K':300.,'impact_ionization':{'model':'van_overstraeten'}},
                'mesh_geometry':{'node_volume_policy':'mixed_voronoi','require_non_obtuse':True},
                'contacts':[{'name':'Anode','type':'ohmic','bias':0.}]}
        self.decks = [copy.deepcopy(deck), deck]
        self.q['stage_policy_sha256'] = a.stage_policy_sha256(deck)
        self.provenance = dict(self.q['input_sha256'])
        self.mesh = {'nodes':{'count':1943},'elements':{'count':3680},'qualification':{'non_obtuse':True}}

    def admit(self):
        return a.bv_spatial_qualification(self.contract,'pn2d_bv',self.provenance,
            self.q['reference_version'],self.q['reference_manifest_sha256'],self.decks,self.mesh)

    def test_exact_scope_admitted_but_static_limits_cannot_qualify_avalanche(self):
        self.assertEqual(self.admit()[0]['id'],self.q['id'])
        del self.contract['biased_spatial_qualification']
        qualification, reason = self.admit()
        self.assertIsNone(qualification)
        self.assertIn('No independent',reason)

    def test_input_reference_policy_and_mesh_drift_fail_closed(self):
        for key in self.provenance:
            previous = self.provenance[key]; self.provenance[key] = 'wrong'
            self.assertIsNone(self.admit()[0]); self.provenance[key] = previous
        for key in ('reference_version','reference_manifest_sha256'):
            args = [self.contract,'pn2d_bv',self.provenance,self.q['reference_version'],
                    self.q['reference_manifest_sha256'],self.decks,self.mesh]
            args[3 if key == 'reference_version' else 4] = 'wrong'
            self.assertIsNone(a.bv_spatial_qualification(*args)[0])
        for key, value in [('temperature_K',310.),('impact_ionization',{'model':'none'})]:
            original = copy.deepcopy(self.decks[1]['solver'])
            self.decks[1]['solver'][key] = value
            self.assertIsNone(self.admit()[0]); self.decks[1]['solver'] = original
        self.decks[1]['mesh_geometry']['node_volume_policy'] = 'barycentric'
        self.assertIsNone(self.admit()[0])
        self.decks[1]['mesh_geometry']['node_volume_policy'] = 'mixed_voronoi'
        self.mesh['qualification']['non_obtuse'] = False
        self.assertIsNone(self.admit()[0])

    def test_partial_pins_and_duplicate_anchors_are_rejected(self):
        del self.q['input_sha256']['parameter_sha256']
        self.assertIsNone(self.admit()[0])
        self.q['input_sha256'] = dict(self.provenance)
        junction_limit = self.q['spatial'].pop('junction_field_relative_l2')
        with self.assertRaisesRegex(ValueError,'incomplete independent'): self.admit()
        self.q['spatial']['junction_field_relative_l2'] = junction_limit
        self.q['anchors'].append(self.q['anchors'][0])
        with self.assertRaisesRegex(ValueError,'anchors'): self.admit()

    def test_independent_budget_rejects_potential_density_and_vector_errors(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); mesh, vtu = field_fixture(root)
            original = vtu.read_text(); limits = dict(self.q['spatial'],carrier_active_min_nodes=3)
            self.assertEqual(a.field_metrics(root,vtu,mesh,limits)['status'],'pass')
            for name, values in [('ElectrostaticPotential','0.002 0.002 0.002'),
                                 ('eDensity','1100000 1100000 1100000'),
                                 ('hDensity','1100000 1100000 1100000'),
                                 ('ElectricField','3.15 4.2 3.15 4.2 3.15 4.2')]:
                tree = ET.fromstring(original)
                tree.find(f'.//PointData/DataArray[@Name="{name}"]').text = values
                ET.ElementTree(tree).write(vtu)
                self.assertIn(name,a.field_metrics(root,vtu,mesh,limits)['failed_fields'])
            # The junction gate remains binding when the global budget permits it.
            limits['field_relative_l2'] = 1.
            self.assertIn('ElectricField',a.field_metrics(root,vtu,mesh,limits)['failed_fields'])
            vtu.write_text(original)
            with self.assertRaisesRegex(ValueError,'active carrier support'):
                a.field_metrics(root,vtu,mesh,dict(limits,carrier_active_min_nodes=4))

    def test_required_anchor_missing_duplicate_and_native_tamper_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); work = root/'work'; work.mkdir(); base = root/'reference'; (base/'tdr').mkdir(parents=True)
            mesh, vtu = field_fixture(work,-13.2)
            native = base/'tdr'/'anchor.tdr'; native.write_bytes(b'pinned native field')
            qualification = dict(self.q,spatial=dict(self.q['spatial'],carrier_active_min_nodes=3),
                anchors=[{'bias_V':-13.2,'native_tdr':'anchor.tdr','native_tdr_sha256':a.sha(native)}])
            def export(args, **kwargs):
                target = Path(args[-1]); target.mkdir(parents=True)
                shutil.copy2(work/'nodes.csv',target/'nodes.csv')
                shutil.copy2(work/'field_manifest.json',target/'field_manifest.json')
                shutil.copytree(work/'fields',target/'fields')
                return SimpleNamespace(returncode=0,stderr='')
            def measure(output, q=qualification):
                return a.biased_field_metrics('pn2d_bv',{'field_bias_V':-13.2},base,root/output,'importer',
                    work,mesh,{'cmd_sha256':None},q['spatial'],q)
            with patch.object(a.subprocess,'run',side_effect=export):
                self.assertEqual(measure('positive')['status'],'pass')
                extended = dict(qualification,required_biases_V=[-13.2,-20.],anchors=qualification['anchors']+[{'bias_V':-20.,'native_tdr':'anchor.tdr','native_tdr_sha256':a.sha(native)}])
                self.assertEqual(measure('missing',extended)['status'],'fail')
                diagnostic = dict(extended,required_biases_V=[-13.2])
                measured = measure('diagnostic',diagnostic)
                self.assertEqual(measured['status'],'pass')
                self.assertEqual(measured['diagnostic_status'],'fail')
                self.assertFalse(measured['anchors'][-1]['acceptance_required'])
                with self.assertRaisesRegex(ValueError,'original field bias'):
                    measure('excluded_primary',dict(extended,required_biases_V=[-20.]))
                shutil.copy2(vtu,work/'stage01_duplicate.vtu')
                self.assertEqual(measure('duplicate')['status'],'fail')
                (work/'stage01_duplicate.vtu').unlink()
                native.write_bytes(b'changed field')
                with self.assertRaisesRegex(ValueError,'fingerprint'): measure('tampered')


class AdmissionTests(unittest.TestCase):
    def test_numerical_sidecar_is_pinned_when_reusing_results(self):
        import tempfile
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            case = {'cmd':'pn2d_bv.cmd','grid':'mesh.tdr'}
            for name in ('pn2d_bv.cmd','mesh.tdr','models.par'):
                (root/name).write_bytes(name.encode())
            sidecar = root/'pn2d_bv.vela.json'
            sidecar.write_text('{"avalanche_current_support_profile":"element_edge_sg_gss_laux"}')
            provenance = {'cmd_sha256':a.sha(root/case['cmd']), 'grid_sha256':a.sha(root/case['grid']),
                          'parameter_sha256':a.sha(root/'models.par'), 'sidecar_sha256':a.sha(sidecar)}
            a.check_candidate_inputs(root,case,provenance)
            with self.assertRaisesRegex(ValueError,'sidecar'):
                a.check_candidate_inputs(root,case,{k:v for k,v in provenance.items() if k!='sidecar_sha256'})
            sidecar.write_text('{}')
            with self.assertRaisesRegex(ValueError,'sidecar'):
                a.check_candidate_inputs(root,case,provenance)
            sidecar.unlink()
            with self.assertRaisesRegex(ValueError,'sidecar'):
                a.check_candidate_inputs(root,case,provenance)

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

    def test_knee_gate_uses_exact_biases_and_rejects_premature_fitted_knee(self):
        settings = a.read_json(a.CONTRACT)['cases']['pn2d_bv']['breakdown_check']
        def curve(knee):
            return [{'bias_V':b,'I':-10**(-19+.1*abs(b)+.4*max(0.,abs(b)-knee))} for b in settings['biases_V']]
        reference = curve(19.6)
        self.assertEqual(a.breakdown_knee_metrics(reference,reference,settings)['status'],'pass')
        self.assertEqual(a.breakdown_knee_metrics(reference,curve(19.1),settings)['status'],'fail')
        with self.assertRaisesRegex(ValueError,'exact breakdown-knee bias'):
            a.breakdown_knee_metrics(reference,reference[:-1],settings)
        with self.assertRaisesRegex(ValueError,'duplicate'):
            a.breakdown_knee_metrics(reference,reference+[reference[0]],settings)

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
