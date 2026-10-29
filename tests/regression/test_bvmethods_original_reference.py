import copy
from pathlib import Path
import sys
import unittest
import hashlib
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from validate_bvmethods_original_reference import source_checks, stop_record_checks, validate, STOP
from validate_wp14_zero_qs import BASELINE, trajectory_checks
from vela_sdevice.cmd_parser import parse_cmd_text


def state(current, inner):
    row = {'time': (inner+1e7*current)/100000}
    for c in ['drain','source','gate','substrate']:
        i = current if c=='drain' else -current if c=='source' else 0
        row[c] = dict(total_current_A_per_um=i,e_current_A_per_um=i,
                      h_current_A_per_um=0,DisplacementCurrent=0,
                      outer_voltage_V=inner+1e7*current if c=='drain' else 0,
                      inner_voltage_V=inner if c=='drain' else 0)
    return row


class OriginalReferenceTest(unittest.TestCase):
    def setUp(self):
        # Repository fixture transformed explicitly to the distinct profile;
        # no commercial source or generated outputs are checked in here.
        self.old = (BASELINE/'external_resistor_sdevice.cmd').read_text()
        self.source = self.old.replace('Band2Band(E2)','Band2Band(Model=NonLocalPath)').replace('Avalanche(Eparallel)','Avalanche(GradQuasiFermi)').replace('1.443e-3','1.3e-3')
        self.rows = [state(0,0),state(1e-4,6.6),state(1.35e-3,7.3)]
        self.log = 'Contact current exceeded limit:\nelectrode: drain\nvalue: 1.3500e-03\nlimit: 1.3000e-03\nGood Bye !'

    def test_original_profile_and_equivalent_numeric_spelling(self):
        self.assertTrue(all(source_checks(self.source).values()))
        source=self.source.replace('1.3e-3','0.001300').replace('1e-7','0.0000001')
        self.assertTrue(all(source_checks(source).values()))

    def test_local_E2_or_alternative_avalanche_cannot_pass(self):
        for replacement in [('NonLocalPath','E2'),('Avalanche(GradQuasiFermi)','Avalanche(Eparallel)')]:
            self.assertFalse(source_checks(self.source.replace(*replacement))['original_physics'])
        self.assertFalse(all(source_checks(self.old).values()))

    def test_changed_stop_and_resistor_are_distinct_contracts(self):
        self.assertFalse(source_checks(self.source.replace('1.3e-3','1.443e-3'))['original_math'])
        self.assertFalse(source_checks(self.source.replace('Resistor=1e7','Resistor=1e6'))['original_electrode'])

    def test_Load_or_auxiliary_gate_cannot_prove_from_zero(self):
        for addition in ['Load(FilePrefix="seed")','Quasistationary(Goal{Name="gate" Voltage=.2}){Coupled{Poisson}}']:
            self.assertFalse(source_checks(self.source.replace('Solve {','Solve {\n'+addition))['original_solve'])

    def test_duplicate_physics_or_extra_precision_not_original(self):
        self.assertFalse(source_checks(self.source+'\nPhysics(Material="Silicon"){Fermi}')['single_device_sections'])
        self.assertFalse(source_checks(self.source.replace('Math {','Math { ExtendedPrecision'))['original_math'])

    def test_comments_do_not_activate_diagnostic_features(self):
        self.assertTrue(all(source_checks('* Load() ExtendedPrecision\n'+self.source).values()))

    def test_original_stop_is_explicit_and_does_not_change_old_default(self):
        self.assertEqual(trajectory_checks(self.rows,self.log,stop=STOP)['status'],'pass')
        self.assertFalse(trajectory_checks(self.rows,self.log)['checks']['original_stop_reached'])

    def test_original_stop_report_is_bound_to_last_point(self):
        self.assertTrue(all(stop_record_checks(self.log,.00135).values()))
        self.assertFalse(stop_record_checks(self.log.replace('1.3000e-03','1.4430e-03'),.00135)['declared_stop_matches'])
        self.assertFalse(stop_record_checks(self.log,.0014)['reported_value_matches_last'])
        self.assertFalse(stop_record_checks('Good Bye !',.00135)['single_drain_stop_record'])

    def test_intermediate_load_or_component_failure_is_not_hidden(self):
        for key in ['outer_voltage_V','h_current_A_per_um']:
            rows=copy.deepcopy(self.rows);rows[1]['drain'][key]+=.1
            self.assertEqual(trajectory_checks(rows,self.log,stop=STOP)['status'],'fail')

    def test_nonlocal_path_is_explicitly_unsupported_by_frontend(self):
        ir=parse_cmd_text(self.source,path_name='original.cmd')
        entries=[x for x in ir['physics']['effective'] if x['key']=='Recombination.Band2Band']
        self.assertEqual(len(entries),1)
        self.assertEqual(entries[0]['status'],'unsupported')
        self.assertTrue(any(x['id']=='E106' for x in ir['errors']))
        self.assertIsNone(entries[0].get('vela'))

    def test_bound_input_files_and_parameter_tampering(self):
        parent=Path(__file__).resolve().parents[2]/'build'/'test_bvmethods_original_reference'
        parent.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(dir=parent) as tmp:
            root=Path(tmp);cmd=root/'case.cmd';par=root/'models.par';mesh=root/'bvmethods_nmos_msh.tdr'
            plt=root/'external_resistor_des.plt';log=root/'external_resistor_des.log'
            refcmd=root/'frozen.cmd';refpar=root/'frozen.par'
            cmd.write_text(self.source);refcmd.write_text(self.source)
            par.write_text('Scharfetter {taumax=1e-7,1e-7}');refpar.write_bytes(par.read_bytes())
            mesh.write_bytes(b'hash-bound test mesh')
            functions=['OuterVoltage','InnerVoltage','eCurrent','hCurrent','TotalCurrent','DisplacementCurrent']
            keys=['outer_voltage_V','inner_voltage_V','e_current_A_per_um','h_current_A_per_um','total_current_A_per_um','DisplacementCurrent']
            contacts=['drain','source','gate','substrate']
            names=['time']+[c+' '+f for c in contacts for f in functions]
            values=[[r['time']]+[r[c][k] for c in contacts for k in keys] for r in self.rows]
            plt.write_text('DF-ISE text\nInfo{datasets=['+' '.join('"'+n+'"' for n in names)+']}\nData{\n'+'\n'.join(' '.join(map(str,row)) for row in values)+'\n}')
            log.write_text('Version T-2022.03-SP2\n'+self.log)
            sha=hashlib.sha256(mesh.read_bytes()).hexdigest()
            args=(cmd,par,mesh,plt,log,refcmd,refpar,sha)
            self.assertEqual(validate(*args)['status'],'pass')
            par.write_text('Scharfetter {taumax=2e-7,1e-7}')
            self.assertFalse(validate(*args)['source_and_binding_checks']['parameters_match_frozen_reference'])
            par.write_bytes(refpar.read_bytes());mesh.write_bytes(b'other mesh')
            self.assertFalse(validate(*args)['source_and_binding_checks']['mesh_matches_frozen_hash'])
            mesh.write_bytes(b'hash-bound test mesh');log.write_text('Version O-2018.06-SP2\n'+self.log)
            self.assertFalse(validate(*args)['source_and_binding_checks']['native_version_matches'])
            log.write_text('Version T-2022.03-SP2\n'+self.log);cmd.write_text(self.source.replace('Current   = "external_resistor_des.plt"','Current   = "other.plt"'))
            self.assertFalse(validate(*args)['source_and_binding_checks']['file_current_bound'])


if __name__=='__main__':
    unittest.main()
