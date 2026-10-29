from pathlib import Path
import argparse
import copy
import json
import math
import sys
import tempfile
import unittest
from unittest.mock import patch

from tests.regression import slot_state_fixture as fixture
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
import slot_reference_frame as frame
import run_slot_ldmos_reference_frame as entry
import slot_state_archive as states


class ReferenceFrameTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        mesh = dict(nodes=[dict(id=i, x=float(i%3), y=float(i//3)) for i in range(6)],
            triangles=[dict(id=0,node_ids=[0,1,2],region_id=0),dict(id=1,node_ids=[3,4,5],region_id=1)],
            regions=[dict(id=0,material='Si'),dict(id=1,material='SiO2')],
            contacts=[dict(name='source',node_ids=[0]),dict(name='drain',node_ids=[1]),dict(name='SLOT',node_ids=[4])])
        rows=[dict(node_id=i,psi=2.25,phin=2. if i<3 else 0.,phip=2. if i<3 else 0.,electrons_m3=1. if i<3 else 0.,holes_m3=2. if i<3 else 0.,
                   electron_qf_reference_V=2. if i<3 else 0.,electron_qf_increment_V=1e-22 if i<3 else 0.,
                   hole_qf_reference_V=2. if i<3 else 0.,hole_qf_increment_V=-1e-22 if i<3 else 0.) for i in range(6)]
        self.seed=self.root/'seed.h5'
        fixture.write(self.seed, rows, geometry=mesh, qualified=True,state_role='accepted',bias_contact='drain',bias_V=2.,contact_biases_V={'source':0.,'drain':2.,'SLOT':0.})
        self.doc=dict(mesh_file='mesh.json',simulation_type='dc_sweep',
            solver=dict(method='newton',abstol=1e-9,carrier_row_convergence=dict(mode='enforce'),global_continuity_closure=dict(mode='enforce')),
            contacts=[dict(name=n,type='ohmic',bias=0.) for n in ('source','drain','SLOT')],
            sweep=dict(contact='drain',bias_points=[2.,3.],start=2.,stop=3.,initial_step=.1,external_circuit=dict(enabled=False)))

    def convert(self, origin=1., name='out.h5', source=None, document=None):
        doc=document or self.doc;target=frame.translate_document(doc,origin);dest=self.root/name
        report=frame.translate_candidate(source or self.seed,dest,doc,self.root,target,2.)
        return dest,target,report

    def test_preserves_physical_voltage_tiny_qf_and_density_but_not_qualification(self):
        before=self.seed.read_bytes();dest,doc,r=self.convert()
        f,m=states.read(dest,doc,self.root);old,_=states.read(self.seed,self.doc,self.root)
        self.assertEqual(m['potential_origin_V'],1.)
        self.assertEqual(m['contact_biases_V'],{'source':-1.,'drain':1.,'SLOT':-1.})
        self.assertEqual(r['Vds_V'],2.);self.assertFalse(m['qualified'])
        for i in range(3):
            self.assertEqual(math.fsum((f['electron_qf_reference_V'][i],f['electron_qf_increment_V'][i],1.,-old['electron_qf_reference_V'][i],-old['electron_qf_increment_V'][i])),0.)
        self.assertEqual(f['phin'][4],-1.);self.assertEqual(f['phin'][3],0.)
        self.assertEqual(f['electron_qf_reference_V'][4],-1.)
        self.assertEqual(f['electron_qf_increment_V'][4],0.)
        self.assertEqual(r['inactive_contact_nodes'],1)
        self.assertEqual(f['electrons_m3'].tolist(),old['electrons_m3'].tolist())
        self.assertEqual(self.seed.read_bytes(),before)
        with self.assertRaises(FileExistsError):self.convert()

    def test_cross_frame_translation_and_immutable_physics(self):
        first,first_doc,_=self.convert(5.,'first.h5')
        second,second_doc,_=self.convert(1.,'second.h5',first,first_doc)
        f,m=states.read(second,second_doc,self.root)
        self.assertEqual(m['bias_V']+m['potential_origin_V'],2.)
        self.assertEqual(second_doc['solver'],self.doc['solver'])
        self.assertEqual(second_doc['sweep']['initial_step'],.1)
        self.assertEqual(second_doc['sweep']['bias_points'],[1.,2.])
        self.assertEqual(f['electron_qf_increment_V'][0],1e-22)

    def test_rejects_wrong_bias_origin_mesh_and_material(self):
        target=frame.translate_document(self.doc,1.)
        with self.assertRaisesRegex(ValueError,'terminal voltage'):
            frame.translate_candidate(self.seed,self.root/'bad.h5',self.doc,self.root,target,3.)
        with self.assertRaisesRegex(ValueError,'origin'):
            self.convert(document=dict(self.doc,potential_origin_V=4.))
        other=copy.deepcopy(target);other['mesh_file']='other.json'
        mesh=json.loads((self.root/'mesh.json').read_text());mesh['nodes'][0]['x']=50.
        (self.root/'other.json').write_text(json.dumps(mesh))
        with self.assertRaisesRegex(ValueError,'mesh identity'):
            frame.translate_candidate(self.seed,self.root/'bad.h5',self.doc,self.root,other,2.)
        for v in (float('nan'),float('inf'),True):
            with self.assertRaises(ValueError):frame.translate_document(self.doc,v)
        bad=copy.deepcopy(self.doc);bad['contacts'][0]['type']='metal_gate'
        with self.assertRaisesRegex(ValueError,'ohmic'):frame.translate_document(bad,1.)
        bad=copy.deepcopy(self.doc);bad['sweep']['external_circuit']=dict(enabled=True,mode='series_resistor')
        with self.assertRaisesRegex(ValueError,'direct'):frame.translate_document(bad,1.)

    def args(self):
        config=self.root/'base.json';config.write_text(json.dumps(self.doc))
        return argparse.Namespace(base_config=config,initial_state=self.seed,initial_drain_voltage_V=2.,
            target_origin_V=1.,previous_state=None,previous_drain_voltage_V=None,previous_source_origin_V=None,
            targets=[3.],output=self.root/'bundle',runner=self.root/'runner.exe',repo_root=self.root,
            max_newton_iterations=24,carrier_corrector_passes=4,max_attempts=48,min_step_V=.00625,timeout_seconds=300.)

    def test_preparation_never_claims_completed_and_freezes_nested_resources(self):
        geometry=self.root/'geometry.json';geometry.write_text('{}')
        self.doc['solver']['mobility']={'ialmob':{'avalanche_geometry_file':'geometry.json'}}
        args=self.args();doc,m=entry.prepare(args)
        self.assertFalse(m['completed']);self.assertFalse(m['seeds'][0]['qualified'])
        self.assertIn(str(geometry),m['resources'])
        self.assertEqual(doc['solver']['abstol'],1e-9)

    def test_failed_previous_anchor_prevents_current_and_target_dispatch(self):
        args=self.args();args.previous_state=self.seed;args.previous_drain_voltage_V=2.;args.initial_drain_voltage_V=3.
        # A separate current seed with matching provenance is required.
        fields,meta=states.read(self.seed,self.doc,self.root);meta.update(bias_V=3.,contact_biases_V={'source':0.,'drain':3.,'SLOT':0.})
        current=self.root/'current.h5';fixture.state_archive.write(current,fields,meta);args.initial_state=current
        args.targets=[4.];args.runner.write_bytes(b'runner')
        doc,m=entry.prepare(args)
        with patch.object(entry.subprocess,'run',return_value=argparse.Namespace(returncode=1)) as run:
            self.assertEqual(entry.execute(args,doc,m),1)
        self.assertEqual(run.call_count,1)
        self.assertFalse((args.output/'current').exists());self.assertFalse(m['completed'])
        command=run.call_args.args[0];index=command.index('--targets')
        self.assertEqual(command[index+1],'1.0')
        self.assertIn('--audit-local-continuity',command)

    def test_external_transport_resource_is_resolved_and_change_blocks_dispatch(self):
        couples=self.root/'couples.csv';couples.write_text('node0,node1,couple_m\n0,1,1e-7\n')
        self.doc['mesh_geometry']={'external_averagebox_couples_file':'couples.csv'}
        args=self.args();doc,m=entry.prepare(args)
        self.assertEqual(doc['mesh_geometry']['external_averagebox_couples_file'],str(couples))
        self.assertEqual(m['resources'][str(couples)],entry.sha(couples))
        couples.write_text('node0,node1,couple_m\n0,1,2e-7\n')
        with patch.object(entry.subprocess,'run') as run:
            with self.assertRaisesRegex(ValueError,'physical resource changed'):
                entry.execute(args,doc,m)
        run.assert_not_called()

    def test_success_reports_physical_voltage_and_requalifies_before_targets(self):
        args=self.args();args.runner.write_bytes(b'runner');doc,m=entry.prepare(args)
        def solve(command, **kwargs):
            first=command.index('--targets')+1
            self.assertEqual(command[first:first+2],['1.0','2.0'])
            folder=args.output/'current/reports';folder.mkdir()
            candidate=m['seeds'][0]['candidate']
            summary=dict(completed=True,last_successful_bias_V=2.,last_successful_state=candidate,
                attempts=[dict(converged=True,target_bias_V=v,current_A_per_um='1e-15',state_file=candidate) for v in (1.,2.)])
            (folder/'adaptive_noavalanche_summary.json').write_text(json.dumps(summary))
            return argparse.Namespace(returncode=0)
        with patch.object(entry.subprocess,'run',side_effect=solve):
            self.assertEqual(entry.execute(args,doc,m),0)
        self.assertTrue(m['completed'])
        self.assertEqual([p['physical_drain_voltage_V'] for p in m['phases'][0]['points']],[2.,3.])
        self.assertEqual([p['Vds_V'] for p in m['phases'][0]['points']],[2.,3.])

    def test_rejects_missing_provenance_and_nontransport_density(self):
        f,m=states.read(self.seed,self.doc,self.root)
        del m['contact_biases_V'];bad=self.root/'missing.h5';fixture.state_archive.write(bad,f,m)
        with self.assertRaisesRegex(ValueError,'provenance'):self.convert(source=bad)
        f,m=states.read(self.seed,self.doc,self.root);f['electrons_m3'][4]=1.
        bad=self.root/'density.h5';fixture.state_archive.write(bad,f,m)
        with self.assertRaisesRegex(ValueError,'nontransport'):self.convert(source=bad)

    def test_poly_option_requires_corrector_and_is_forwarded_explicitly(self):
        args=self.args();args.restore_poly_equilibrium=True;args.carrier_corrector_passes=0
        with self.assertRaisesRegex(ValueError,'require a carrier corrector'):entry.prepare(args)
        args.carrier_corrector_passes=4;args.runner.write_bytes(b'runner')
        with patch.object(entry.poly,'_profile') as profile, patch.object(entry.poly,'_components') as components:
            doc,m=entry.prepare(args)
        profile.assert_called_once();components.assert_called_once()
        self.assertTrue(m['restore_poly_equilibrium'])
        with patch.object(entry.subprocess,'run',return_value=argparse.Namespace(returncode=1)) as run:
            self.assertEqual(entry.execute(args,doc,m),1)
        self.assertIn('--restore-poly-equilibrium',run.call_args.args[0])
        self.assertFalse(m['completed'])

    def test_density_corrector_is_explicit_and_preserves_physical_configuration(self):
        args=self.args();args.carrier_corrector_mode='carrier_density';args.runner.write_bytes(b'runner')
        doc,m=entry.prepare(args)
        self.assertEqual(m['carrier_corrector_mode'],'carrier_density')
        self.assertEqual(doc['solver'],self.doc['solver'])
        with patch.object(entry.subprocess,'run',return_value=argparse.Namespace(returncode=1)) as run:
            self.assertEqual(entry.execute(args,doc,m),1)
        cmd=run.call_args.args[0]
        self.assertEqual(cmd[cmd.index('--carrier-corrector-mode')+1],'carrier_density')
        self.assertIn('--audit-local-continuity',cmd)
        self.assertFalse(m['completed'])

    def test_invalid_corrector_mode_is_rejected_before_preparation(self):
        args=self.args();args.carrier_corrector_mode='unknown'
        with self.assertRaisesRegex(ValueError,'corrector mode'):entry.prepare(args)
        self.assertFalse(args.output.exists())


if __name__=='__main__':unittest.main()
