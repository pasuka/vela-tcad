"""Independently recompute qualification and finite differences before sealing."""
from pathlib import Path
import ast
import copy
import math
import re
import subprocess
import numpy as np
import validate_simplemos_masetti_channel_20260907 as v

p=v.p;a=p.a;d=p.d;OUT=p.OUT;LOCAL=p.LOCAL


def same(x,y):
    assert math.isclose(x,y,rel_tol=1e-11,abs_tol=1e-26),(x,y)


def record_csv(path, rows):
    if path.exists():
        assert a.rows(path)==[{k:str(value) for k,value in row.items()} for row in rows]
    else:
        a.write_csv(path,rows)


def main():
    manifests=[p.prev.OUT/'validation_evidence.json']+[OUT/n for n in (
        'native_freeze.json','vela_freeze.json','analysis_scope.json','native_geometry_audit.json',
        'poisson_constants_addendum.json','free_row_screening_scope.json','native_parameter_audit.json','channel_response/freeze.json')]
    for path in manifests:a.verify(path)
    checks=[];fd=[];dc=a.rows(v.OUT/'dc.csv');cal=a.rows(v.OUT/'calibration.csv')
    assert len(dc)==20 and len(cal)==8
    for c in a.read(v.OUT/'contract.json')['cases']:
        geo=d.matrix.spatial.m73.Geometry(c['device']);mask=d.matrix.spatial.old.m78.supports(c['device'],geo,.05)[0]['all_si'].copy();mask[geo.contact_nodes]=False
        assert 2*int(sum(mask))==1814
        baseline=a.read(Path(c['base'])/'config.json');currents={}
        for job in c['jobs']:
            path=Path(job['config']);root=path.parent;cfg=a.read(path)
            x=copy.deepcopy(cfg);y=copy.deepcopy(baseline)
            for z in (x,y):
                z.pop('state_file');z.pop('output_state_file');z['solver'].pop('local_update_diagnostics')
            assert x==y
            assert cfg['state_file']==str(Path(c['base'])/'state.csv')
            post=a.read(root/'all_row.json');assert post['solver']['carrier_row_convergence']['mode']=='report'
            assert post['solver']['global_continuity_closure']==dict(mode='enforce',tolerance=1e-6,source_floor=1e-10)
            status=a.read(root/'config.status.json');q=a.read(root/'all_row.status.json')
            terms=d.ordered(root/'all_row.csv',geo.count);ratios=[]
            for row,keep in zip(terms,mask):
                if not keep:continue
                for carrier in ('electron','hole'):
                    scale=max(float(row[carrier+'_flux_abs_sum']),abs(float(row[carrier+'_recombination'])),abs(float(row[carrier+'_impact'])))
                    assert scale>0
                    ratios.append(abs(float(row[carrier+'_residual']))/scale)
            assert len(ratios)==q['carrier_row_convergence']['qualified_row_count']==1814
            bad=sum(r>1e-6 for r in ratios);assert bad==q['carrier_row_convergence']['violation_count']
            same(max(ratios),q['carrier_row_convergence']['max_ratio'])
            cc=status['contact_currents_A_per_um'];current=cc['drain'];kcl=abs(math.fsum(cc.values()))/abs(current)
            qualified=status['exit_code']==q['exit_code']==0 and status['converged'] and not bad and q['global_continuity_closure']['satisfied'] and kcl<=1e-8
            reported=next(r for r in dc if r['key']==c['key'] and r['label']==job['label'])
            assert qualified==(reported['qualified']=='True');same(current,float(reported['current_A_per_um']));assert qualified
            linear=a.rows(root/'linear_summary.csv')
            for iteration in {r['iteration'] for r in linear}:assert [int(r['correction']) for r in linear if r['iteration']==iteration]==list(range(5))
            currents[job['label']]=current
            checks.append(dict(key=c['key'],label=job['label'],active_rows=len(ratios),violations=bad,zero_scale_rows=0,
                max_row_ratio=max(ratios),kcl_over_Id=kcl,qualified=qualified,frozen_physics_and_solver_unchanged=True))
        root=v.LOCAL/c['key'];zero=d.ordered(root/'preflight/zero/residual.csv',geo.count);full=d.ordered(root/'preflight/full/residual.csv',geo.count)
        old=d.ordered(LOCAL/'vela'/c['key']/'strict/residual.csv',geo.count)
        for i,(z,f,o) in enumerate(zip(zero,full,old)):
            for key in ('psi_residual','phin_residual','phip_residual'):assert z[key]==o[key]
            for key in ('psi_residual','phip_residual'):assert z[key]==f[key]
            source=float(f['phin_residual'])-float(z['phin_residual'])
            expect=c['scaled_amplitude']*(int(i==320)-int(i==324))
            assert abs(source-expect)<=1e-8*abs(c['scaled_amplitude'])
        assert a.read(root/'preflight/zero/config.status.json')['current_A_per_um']==a.read(root/'preflight/full/config.status.json')['current_A_per_um']
        edge=next(e for e in a.rows(LOCAL/'vela'/c['key']/'strict/edges.csv') if int(e['edge_id'])==c['edge'])
        adj=d.ordered(LOCAL/'vela'/c['key']/'adjoint.csv',geo.count)
        conversion=d.fixed.Q*1e-6*float(edge['electron_particle_line_flux_per_m_s'])/float(edge['electron_flux'])
        prediction=-(float(adj[320]['lambda_electron'])-float(adj[324]['lambda_electron']))*c['current_amplitude_A_per_um']/conversion
        same(prediction,c['prediction_A_per_um'])
        same(c['native_amplitude']*d.fixed.Q*1e-8,c['current_amplitude_A_per_um'])
        center=currents['zero'];odd={n:(currents['plus_'+n]-currents['minus_'+n])/2 for n in ('full','half')}
        linearity=abs(odd['full']/(2*odd['half'])-1);drift=abs(center-c['base_Id_A_per_um'])
        for name,scale in [('full',1.),('half',.5)]:
            plus=currents['plus_'+name];minus=currents['minus_'+name]
            error=abs(odd[name]/(scale*prediction)-1);even=abs(((plus+minus)/2-center)/odd[name])
            snr=abs(odd[name])/max(drift,1e-300);driftdex=abs(math.log10(center/c['base_Id_A_per_um']))
            sign=(plus-center)*prediction>0 and (minus-center)*prediction<0
            passed=error<=.001 and linearity<=.001 and even<=.01 and snr>=100 and driftdex<=1e-5 and sign
            reported=next(x for x in cal if x['key']==c['key'] and x['amplitude']==name)
            # Reordering the independent prediction's products shifts a nearly
            # zero relative-error metric by a double ulp. Compare that derived
            # metric at machine precision, without altering its 0.001 gate.
            assert abs(error-float(reported['prediction_relative_error']))<=8*np.finfo(float).eps
            assert abs(linearity-float(reported['two_amplitude_relative']))<=8*np.finfo(float).eps
            assert passed==(reported['qualified']=='True');assert passed
            fd.append(dict(key=c['key'],amplitude=name,prediction_relative_error=error,two_amplitude_relative=linearity,
                even_over_odd=even,signal_over_zero_drift=snr,qualified=passed))
        print('Independently verified',c['key'],'5 states and 2 amplitudes',flush=True)
    # Probe status, exact edge-to-port identity, and diagonal model invariance.
    for c in a.read(OUT/'vela_contract.json')['cases']:
        for path in map(Path,c['probes']):assert a.read(path.with_suffix('.status.json'))['exit_code']==0
        assert a.read(LOCAL/'vela'/c['key']/'adjoint.status.json')['adjoint_relative_residual']<=1e-10
        cfg=a.read(Path(c['config']));mesh=a.read(Path(cfg['mesh_file']))
        for role in ('strict','mapped'):
            root=LOCAL/'vela'/c['key']/role;edges=a.rows(root/'edges.csv')
            cc=d.fixed.sum_contacts(edges,mesh,np.array([float(e['electron_particle_line_flux_per_m_s']) for e in edges],dtype=np.longdouble),np.array([float(e['hole_particle_line_flux_per_m_s']) for e in edges],dtype=np.longdouble))
            assert abs(cc['drain']/a.read(root/'functional.status.json')['current_A_per_um']-1)<=1e-8
    for path in (LOCAL/'bundle').rglob('*'):
        if path.is_file():assert a.sha(path)==a.sha(LOCAL/'native_raw/bundle'/path.relative_to(LOCAL/'bundle'))
    for device in ('n19','n23'):
        root=LOCAL/'native_raw/bundle'/device;assert (root/'exit_code.txt').read_text().strip()=='0'
        assert 'Good Bye' in (root/'console.log').read_text()
    assert len(a.rows(OUT/'native_parameter_check.csv'))==18
    assert all(r['exactly_equal']=='True' for r in a.rows(OUT/'native_parameter_check.csv'))
    assert all(float(r['native_constant_replay_over_charge'])<=1e-10 for r in a.rows(OUT/'native_poisson_constant_replay.csv'))
    assert all(r['passed']=='True' for r in a.rows(OUT/'native_element_support_check.csv'))
    assert all(float(r['projection_sum_relative_error'])<=1e-10 for r in a.rows(OUT/'free_row_screening_checks.csv'))
    # There is still no accepted native mobility identity; do not erase failures.
    assert not any(r['pointwise_identity_at_1e_minus_8']=='True' for r in a.rows(OUT/'element_mobility_candidates.csv'))
    record_csv(OUT/'independent_qualification.csv',checks);record_csv(OUT/'independent_channel_calibration.csv',fd)
    report=p.REPO/'docs/validation/simplemos_masetti_local_discretization_2026-09-07.md'
    scripts=[p.REPO/'scripts'/n for n in (
        'prepare_simplemos_masetti_local_20260907.py','analyze_simplemos_masetti_local_20260907.py',
        'validate_simplemos_masetti_channel_20260907.py','audit_simplemos_masetti_native_geometry_20260907.py',
        'audit_simplemos_masetti_poisson_constants_20260907.py','qualify_simplemos_masetti_free_row_ledger_20260907.py',
        'audit_simplemos_masetti_parameters_20260907.py','seal_simplemos_masetti_local_20260907.py')]
    for path in scripts:ast.parse(path.read_text(encoding='utf-8-sig'),filename=str(path))
    files=set(manifests+scripts+[report,p.REPO/'src/tools/vela_example_runner.cpp',
        p.REPO/'src/equation/CoupledDDAssembler.cpp',p.REPO/'src/physics/MobilityModel.cpp'])
    for root in (LOCAL,OUT):files.update(x for x in root.rglob('*') if x.is_file() and x.name!='validation_evidence.json')
    for target in re.findall(r'\]\(([^)]+)\)',report.read_text(encoding='utf-8')):
        if '://' in target:continue
        path=(report.parent/target.split('#')[0]).resolve()
        if path.name!='validation_evidence.json':assert path.exists(),path;files.add(path)
    # MSYS2 and the host Git differ in their CRLF config discovery. This read-only
    # comparison ignores line-ending whitespace, and actual file bytes are hashed.
    diff=subprocess.run(['git','-c','core.fsmonitor=false','diff','--numstat','--ignore-space-at-eol'],cwd=p.REPO,capture_output=True,text=True)
    assert diff.returncode==0 and diff.stdout.strip()=='127\t1\tsrc/tools/vela_example_runner.cpp'
    a.write(OUT/'validation_evidence.json',dict(schema_version=1,completed_date='2026-09-07',status='completed_bounded_masetti_local_audit',
        input_hashes={a.rel(x):a.sha(x) for x in sorted(files)},verified_manifests=len(manifests),
        native_additional_DC=2,native_parameter_exports=1,local_read_only_probes=56,source_preflight_probes=8,
        new_self_consistent_states=20,post_row_probes=20,independently_recomputed_rows=len(checks)*1814,
        qualified_states=20,qualified_channel_amplitudes=8,python_syntax_checked=len(scripts),
        observed_existing_conservative_flux_unit_tests=dict(command='python tests/regression/test_simplemos_local_conservative_flux.py',tests=5,result='passed',source='Executed in this task before sealing; tool output recorded in conversation.'),
        production_changes=False,acceptance_changes=False,m82_released=False,m83_released=False,
        limitations=['Eight operating points, four channel calibration points; no complete 0-1 V field-curve qualification.',
            'Native Masetti parameters match, but no tested cell mobility identity passed.',
            'Native Poisson expression was replayed, not its row residual directly exported.',
            'Distributed Poisson/transport projections remain screening; the channel calibration does not certify arbitrary sources.',
            'Raw contact-row component screening is retained but superseded by the explicit free-row ledger.',
            'No production candidate replacement or M82/M83 regression was released.']))
    a.verify(OUT/'validation_evidence.json');print('Sealed',len(files),'files; independently checked',len(checks)*1814,'carrier rows and 8 amplitude responses',flush=True)


if __name__=='__main__':main()
