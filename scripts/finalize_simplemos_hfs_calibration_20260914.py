"""Seal scoped HFS calibration without promoting it to production acceptance."""
from pathlib import Path
import check_simplemos_hfs_native_identity_20260912 as n

a,d,O,R=n.a,n.d,n.O,n.calibration.R


def main():
    assert not (O/'completion_20260914.json').exists()
    evidence=[]
    for stage in ('pilot','rest'):
        evidence += [O/f'{stage}_{name}_evidence.json' for name in ('native','export','identity','observer')]
        evidence += [O/folder/f'{stage}_evidence.json' for folder in
                     ('formation_20260914/v2','derivatives_20260914','chain_20260914','prior_floor_control_20260914')]
    probe=O/'cutoff_probe_20260914/v4'
    evidence += [probe/f'{name}_evidence.json' for name in ('input','native','export','integrity','analysis')]
    evidence += [O/'cutoff_probe_20260914/v3/integrity_evidence.json',O/'input_evidence.json']
    for e in evidence:a.verify(e)
    points=sum([a.rows(O/f'{s}_points.csv') for s in ('pilot','rest')],[])
    identities=sum([a.rows(O/f'{s}_identity.csv') for s in ('pilot','rest')],[])
    formation=sum([a.rows(O/'formation_20260914/v2'/f'{s}_summary.csv') for s in ('pilot','rest')],[])
    derivatives=[a.read(O/'derivatives_20260914'/f'{s}_summary.json') for s in ('pilot','rest')]
    chain=[a.read(O/'chain_20260914'/f'{s}_summary.json') for s in ('pilot','rest')]
    oldfloor=[a.read(O/'prior_floor_control_20260914'/f'{s}_summary.json') for s in ('pilot','rest')]
    final_probe=a.read(probe/'summary.json')
    assert len(points)==16 and all(r['qualified']=='True' for r in points)
    assert len(identities)==8 and all(r['qualified']=='True' for r in identities)
    assert all(r['all_local_partials_qualified'] for r in derivatives)
    assert all(r['all_prototype_columns_qualified'] for r in chain)
    assert final_probe['all_samples_qualified']
    selected=[r for r in formation if r['mode']=='boundary_any_contact_all_cut1']
    assert all(int(r['samples_without_field_match'])==0 for r in selected)
    enormal=R/'reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912'
    archive=a.read(enormal/'qualified_source_archive.json')
    changed=[name for name,sha in archive['files'].items() if a.sha(R/name)!=sha]
    assert not changed
    a.verify(enormal/'completion_evidence.json')
    base={r['name']:r for r in a.rows(n.calibration.O/'native_points.csv') if r['arm']=='baseline'}
    effects=[]
    for row in points:
        if row['arm']!='baseline':continue
        before=float(base[row['name']]['Id_A_per_um']);after=float(row['Id_A_per_um'])
        effects.append(dict(key=row['key'],device=row['device'],vd=row['vd'],vg=row['vg'],
            Enormal_native_Id=before,HFS_native_Id=after,native_model_change_percent=100*(after/before-1)))
    effect_path=O/'native_model_effect_20260914.csv'
    if effect_path.exists():
        assert a.rows(effect_path)==[{k:str(v) for k,v in r.items()} for r in effects]
    else:
        a.write_csv(effect_path,effects)
    result=dict(date='2026-09-14',completed_scope='Native HFS formation and scoped derivative calibration',
        native_DC_points=16,native_identity_pairs=8,
        max_kcl_over_Id=max(float(r['kcl_over_Id']) for r in points),
        native_observer_samples=sum(int(r['samples']) for r in selected),
        field_max_absolute_V_cm=max(float(r['max_sample_field_V_cm']) for r in selected),
        electron_cell_max_relative=max(float(r['max_cell_relative']) for r in selected if r['carrier']=='e'),
        hole_cell_max_relative=max(float(r['max_cell_relative']) for r in selected if r['carrier']=='h'),
        original_hole_cell_gate_passed=False,
        prior_diagnostic_floor_max_relative=max(r['left_max_relative'] for r in oldfloor),
        mean_first_electron_max_relative=max(float(r['max_mean_first_relative']) for r in selected if r['carrier']=='e'),
        scalar_derivative_samples=sum(r['samples'] for r in derivatives),
        scalar_derivative_max_relative=max(r['max_error'] for r in derivatives),
        composite_derivative_columns=sum(r['columns'] for r in chain),
        composite_derivative_max_relative=max(r['max_nonzero_relative'] for r in chain),
        cutoff_probe_samples=final_probe['samples'],cutoff_probe_eligible_nodes=final_probe['eligible_nodes'],
        cutoff_bracket_V_cm=final_probe['native_threshold_bracket_V_cm'],
        prior_v3_probe_loaded_state_failure_preserved=True,
        preserved_Enormal_files=len(archive['files']),preserved_Enormal_changed_files=changed,
        production_changed=False,production_HFS_accepted=False,
        native_matrix_derivatives_qualified=False,production_matrix_derivatives_qualified=False,
        full_HFS_curves_completed=False,new_fitted_parameters=0,acceptance_changed=False)
    a.write(O/'summary_20260914.json',result)
    report=R/'docs/validation/simplemos_hfs_native_formation_and_derivatives_2026-09-14.md'
    scripts=list((R/'scripts').glob('*hfs*20260914.py'))
    d.matrix.freeze(O/'completion_20260914.json',evidence+scripts+[
        report,O/'summary_20260914.json',O/'native_model_effect_20260914.csv',O/'preservation_20260914.json',
        enormal/'completion_evidence.json',enormal/'qualified_source_archive.json',
        R/'docs/validation/simplemos_constants_srh_response_validation_2026-09-07.md'])
    print(result,flush=True)


if __name__=='__main__':main()
