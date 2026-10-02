"""Report all reduced-model outcomes without promoting failed states or changing gates."""
from pathlib import Path
import math
import numpy as np
import run_simplemos_reference_subset_20260906 as r

a=r.a;d=r.d;OUT=r.OUT;LOCAL=r.LOCAL


def yes(value):
    return value is True or value=='True'


def main():
    a.verify(OUT/'vela_freeze.json');a.verify(OUT/'export_contract.json')
    rows=a.rows(OUT/'vela_runs.csv');assert len(rows)==24
    keyed={(x['arm'],x['device'],float(x['vd']),float(x['vg'])):x for x in rows}
    assert len(keyed)==24
    summary={}
    for arm in r.ARMS:
        xs=[x for x in rows if x['arm']==arm];good=[x for x in xs if yes(x['comparison_qualified'])]
        summary[arm]=dict(total=len(xs),qualified=len(good),failed=len(xs)-len(good),
            qualified_max_abs_Id_error_relative=max((abs(float(x['signed_Id_error_relative'])) for x in good),default=None),
            diagnostic_all_max_abs_Id_error_relative=max(abs(float(x['signed_Id_error_relative'])) for x in xs),
            failure_reasons=sorted(set(x['failure'] for x in xs if not yes(x['comparison_qualified']))))
    pairs=[];models=[];numerics=[]
    for vd in (.05,1.):
        for vg in (.8,1.):
            for arm in r.ARMS:
                lo=keyed[arm,'n19',vd,vg];hi=keyed[arm,'n23',vd,vg]
                native=math.log10(float(hi['native_Id_A_per_um'])/float(lo['native_Id_A_per_um']))
                vela=math.log10(float(hi['Id_A_per_um'])/float(lo['Id_A_per_um']))
                pairs.append(dict(arm=arm,vd=vd,vg=vg,qualified=yes(lo['comparison_qualified']) and yes(hi['comparison_qualified']),
                    native_log10_high_over_low=native,vela_log10_high_over_low=vela,pair_error_dex=vela-native))
            for device in ('n19','n23'):
                b=keyed['baseline_refined',device,vd,vg];m=keyed['masetti_refined',device,vd,vg];p=keyed['masetti_plain',device,vd,vg]
                models.append(dict(device=device,vd=vd,vg=vg,baseline_qualified=yes(b['comparison_qualified']),masetti_qualified=yes(m['comparison_qualified']),
                    qualified=yes(b['comparison_qualified']) and yes(m['comparison_qualified']),
                    baseline_error_relative=float(b['signed_Id_error_relative']),masetti_error_relative=float(m['signed_Id_error_relative']),
                    absolute_error_change_percentage_points=100*(abs(float(m['signed_Id_error_relative']))-abs(float(b['signed_Id_error_relative']))),
                    vela_Id_model_ratio=float(m['Id_A_per_um'])/float(b['Id_A_per_um']),
                    native_Id_model_ratio=float(m['native_Id_A_per_um'])/float(b['native_Id_A_per_um'])))
                plain=d.ordered(Path(p['config']).parent/'state.csv',d.matrix.spatial.m73.Geometry(device).count)
                refined=d.ordered(Path(m['config']).parent/'state.csv',len(plain))
                delta={f'{k}_max_absolute_delta_V':float(max(abs(float(x[k])-float(y[k])) for x,y in zip(plain,refined))) for k in ('psi','phin','phip')}
                numerics.append(dict(device=device,vd=vd,vg=vg,plain_qualified=yes(p['comparison_qualified']),refined_qualified=yes(m['comparison_qualified']),
                    both_qualified=yes(p['comparison_qualified']) and yes(m['comparison_qualified']),plain_iterations=int(p['iterations']),refined_iterations=int(m['iterations']),
                    Id_relative_change=float(m['Id_A_per_um'])/float(p['Id_A_per_um'])-1,**delta))
    a.write_csv(OUT/'nwell_pairs.csv',pairs);a.write_csv(OUT/'model_comparison.csv',models);a.write_csv(OUT/'numerical_comparison.csv',numerics)
    a.write(OUT/'summary.json',dict(arms=summary,native_states=16,native_qualified=sum(yes(x['native_qualified']) for x in a.rows(OUT/'native_points.csv')),
        strict_qualified_model_pairs=sum(x['qualified'] for x in models),qualified_nwell_pairs=sum(x['qualified'] for x in pairs),
        scope='Eight independent working points; composite bulk/surface/high-field control. Failed outcomes remain diagnostic; no full-curve or initialization-invariance qualification.',
        production_changes=False,acceptance_changes=False,m82_released=False,m83_released=False))
    print(summary,flush=True)


if __name__=='__main__':main()
