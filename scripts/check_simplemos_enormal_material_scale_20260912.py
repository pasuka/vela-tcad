"""Qualify corrected native material scaling, rejecting an unresolvable signal."""
import argparse, math, tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import simplemos_enormal_material_scale_20260912 as r

a,d,p,L,O=r.a,r.d,r.p,r.L,r.O

def collect(stage, expected):
    a.verify(O/'native_freeze.json')
    archive=L/f'{stage}_results.tgz'
    assert a.sha(archive)==expected
    raw=L/f'{stage}_raw'; assert not raw.exists()
    with tarfile.open(archive) as t:
        for m in t.getmembers():
            assert (raw/m.name).resolve().is_relative_to(raw.resolve()) and not m.issym() and not m.islnk()
        t.extractall(raw,filter='data')
    contract=a.read(O/'native_contract.json')
    jobs=[j for j in contract['jobs'] if (j['case']==contract['pilot_case'])==(stage=='pilot')]
    rows=[]; exports=[]
    for j in jobs:
        root=raw/'bundle'/j['name']
        for f in (L/'bundle'/j['name']).iterdir(): assert a.sha(f)==a.sha(root/f.name)
        code=int((root/'exit_code.txt').read_text());log=(root/'console.log').read_text(errors='replace')
        values=p.c.n.exporter.pltrows(root/'native_des.plt'); assert len(values)==1
        v=values[0]; currents=[v[k+' TotalCurrent'] for k in ('drain','source','gate','substrate')]
        Id=currents[0]; assert math.isfinite(Id) and Id!=0
        kcl=abs(math.fsum(currents))/abs(Id)
        bias=max(abs(v['gate OuterVoltage']-j['vg']),abs(v['drain OuterVoltage']-j['vd']))
        good=code==0 and 'T-2022.03-SP2' in log and 'Good Bye' in log and kcl<=1e-8 and bias<=1e-10
        rows.append(dict(**j,exit_code=code,Id_A_per_um=Id,kcl_over_Id=kcl,bias_error_V=bias,native_qualified=good))
        exports.append(dict(case=j['case'],index=j['index'],tdr=str(root/'final_des.tdr'),export=str(L/(stage+'_exports')/j['name'])))
    a.write_csv(O/f'{stage}_points.csv',rows)
    assert all(j['native_qualified'] for j in rows), rows
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(p.c.n.exporter.export_one,exports))
    old={j['name']:j for j in a.rows(p.O/'native_points.csv')}
    zero_checks=[];identities=[];checks=[]
    for case in sorted({j['case'] for j in rows}):
        subset={j['label']:j for j in rows if j['case']==case};z=subset['zero']
        ident=p.identity(p.L/'exports'/z['baseline'],L/(stage+'_exports')/z['name'])
        identities.extend(dict(name=z['name'],**v) for v in ident)
        potential=max(v['max_absolute'] for v in ident if 'Potential' in v['field'])
        density=max(v['max_relative'] for v in ident if v['field'].startswith(('eDensity_','hDensity_')))
        mobility=max(v['max_relative'] for v in ident if v['field'].startswith(('eMobility_','hMobility_')))
        I0=z['Id_A_per_um']; oldI=float(old[z['baseline']]['Id_A_per_um']);drift=abs(I0-oldI)
        current=abs(I0/oldI-1)
        zero_checks.append(dict(case=case,current_relative=current,potential_max_V=potential,density_relative=density,mobility_relative=mobility,qualified=current<=1e-12 and potential<=1e-12 and density<=1e-10 and mobility<=1e-12))
        derivatives=[]; local=[]
        for label,amp in [('large',.001),('small',.0005)]:
            ip=subset['plus_'+label]['Id_A_per_um'];im=subset['minus_'+label]['Id_A_per_um']
            derivative=(ip-im)/(2*amp);derivatives.append(derivative)
            odd=abs(ip-im)/2;even=abs((ip-I0)+(im-I0))/2
            resolvable=odd>0 and odd>=100*drift and odd>=100*math.ulp(I0)
            local.append(dict(case=case,amplitude=amp,Id_zero=I0,derivative_A_per_um=derivative,normalized_derivative=derivative/I0,odd_A_per_um=odd,even_over_odd=even/odd if odd else None,signal_over_zero_drift=odd/drift if drift else None,nonzero_resolvable_signal=resolvable,qualified=resolvable and even<=.01*odd and zero_checks[-1]['qualified']))
        convergence=abs(derivatives[1]/derivatives[0]-1) if derivatives[0]!=0 else None
        for j in local:
            j['two_amplitude_relative']=convergence;j['qualified'] &= convergence is not None and convergence<=1e-3
        checks.extend(local)
    a.write_csv(O/f'{stage}_zero_identity.csv',zero_checks);a.write_csv(O/f'{stage}_zero_fields.csv',identities)
    a.write_csv(O/f'{stage}_native_response.csv',checks)
    summary=dict(stage=stage,native_points=len(rows),zero_passed=sum(j['qualified'] for j in zero_checks),response_checks=len(checks),response_passed=sum(j['qualified'] for j in checks),all_qualified=all(j['qualified'] for j in checks),acceptance_relaxed=False)
    a.write(O/f'{stage}_summary.json',summary)
    files=[Path(__file__).resolve(),O/'native_freeze.json',archive]+[O/f'{stage}_{n}' for n in ('points.csv','zero_identity.csv','zero_fields.csv','native_response.csv','summary.json')]
    files += [f for f in (L/(stage+'_exports')).rglob('*') if f.is_file()]
    d.matrix.freeze(O/f'{stage}_evidence.json',files)
    print(summary,checks,flush=True);assert summary['all_qualified']

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=('pilot','rest'));ap.add_argument('--sha',required=True)
    args=ap.parse_args();collect(args.stage,args.sha)
