"""Qualify the native zero parameter controls before the signed response batch."""
import argparse,math,tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import simplemos_enormal_response_20260912 as r
p=r.p;a,d=r.a,r.d;L,O=r.L,r.O

def collect(stage,sha):
    a.verify(O/'native_freeze.json');archive=L/f'{stage}_results.tgz';assert a.sha(archive)==sha
    raw=L/f'{stage}_raw';assert not raw.exists()
    with tarfile.open(archive) as t:
        for m in t.getmembers():assert (raw/m.name).resolve().is_relative_to(raw.resolve()) and not m.issym() and not m.islnk()
        t.extractall(raw,filter='data')
    for f in (L/'bundle').rglob('*'):
        if f.is_file():assert a.sha(f)==a.sha(raw/'bundle'/f.relative_to(L/'bundle'))
    jobs=[j for j in a.read(O/'native_contract.json')['jobs'] if stage!='zero' or j['label']=='zero'];rows=[];exports=[]
    for j in jobs:
        root=raw/'bundle'/j['name'];code=int((root/'exit_code.txt').read_text());log=(root/'console.log').read_text(errors='replace')
        values=p.c.n.exporter.pltrows(root/'native_des.plt');assert len(values)==1;v=values[0]
        currents=[v[k+' TotalCurrent'] for k in ('drain','source','gate','substrate')];Id=currents[0]
        kcl=abs(math.fsum(currents))/abs(Id);bias=max(abs(v['gate OuterVoltage']-j['vg']),abs(v['drain OuterVoltage']-j['vd']))
        good=code==0 and 'T-2022.03-SP2' in log and 'Good Bye' in log and kcl<=1e-8 and bias<=1e-10
        rows.append(dict(**j,Id_A_per_um=Id,kcl_over_Id=kcl,bias_error_V=bias,native_qualified=good));assert good,rows[-1]
        exports.append(dict(case=j['case'],index=j['index'],tdr=str(root/'final_des.tdr'),export=str(L/(stage+'_exports')/j['name'])))
    a.write_csv(O/f'{stage}_points.csv',rows)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(p.c.n.exporter.export_one,exports))
    if stage=='zero':
        old={j['name']:j for j in a.rows(p.O/'native_points.csv')};identities=[];checks=[]
        for j in rows:
            src=L/'zero_exports'/j['name'];base=p.L/'exports'/j['baseline'];ident=p.identity(base,src)
            identities.extend(dict(name=j['name'],**v) for v in ident)
            potential=max(v['max_absolute'] for v in ident if 'Potential' in v['field']);density=max(v['max_relative'] for v in ident if v['field'].startswith(('eDensity_','hDensity_')))
            mobility=max(v['max_relative'] for v in ident if v['field'].startswith(('eMobility_','hMobility_')))
            current=abs(j['Id_A_per_um']/float(old[j['baseline']]['Id_A_per_um'])-1)
            checks.append(dict(name=j['name'],current_relative=current,potential_max_V=potential,density_relative=density,mobility_relative=mobility,qualified=current<=1e-12 and potential<=1e-12 and density<=1e-10 and mobility<=1e-12))
        a.write_csv(O/'zero_identity.csv',checks);a.write_csv(O/'zero_fields.csv',identities);assert all(j['qualified'] for j in checks),checks
    else:
        zero={j['case']:j for j in a.rows(O/'zero_identity.csv')} if False else {j['name']:j for j in a.rows(O/'zero_identity.csv')}
        checks=[]
        for case in sorted({j['case'] for j in rows}):
            subset={j['label']:j for j in rows if j['case']==case};I0=subset['zero']['Id_A_per_um'];derivatives=[]
            for label,amp in [('large',.001),('small',.0005)]:
                ip=subset['plus_'+label]['Id_A_per_um'];im=subset['minus_'+label]['Id_A_per_um'];derivative=(ip-im)/(2*amp);derivatives.append(derivative)
                odd=abs(ip-im)/2;even=abs((ip-I0)+(im-I0))/2;drift=zero[subset['zero']['name']]['current_relative'];drift=float(drift)*abs(I0)
                checks.append(dict(case=case,amplitude=amp,Id_zero=I0,derivative_A_per_um=derivative,normalized_derivative=derivative/I0,even_over_odd=even/max(odd,1e-300),signal_over_zero_drift=odd/drift if drift else None,qualified=even<=.01*odd and odd>=100*drift))
            convergence=abs(derivatives[1]/derivatives[0]-1)
            for j in checks[-2:]:j['two_amplitude_relative']=convergence;j['qualified'] &= convergence<=1e-3
        a.write_csv(O/'native_response.csv',checks);assert all(j['qualified'] for j in checks),checks
    files=[Path(__file__).resolve(),O/'native_freeze.json',archive,O/f'{stage}_points.csv']+[f for f in (L/(stage+'_exports')).rglob('*') if f.is_file()]
    files += [O/n for n in (('zero_identity.csv','zero_fields.csv') if stage=='zero' else ('native_response.csv',))]
    d.matrix.freeze(O/f'{stage}_evidence.json',files)
    print(checks,flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=('zero','signed'));parser.add_argument('--sha',required=True);args=parser.parse_args();collect(args.stage,args.sha)
