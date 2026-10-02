"""Check actual Sentaurus 2022 default parameter export against Vela values."""
from pathlib import Path
import re
import prepare_simplemos_masetti_local_20260907 as p
import audit_simplemos_subset_mobility_export_20260907 as m


def main():
    path=p.LOCAL/'native_models.par';text=path.read_text()
    blocks={name:re.search(r'(?m)^'+name+r':\s*\{(.*?)\n\}',text,re.S).group(1) for name in ('ConstantMobility','DopingDependence')}
    def values(block,key):
        value=re.search(r'(?m)^\s*'+key+r'\s*=\s*([^#\n]+)',blocks[block]).group(1)
        return tuple(float(x.strip()) for x in value.split(','))
    assert values('DopingDependence','formula')==(1.,1.)
    rows=[]
    for index,carrier in enumerate(('e','h')):
        native=(values('ConstantMobility','mumax')[index],)+tuple(values('DopingDependence',k)[index] for k in ('mumin1','mumin2','mu1','Pc','Cr','Cs','alpha','beta'))
        assert native==m.PARAMETERS[carrier]
        for key,x,y in zip(('mumax','mumin1','mumin2','mu1','Pc','Cr','Cs','alpha','beta'),native,m.PARAMETERS[carrier]):
            rows.append(dict(carrier=carrier,parameter=key,native=x,vela=y,exactly_equal=x==y))
    p.a.write_csv(p.OUT/'native_parameter_check.csv',rows)
    p.a.write(p.OUT/'native_parameter_audit.json',dict(status='all_18_Masetti_parameters_match',
        command='sdevice -P:Silicon',remote_directory=p.REMOTE+'/parameter_audit',
        formula_electron=1,formula_hole=1,temperature_K=300,
        limitations='This checks built-in numerical parameters for the selected formula. It does not identify internal vertex/element averaging or prove an exported plot is native edge mobility.',
        input_hashes={p.a.rel(x):p.a.sha(x) for x in [Path(__file__).resolve(),path,Path(m.__file__).resolve()]}))
    print('18 native Masetti parameters exactly match Vela; both select formula 1',flush=True)


if __name__=='__main__':main()
