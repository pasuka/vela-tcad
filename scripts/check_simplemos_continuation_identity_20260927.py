"""Compare new same-bias controls with the frozen SRH-qualified binary results."""
import argparse
from pathlib import Path
import simplemos_hfs_cloud_20260926 as h

def compare(base,revision=''):
    old=h.read(base/'srh_validation/finite/summary/selected.json')
    root=base/('continuation_validation_20260927'+revision)
    new=h.read(root/'controls/finite/summary/selected.json')
    index={(r['case'],r['index'],r['arm']):r for r in old};out=[]
    for r in new:
        prior=index[r['case'],r['index'],r['arm']]
        geo=h.read(base/'replay_v4/inputs'/r['case']/'geometry.json')
        d=h.dual(prior,r,geo)
        a=h.ordered(Path(prior['dest'])/'state.csv',geo['count'])
        b=h.ordered(Path(r['dest'])/'state.csv',geo['count'])
        keys=[k for k in a[0] if k.startswith('packed_') or k in ('electron_qf_reference_V','hole_qf_reference_V')]
        d.update(case=r['case'],index=r['index'],arm=r['arm'],identical_coordinate_values=all(float(x[k])==float(y[k]) for x,y in zip(a,b) for k in keys))
        out.append(d)
    h.csvout(root/'same_bias_identity.csv',out)
    report=dict(states=len(out),qualified=sum(r['qualified'] for r in out),identical_coordinate_states=sum(r['identical_coordinate_values'] for r in out))
    h.write(root/'same_bias_identity.json',report);print(report)
    assert len(out)==32 and all(r['qualified'] for r in out)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--revision',default='');a=p.parse_args();compare(a.base,a.revision)
