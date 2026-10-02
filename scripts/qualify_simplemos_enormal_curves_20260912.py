"""Supplementary source identity and overlap gates for the Enormal curves."""
import argparse
from pathlib import Path
import simplemos_enormal_curves_20260912 as e
a,d=e.a,e.d


def overlap():
    a.verify(e.O/'native_evidence.json')
    old=a.rows(e.q.p.O/'native_points.csv')
    current=a.rows(e.O/'native_points.csv');rows=[]
    for cc in a.read(e.q.O/'inputs.json')['cases']:
        previous=next(r for r in old if r.get('name',r.get('job',''))==cc['key']+'_baseline')
        now=next(r for r in current if r['case']==cc['case'] and int(r['index'])==cc['index'])
        assert previous['qualified']=='True' and now['native_qualified']=='True'
        delta=abs(float(now['Id_A_per_um'])/float(previous['Id_A_per_um'])-1)
        rows.append(dict(key=cc['key'],relative_Id=delta,qualified=delta<=1e-12))
    a.write_csv(e.O/'native_overlap.csv',rows)
    d.matrix.freeze(e.O/'native_overlap_evidence.json',[Path(__file__).resolve(),e.O/'native_evidence.json',e.q.p.O/'native_points.csv',e.O/'native_overlap.csv'])
    print(rows,flush=True);assert all(r['qualified'] for r in rows)


def freeze():
    e.gate();a.verify(e.O/'vela_freeze.json');a.verify(e.O/'native_overlap_evidence.json')
    assert not (e.O/'supplemental_source_evidence.json').exists()
    archived=a.read(e.q.O/'archive_evidence.json')
    assert a.sha(e.q.L/'candidate_source_and_binary.tgz')==archived['sha256']
    files=[e.O/'vela_freeze.json',e.O/'native_overlap_evidence.json',e.q.O/'archive_evidence.json',e.q.L/'candidate_source_and_binary.tgz']
    files+=list((e.R/'scripts').glob('*simplemos*.py'))
    files+=[f for directory in ('src','include') for f in (e.R/directory).rglob('*') if f.suffix in ('.cpp','.h')]
    d.matrix.freeze(e.O/'supplemental_source_evidence.json',files)
    print('Supplemental source, helper and candidate archive identities frozen.',flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=('overlap','freeze','controls','continuation','native','analyze'));args=ap.parse_args();e.configure()
    if args.action in ('overlap','freeze'):globals()[args.action]()
    else:
        a.verify(e.O/'supplemental_source_evidence.json')
        if args.action=='controls':e.controls()
        elif args.action=='analyze':e.v.analyze()
        else:
            a.verify(e.O/'control_evidence.json');assert a.read(e.O/'control_summary.json')['all_qualified']
            e.v.run(args.action)
