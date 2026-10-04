"""Source and overlap qualification for immutable HFS curve extension."""
import argparse
from pathlib import Path
import simplemos_hfs_curves_20260914 as e
a,d=e.a,e.d

def overlap():
    a.verify(e.O/'native_evidence.json')
    current=a.rows(e.O/'native_points.csv'); rows=[]
    for cc in a.read(e.q.O/'inputs.json')['cases']:
        now=next(r for r in current if r['case']==cc['case'] and int(r['index'])==cc['index'])
        delta=abs(float(now['Id_A_per_um'])/cc['native_Id_A_per_um']-1)
        rows.append(dict(key=cc['key'],relative_Id=delta,qualified=now['native_qualified']=='True' and delta<=1e-12))
    a.write_csv(e.O/'native_overlap.csv',rows)
    d.matrix.freeze(e.O/'native_overlap_evidence.json',[Path(__file__).resolve(),e.O/'native_evidence.json',e.q.O/'inputs.json',e.O/'native_overlap.csv'])
    print(rows,flush=True);assert all(r['qualified'] for r in rows)

def freeze():
    e.gate();a.verify(e.O/'vela_freeze.json');a.verify(e.O/'native_overlap_evidence.json')
    assert not (e.O/'supplemental_source_evidence.json').exists()
    a.verify(e.q.O/'candidate_freeze.json')
    archive=e.q.L/'source_and_binary.tgz'
    assert a.sha(archive)=='a5a3b090dfe1cb86a2db017d389e98f1c0374158afd6d6c12e6f0a33062f1ff9'
    a.verify(e.O/'early_low_evidence.json')
    files=[e.O/'early_low_evidence.json',e.O/'vela_freeze.json',e.O/'native_overlap_evidence.json',e.q.O/'archive.json',archive]
    files+=list((e.R/'scripts').glob('*simplemos*.py'))
    files+=[f for directory in ('src','include') for f in (e.R/directory).rglob('*') if f.suffix in ('.cpp','.h')]
    d.matrix.freeze(e.O/'supplemental_source_evidence.json',files)
    print('Production binary, source, helper and input identities frozen.',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('overlap','freeze','controls','continuation','native','analyze'));args=p.parse_args();e.configure()
    if args.action in ('overlap','freeze'):globals()[args.action]()
    else:
        a.verify(e.O/'supplemental_source_evidence.json')
        if args.action=='controls':
            a.verify(e.O/'early_low_evidence.json')
            early=a.rows(e.O/'early_low_attempts.csv'); solve=e.v.solve_target
            def reuse(cc,index,arm,seed):
                if arm!='enormal_seed':return solve(cc,index,arm,seed)
                group=[r for r in early if r['case']==cc['case'] and int(r['index'])==index]
                assert group
                result=[]
                for rr in group:
                    root=Path(rr['dest']);a.verify(root/'input_freeze.json')
                    cfg=a.read(root/'config.json')
                    expected=a.read(Path(cc['template']))
                    for contact in expected['contacts']:
                        if contact['name']=='gate':contact['bias']=e.n.GRID[index]
                    actual={k:value for k,value in cfg.items() if k not in ('state_file','output_state_file')}
                    assert actual==expected, 'Early-control physics or numerical configuration differs'
                    if int(rr['attempt'])==0:assert cfg['state_file']==str(seed)
                    row=a.read(root/'result.json');assert row['arm']==arm and row['index']==index
                    result.append(row)
                return result
            e.v.solve_target=reuse
            e.controls()
        elif args.action=='analyze':e.v.analyze()
        else:
            a.verify(e.O/'control_evidence.json');assert a.read(e.O/'control_summary.json')['all_qualified'];e.v.run(args.action)
