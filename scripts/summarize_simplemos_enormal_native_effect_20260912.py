"""Native-only finite model restoration effect; never a Vela parity result."""
from pathlib import Path
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import simplemos_enormal_curves_20260912 as e
a,d=e.a,e.d
O=e.O/'native_restoration_effect'


def main():
    a.verify(e.PHUMOB/'completion_evidence.json');a.verify(e.O/'native_evidence.json')
    assert not (O/'evidence.json').exists()
    old=a.rows(e.PHUMOB/'native_points.csv');new=a.rows(e.O/'native_points.csv');rows=[];checks=[]
    files=[Path(__file__).resolve(),e.PHUMOB/'completion_evidence.json',e.O/'native_evidence.json']
    norm=lambda s:re.sub(r'\s+','',s)
    for cc in a.read(e.O/'native_contract.json')['cases']:
        paths=[e.R/'build-release/phumob_curves_20260912/native_raw/bundle'/cc['case'],e.L/'native_raw/bundle'/cc['case']]
        texts=[(p/'native_des.cmd').read_text() for p in paths]
        physics=[t.split('Physics',1)[1].split('Plot',1)[0] for t in texts]
        math=[re.search(r'Math\s*\{[^}]*\}',t,re.S).group() for t in texts]
        identity=dict(case=cc['case'],mesh_identical=a.sha(paths[0]/'input_fps.tdr')==a.sha(paths[1]/'input_fps.tdr'),other_physics_identical=norm(physics[0])==norm(physics[1].replace('Mobility(PhuMob Enormal)','Mobility(PhuMob)')),math_identical=norm(math[0])==norm(math[1]))
        assert all(identity[k] for k in ('mesh_identical','other_physics_identical','math_identical'));checks.append(identity)
        files.extend(p/f for p in paths for f in ('input_fps.tdr','native_des.cmd'))
        for i in range(51):
            before=next(r for r in old if r['case']==cc['case'] and int(r['index'])==i)
            after=next(r for r in new if r['case']==cc['case'] and int(r['index'])==i)
            assert before['native_qualified']=='True' and after['native_qualified']=='True' and float(before['vg'])==float(after['vg'])
            p=float(before['Id_A_per_um']);n=float(after['Id_A_per_um'])
            rows.append(dict(case=cc['case'],device=cc['device'],vd=cc['vd'],vg=float(after['vg']),phumob_native_A_per_um=p,enormal_native_A_per_um=n,delta_A_per_um=n-p,native_model_change_percent=100*(n/p-1)))
    a.write_csv(O/'points.csv',rows);a.write_csv(O/'identity.csv',checks)
    summary=dict(points=204,scope='Native Enormal versus native PhuMob finite model effect. Both native stages qualified. This is not Vela current error, a small-perturbation derivative, or a production acceptance gate.',curves=[])
    fig,axes=plt.subplots(1,2,figsize=(10,4),sharex=True)
    for col,device in enumerate(('n19','n23')):
        for vd,color in ((.05,'#1768ac'),(1.,'#d04b35')):
            group=[r for r in rows if r['device']==device and r['vd']==vd]
            summary['curves'].append(dict(device=device,vd=vd,min_change_percent=min(r['native_model_change_percent'] for r in group),max_change_percent=max(r['native_model_change_percent'] for r in group),endpoints={str(r['vg']):r['native_model_change_percent'] for r in group if r['vg'] in (0,.2,.8,1.)}))
            axes[col].plot([r['vg'] for r in group],[r['native_model_change_percent'] for r in group],color=color,label=f'Vd={vd:g} V')
        axes[col].set_title(device);axes[col].set_xlabel('Vg (V)');axes[col].set_ylabel('Native Enormal / native PhuMob - 1 (%)');axes[col].grid(alpha=.25);axes[col].legend();axes[col].set_xlim(0,1)
    fig.suptitle('Physical effect of restoring Enormal (native solver only)');fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(O/('model_effect.'+ext),dpi=180)
    plt.close(fig);a.write(O/'summary.json',summary)
    d.matrix.freeze(O/'evidence.json',files+[O/name for name in ('points.csv','identity.csv','summary.json','model_effect.png','model_effect.pdf')]);print(summary,flush=True)

if __name__=='__main__':main()
