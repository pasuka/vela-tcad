"""Seal qualified Enormal curve results and scientific plots without promotion."""
from pathlib import Path
import tarfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import simplemos_enormal_curves_20260912 as e
a,d=e.a,e.d


def main():
    for p in ('supplemental_source_evidence.json','comparison_evidence.json','fields/evidence.json'):
        a.verify(e.O/p)
    summary=a.read(e.O/'summary.json')
    assert all(summary[k]==204 for k in ('points','native_qualified','continuation_qualified','independent_qualified','dual_qualified','comparison_qualified'))
    assert a.read(e.O/'fields/summary.json')['all_reconstruction_qualified']
    assert not (e.O/'completion_evidence.json').exists()
    points=a.rows(e.O/'comparison.csv')
    attempts=[r for arm in ('continuation','native') for r in a.rows(e.O/(arm+'_attempts.csv'))]
    final=dict(summary,attempts=len(attempts),failed_attempts=sum(r['qualified']!='True' for r in attempts),reloads=sum(int(r['attempt'])>0 for r in attempts),
        max_row_ratio=max(float(r['max_row_ratio']) for r in attempts if r['qualified']=='True'),
        max_port_relative=max(float(r['port_relative']) for r in attempts if r['qualified']=='True'),
        max_dual_Id_relative=max(float(r['dual_Id_relative']) for r in points),
        max_dual_psi_V=max(float(r['psi_max_V']) for r in points),
        max_dual_phin_V=max(float(r['phin_max_V']) for r in points),
        max_dual_phip_V=max(float(r['phip_max_V']) for r in points),
        max_dual_density_relative=max(float(r['density_max_relative']) for r in points),
        all_qualified=True,acceptance_changed=False,production_default_changed=False,
        scope='300 K PhuMob + Enormal, OldSlotboom, doping SRH; original n19/n23 meshes, Vd .05/1, Vg 0..1 step .02. Explicit split state and geometry bundle; no HFS.',
        limits='PhuMob internal native minimum-search/weak hole cell gate remains unresolved. Source terms below the existing source floor do not gain independent relative-source qualification. Current curve qualification does not certify every local physical field. HFS and original n17-n24 0..2.5 V matrix remain separate gates.',
        runner=str(e.q.RUNNER),runner_sha256=a.sha(e.q.RUNNER))
    a.write(e.O/'final_summary.json',final)
    fig,axes=plt.subplots(2,2,figsize=(11,7),sharex=True)
    for col,device in enumerate(('n19','n23')):
        for vd,color in ((.05,'#1768ac'),(1.,'#d04b35')):
            rs=sorted([r for r in points if r['device']==device and float(r['vd'])==vd],key=lambda r:float(r['vg']))
            x=[float(r['vg']) for r in rs];ref=np.array([float(r['sentaurus_Id_A_per_um']) for r in rs]);cur=np.array([float(r['continuation_Id_A_per_um']) for r in rs])
            axes[0,col].semilogy(x,abs(ref),color=color,label=f'Sentaurus Vd={vd:g} V')
            axes[0,col].semilogy(x,abs(cur),linestyle='none',marker='o',markersize=3,markevery=3,fillstyle='none',color=color,label=f'Vela Vd={vd:g} V')
            axes[1,col].plot(x,[float(r['continuation_error_percent']) for r in rs],color=color,label=f'Vd={vd:g} V')
        axes[0,col].set_title(device);axes[0,col].set_ylabel('|Id| (A/um)');axes[0,col].legend(fontsize=8)
        axes[1,col].set_xlabel('Vg (V)');axes[1,col].set_ylabel('100 (Id Vela / Id Sentaurus - 1) (%)')
        axes[1,col].axhline(0,color='gray',lw=.6)
        for ax in axes[:,col]:ax.grid(alpha=.25);ax.set_xlim(0,1)
    fig.suptitle('PhuMob + Enormal / 204 points with dual initialization');fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(e.O/('curves.'+ext),dpi=180)
    plt.close(fig)
    files=[p for folder in ('src','include','tests') for p in (e.R/folder).rglob('*') if p.is_file() and p.suffix in ('.cpp','.h','.py')]
    files+=list((e.R/'scripts').glob('*simplemos*.py'))
    files += [e.R/'CMakeLists.txt',e.R/'CMakePresets.json',e.R/'docs/config_schema.md',e.q.RUNNER]
    files=sorted(set(files));archive=e.L/'qualified_source_and_binary.tgz'
    with tarfile.open(archive,'w:gz') as tar:
        for p in files:tar.add(p,arcname=p.relative_to(e.R).as_posix())
    a.write(e.O/'qualified_source_archive.json',dict(archive=str(archive),sha256=a.sha(archive),files={a.rel(p):a.sha(p) for p in files},candidate_archive=a.read(e.q.O/'archive_evidence.json')))
    paths=[Path(__file__).resolve(),e.O/'final_summary.json',e.O/'comparison_evidence.json',e.O/'fields/evidence.json',e.O/'supplemental_source_evidence.json',e.O/'qualified_source_archive.json',archive,e.O/'curves.png',e.O/'curves.pdf']
    d.matrix.freeze(e.O/'completion_evidence.json',paths)
    print(final,flush=True)

if __name__=='__main__':main()
