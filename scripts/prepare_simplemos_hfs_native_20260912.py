"""Prepare the next native HFS baseline/observer cohort only after Enormal curves."""
import shutil,tarfile
from pathlib import Path
import simplemos_enormal_curves_20260912 as e
import build_simplemos_hfs_observer_native_20260912 as build
a,d=e.a,e.d
L=e.R/'build-release/hfs_restore_20260912';O=e.R/'reference_tcad/simplemos_sentaurus2022/hfs_restore_20260912'
REMOTE='/tmp/vela_simplemos_hfs_restore_20260912'


def main():
    a.verify(e.O/'completion_evidence.json');assert a.read(e.O/'final_summary.json')['all_qualified']
    a.verify(build.O/'evidence.json');assert a.read(build.O/'summary.json')['compiled']
    assert not (O/'input_evidence.json').exists()
    jobs=[];files=[Path(__file__).resolve(),e.O/'completion_evidence.json',build.O/'evidence.json']
    pmi=L/'pmi';pmi.mkdir(parents=True,exist_ok=False)
    for name in ('canali_diagnostic.h','pmi_vela_hfs_observer.C','pmi_vela_hfs_observer.so.linux64'):
        src=build.L/'raw'/name;shutil.copyfile(src,pmi/name);files.extend([src,pmi/name])
    for c in a.read(e.q.O/'inputs.json')['cases']:
        src=e.L/'native_raw/bundle'/c['case']
        header=(src/'native_des.cmd').read_text().split('Solve {',1)[0]
        assert 'Mobility(PhuMob Enormal)' in header and 'HighFieldSaturation' not in header
        assert 'ComputeGradQuasiFermiAtContacts' not in header
        for arm in ('baseline','observed'):
            name=c['key']+'_'+arm;dest=L/'bundle'/name;dest.mkdir(parents=True,exist_ok=False)
            model='HighFieldSaturation' if arm=='baseline' else 'HighFieldSaturation(pmi_vela_hfs_observer)'
            text=header.replace('Mobility(PhuMob Enormal)',f'Mobility(PhuMob {model} Enormal)')
            if arm=='observed':text=text.replace('File {','File { PMIPath="../../pmi"',1)
            text+='''Solve {
 Load(FilePrefix="anchor")
 Coupled { Poisson Electron Hole }
 Plot(FilePrefix="final")
 Save(FilePrefix="final")
}
'''
            (dest/'native_des.cmd').write_text(text,newline='\n')
            for old,new in [('input_fps.tdr','input_fps.tdr'),(f"vg_{c['index']:03d}_des.sav",'anchor_des.sav'),(f"vg_{c['index']:03d}_circuit_des.sav",'anchor_circuit_des.sav')]:
                shutil.copyfile(src/old,dest/new);files.extend([src/old,dest/new])
            jobs.append(dict(case=c['case'],key=c['key'],device=c['device'],vd=c['vd'],vg=c['vg'],index=c['index'],arm=arm,name=name))
            files.append(dest/'native_des.cmd')
    pilot=jobs[:2]
    for label,group in [('pilot',pilot),('rest',jobs[2:])]:
        shell='''#!/bin/bash
set -u
cd "$(dirname "$0")" || exit 90
printf '%s\\n' "$$" > '''+label+'''_launcher.pid
running=0
for name in '''+' '.join(j['name'] for j in group)+'''; do
 (cd "bundle/$name" || exit 91
  test ! -f exit_code.txt || exit 92
  /atctools/Synopsys/tcad/T-2022.03/bin/sdevice native_des.cmd > console.log 2>&1 &
  pid=$!;printf '%s\\n' "$pid" > solver.pid
  wait "$pid";code=$?;printf '%s\\n' "$code" > exit_code.txt
 ) &
 running=$((running+1));if [ "$running" -eq 2 ]; then wait;running=0;fi
done
wait
tar czf '''+label+'''_results.tgz pmi '''+' '.join('bundle/'+j['name'] for j in group)+'''
printf 'complete; inspect each exit and identity gate\\n' > '''+label+'''_complete.txt
'''
        (L/('run_'+label+'.sh')).write_text(shell,newline='\n');files.append(L/('run_'+label+'.sh'))
    a.write(O/'contract.json',dict(remote_root=REMOTE,jobs=jobs,pilot=[j['name'] for j in pilot],
        scope='Restore HFS only after qualified Enormal curves. Eight targets, native baseline and mathematical Canali replacement observer. Original default contact field switch; no contact-QF override.',
        gates=dict(exit_zero=True,bias_V=1e-10,kcl_over_Id=1e-8,observer_Id_relative=1e-12,observer_psi_V=1e-12,observer_density_relative=1e-10,observer_mobility_relative=1e-12),
        protocol='Run the first pair, export and qualify actual baseline/observer identity before launching the other fourteen jobs. Compilation alone does not qualify samples. Missing or ambiguous cell context stays unqualified.',
        no_parameter_fit=True,production_changed=False,contact_QF_control='Separate diagnostic only after default baseline/observer identity passes. Do not substitute its physics for original deck.'))
    d.matrix.freeze(O/'input_evidence.json',files+[O/'contract.json'])
    with tarfile.open(L/'input.tgz','w:gz') as tar:
        for name in ('bundle','pmi','run_pilot.sh','run_rest.sh'):tar.add(L/name,arcname=name)
    print(L/'input.tgz',a.sha(L/'input.tgz'),flush=True)

if __name__=='__main__':main()
