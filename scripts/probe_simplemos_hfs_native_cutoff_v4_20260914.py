"""Native fixed-state threshold probe, explicitly NOT a DC qualification.

Modify copies of a checkpoint coherently in phi_n/phi_p/n/p, then Load/Plot
without Newton. Disable boundary projection and the contact E substitution
only in this diagnostic. Exported loaded values must match the requested
states; otherwise the experiment is invalid. Original checkpoints unchanged.
"""
import argparse
import math
import shutil
import tarfile
from pathlib import Path
import numpy as np
import h5py
import check_simplemos_hfs_native_identity_20260912 as n

a,d,c=n.a,n.d,n.calibration
L=n.L/'cutoff_probe_20260914/v4';O=n.O/'cutoff_probe_20260914/v4'
REMOTE='/tmp/vela_simplemos_hfs_cutoff_v4_20260914'
FORCES=(.5,.999,.999999,1.,1.000001,1.001,1.5)


def prepare():
    assert not (O/'input_evidence.json').exists()
    a.verify(n.O/'pilot_identity_evidence.json')
    assert a.read(n.O/'pilot_identity_summary.json')['all_identity_qualified']
    original=n.L/'pilot_raw/bundle/m65_n19_vd_0p050000_endpoint_vg_040_baseline'
    exported=n.L/'pilot_exports/m65_n19_vd_0p050000_endpoint_vg_040_baseline'
    xy={int(r['id']):float(r['x_um']) for r in a.rows(exported/'nodes.csv')}
    fields={key:c.scalar(exported,name+'_region0.csv') for key,name in
            [('psi','ElectrostaticPotential'),('phin','eQuasiFermiPotential'),('phip','hQuasiFermiPotential'),('n','eDensity'),('p','hDensity')]}
    contacts=set(c.g.geometry('n19')[0].contact_nodes)
    native_vt=1.380662e-23*300/1.602192e-19
    ids=sorted(fields['n']);(L/'states').mkdir(parents=True,exist_ok=False)
    requested=[]
    for index,F in enumerate(FORCES):
        dest=L/'states'/f's{index:02d}_des.sav';shutil.copyfile(original/'final_des.sav',dest)
        shutil.copyfile(original/'final_circuit_des.sav',L/'states'/f's{index:02d}_circuit_des.sav')
        with h5py.File(dest,'r+') as h:
            group=h['collection/geometry_0/state_0']
            named={z.attrs['name'].decode():z for z in group.values() if 'name' in z.attrs and int(z.attrs.get('region',-1))==0}
            # Assert the region-local ordering against the independently exported native state.
            for key,name in [('psi','ElectrostaticPotential'),('phin','eQuasiFermiPotential'),('phip','hQuasiFermiPotential')]:
                assert np.array_equal(named[name]['values'][:],np.array([fields[key][k] for k in ids]))
            phi=np.array([fields['phin'][k] if k in contacts else .02+F*xy[k]*1e-4 for k in ids])
            phi_p=np.array([fields['phip'][k] if k in contacts else .02+F*xy[k]*1e-4 for k in ids])
            nn=np.array([fields['n'][k]*math.exp((fields['phin'][k]-v)/native_vt) for k,v in zip(ids,phi)])
            pp=np.array([fields['p'][k]*math.exp((v-fields['phip'][k])/native_vt) for k,v in zip(ids,phi_p)])
            # Native checkpoints store quasi-Fermi potentials, not density state
            # datasets. Density is reconstructed by Load; verify it after export.
            for name,value in [('eQuasiFermiPotential',phi),('hQuasiFermiPotential',phi_p)]:
                named[name]['values'][:]=value
            requested += [dict(index=index,F_requested=F,node=k,phin=v,phip=vp,n=ne,p=hp,psi=fields['psi'][k]) for k,v,vp,ne,hp in zip(ids,phi,phi_p,nn,pp)]
    (L/'pmi/gcc_compat').mkdir(parents=True)
    shutil.copyfile(n.L/'pmi/canali_diagnostic.h',L/'pmi/canali_diagnostic.h')
    shutil.copyfile(c.L/'native_raw/pmi/gcc_compat/g++',L/'pmi/gcc_compat/g++')
    code=(n.L/'pmi/pmi_vela_hfs_observer.C').read_text()
    code=code.replace('#include <map>','#include <map>\n#include <vector>')
    code=code.replace('std::map<std::tuple<double,double,double>,std::array<double,14> > samples_;','std::vector<std::array<double,14> > samples_;')
    code=code.replace('row.second','row')
    code=code.replace('samples_[std::make_tuple(x,y,z)]={{','samples_.push_back({{')
    code=code.replace('static_cast<double>(++sequence_),distance}};','static_cast<double>(++sequence_),distance}});')
    code=code.replace('std::ofstream out(filename_.c_str());','std::ofstream out(filename_.c_str(),std::ios::app);')
    code=code.replace('out<<\"x_um','if(out.tellp()==0) out<<\"x_um')
    code=code.replace('for(const auto& row:samples_)','out<<std::setprecision(17);\n        for(const auto& row:samples_)')
    (L/'pmi/pmi_vela_hfs_observer.C').write_text(code,newline='\n')
    for arm in ('baseline','observed'):
        dest=L/'bundle'/arm;dest.mkdir(parents=True)
        shutil.copyfile(original/'input_fps.tdr',dest/'input_fps.tdr')
        text=(original/'native_des.cmd').read_text().split('Solve {')[0]
        text=text.replace('Math {','Math {\n ComputeGradQuasiFermiAtContacts=UseQuasiFermi\n -ParallelToInterfaceInBoundaryLayer\n',1)
        if arm=='observed':
            text=text.replace('File {','File { PMIPath="../../pmi"',1)
            text=text.replace('HighFieldSaturation','HighFieldSaturation(pmi_vela_hfs_observer)',1)
        text+='Solve {\n'+''.join(f' Load(FilePrefix="../../states/s{i:02d}")\n Plot(FilePrefix="s{i:02d}")\n' for i in range(len(FORCES)))+'}\n'
        (dest/'native_des.cmd').write_text(text,newline='\n')
    run='''#!/bin/bash
set -u
cd "$(dirname "$0")" || exit 90
(cd pmi || exit 91
 chmod +x gcc_compat/g++
 export PATH="$PWD/gcc_compat:$PATH"
 /atctools/Synopsys/tcad/T-2022.03/bin/cmi -v pmi_vela_hfs_observer.C > compile.log 2>&1
 printf '%s\\n' "$?" > compile.exit_code.txt
)
if [ "$(cat pmi/compile.exit_code.txt)" != 0 ]; then tar czf results.tgz pmi; exit 92; fi
for arm in baseline observed; do
 (cd "bundle/$arm" || exit 91
  /atctools/Synopsys/tcad/T-2022.03/bin/sdevice native_des.cmd > console.log 2>&1
  printf '%s\\n' "$?" > exit_code.txt
 ) &
done
wait
tar czf results.tgz bundle states pmi
printf 'complete; fixed-state probe only\\n' > complete.txt
'''
    (L/'run.sh').write_text(run,newline='\n')
    a.write_csv(O/'requested.csv',requested)
    a.write(O/'contract.json',dict(remote_root=REMOTE,forces_V_cm=FORCES,
        scope='Load/Plot on modified copies, no Coupled solve. Not a physical or DC accepted state. Both model arms must export the requested loaded state.',
        gates=dict(psi_V=1e-12,phif_V=1e-10,density_relative=1e-8,observer_identity_relative=1e-12),
        diagnostic_math=['ComputeGradQuasiFermiAtContacts=UseQuasiFermi','-ParallelToInterfaceInBoundaryLayer'],
        original_state_modified=False,native_thermal_voltage_V=native_vt,
        prior_failure='v3 requested changed contact values and predicted densities with modern SI Vt; both assumptions invalid. v4 preserves contact values and uses previously documented native constants. No gate relaxation.',
        sample_scope='Only nodes with no contact-adjacent cell enter uniform-field checks; all other calls retained separately.'))
    d.matrix.freeze(O/'input_evidence.json',[Path(__file__).resolve(),n.O/'pilot_identity_evidence.json',
        original/'final_des.sav',original/'final_circuit_des.sav',O/'requested.csv',O/'contract.json']+[p for p in L.rglob('*') if p.is_file()])
    with tarfile.open(L/'input.tgz','w:gz') as t:
        for name in ('bundle','states','pmi','run.sh'):t.add(L/name,arcname=name)
    print(L/'input.tgz',a.sha(L/'input.tgz'),flush=True)


def collect(sha):
    a.verify(O/'input_evidence.json');assert sha==a.sha(L/'results.tgz')
    raw=L/'raw';assert not raw.exists()
    with tarfile.open(L/'results.tgz') as t:
        for m in t.getmembers():assert (raw/m.name).resolve().is_relative_to(raw.resolve()) and (m.isfile() or m.isdir())
        t.extractall(raw,filter='data')
    for folder in ('bundle','states','pmi'):
        for p in (L/folder).rglob('*'):
            if p.is_file():assert a.sha(p)==a.sha(raw/p.relative_to(L))
    codes={arm:int((raw/'bundle'/arm/'exit_code.txt').read_text()) for arm in ('baseline','observed')}
    a.write(O/'native_summary.json',dict(exit_codes=codes,all_exit_zero=all(v==0 for v in codes.values()),DC_qualified=False))
    d.matrix.freeze(O/'native_evidence.json',[Path(__file__).resolve(),O/'input_evidence.json',O/'native_summary.json',L/'results.tgz']+[p for p in raw.rglob('*') if p.is_file()])
    print(codes,flush=True);assert all(v==0 for v in codes.values())
    for arm in codes:
        for i in range(len(FORCES)):
            n.p.e.c.n.exporter.export_one(dict(case=arm,index=i,tdr=str(raw/'bundle'/arm/f's{i:02d}_des.tdr'),export=str(L/'exports'/arm/f's{i:02d}')))
    d.matrix.freeze(O/'export_evidence.json',[O/'native_evidence.json']+[p for p in (L/'exports').rglob('*') if p.is_file()])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','collect'));p.add_argument('--sha');args=p.parse_args()
    if args.action=='prepare':prepare()
    else:collect(args.sha)
