"""Sentaurus PMI same-integrated-source response; no lifetime replacement."""
import argparse
import math
import shutil
import tarfile
from pathlib import Path
import calibrate_simplemos_phumob_local_source_20260909 as s

a,d=s.a,s.d
LOCAL=s.LOCAL/'native'
OUT=s.OUT/'native'
REMOTE='/tmp/vela_simplemos_phumob_pair_source_20260909'
PMI=r'''
#include "PMI.h"
#include <fstream>
#include <vector>
#include <set>
#include <cmath>
#include <iomanip>
#include <stdexcept>
class PairSource : public PMI_Recombination_Base {
  struct Site {double x,y,r;};
  std::vector<Site> sites; std::set<size_t> logged; double alpha;
public:
  PairSource(const PMI_Environment& env):PMI_Recombination_Base(env) {
    std::ifstream f("pair_source.txt");size_t n;
    if(!(f>>alpha>>n)||n!=8||!std::isfinite(alpha)||std::abs(alpha)>.001)throw std::runtime_error("Invalid pair source header");
    for(size_t k=0;k<n;++k){Site p;if(!(f>>p.x>>p.y>>p.r)||!std::isfinite(p.r))throw std::runtime_error("Invalid pair source record");sites.push_back(p);}
    std::string extra;if(f>>extra)throw std::runtime_error("Extra pair source data");
  }
  void compute(const Input& in,Output& out) {
    double x,y,z;in.ReadCoordinate(x,y,z);out.r=0.;
    for(size_t k=0;k<sites.size();++k)if(std::abs(x-sites[k].x)<1e-10&&std::abs(y-sites[k].y)<1e-10){
      out.r=pmi_float(alpha)*pmi_float(sites[k].r);
      if(logged.insert(k).second){std::ofstream f("pair_source_seen.csv",std::ios::app);f<<std::setprecision(17)<<k<<','<<x<<','<<y<<','<<sites[k].r<<','<<alpha<<','<<pmi_float::get_precision()<<'\n';}
      return;
    }
  }
};
extern "C" PMI_Recombination_Base* new_PMI_Recombination_Base(const PMI_Environment& env){return new PairSource(env);}
extern "C" int PMI_Recombination_ElectricField(){return 0;}
'''

def prepare():
    a.verify(s.OUT/'freeze.json');scope=a.rows(s.OUT/'scope.csv');files=[Path(__file__).resolve(),s.OUT/'freeze.json'];jobs=[]
    LOCAL.mkdir(parents=True,exist_ok=False)
    model=LOCAL/'pairsource.C';model.write_text(PMI,newline='\n');files.append(model)
    for c in a.read(s.OUT/'contract.json')['cases']:
        src=s.prior.NLOCAL/'native_raw/bundle/phumob'/(c['key']+'_vg_000')
        for label,alpha in s.AMPLITUDES:
            dest=LOCAL/'bundle'/c['key']/label;dest.mkdir(parents=True,exist_ok=False)
            for old,new in (('input_fps.tdr','input_fps.tdr'),('final_des.sav','seed_des.sav'),('final_circuit_des.sav','seed_circuit_des.sav')):
                shutil.copyfile(src/old,dest/new);files += [src/old,dest/new]
            deck=(src/'native_des.cmd').read_text().replace('File {','File { PMIPath="../../.."',1).replace('Recombination(SRH(DopingDependence))','Recombination(SRH(DopingDependence) pairsource)').replace('Load(FilePrefix="result_000")','Load(FilePrefix="seed")')
            deck=deck.replace('CurrentPlot { Tcl(tcl="source runtime.tcl") }','')
            (dest/'native_des.cmd').write_text(deck,newline='\n')
            rows=[r for r in scope if r['key']==c['key']];assert len(rows)==8
            (dest/'pair_source.txt').write_text(f'{alpha:.17g} 8\n'+''.join(f"{float(r['x_um']):.17g} {float(r['y_um']):.17g} {float(r['native_rate_cm3_s']):.17g}\n" for r in rows),newline='\n')
            files += [src/'native_des.cmd',dest/'native_des.cmd',dest/'pair_source.txt']
            jobs.append(dict(key=c['key'],device=c['device'],vd=c['vd'],vg=0.,label=label,alpha=alpha))
    script='''#!/bin/bash
set -eu
cd "$(dirname "$0")"
/atctools/Synopsys/tcad/T-2022.03/bin/cmi -O pairsource.C > compile.log 2>&1
for path in bundle/*/zero; do
 (cd "$path"; set +e; /atctools/Synopsys/tcad/T-2022.03/bin/sdevice native_des.cmd > console.log 2>&1; printf '%s\\n' "$?" > exit_code.txt) &
done
wait
tar -czf zero_results.tgz bundle compile.log pairsource.so.*
printf 'zero_complete\\n' > zero_complete.txt
'''
    (LOCAL/'run_zero.sh').write_text(script,newline='\n');files.append(LOCAL/'run_zero.sh')
    script='''#!/bin/bash
set -eu
cd "$(dirname "$0")"
for label in plus_full minus_full plus_half minus_half; do
 for path in bundle/*/"$label"; do
  (cd "$path"; set +e; /atctools/Synopsys/tcad/T-2022.03/bin/sdevice native_des.cmd > console.log 2>&1; printf '%s\\n' "$?" > exit_code.txt) &
 done
 wait
done
tar -czf results.tgz bundle compile.log pairsource.so.*
printf 'complete\\n' > complete.txt
'''
    (LOCAL/'run_response.sh').write_text(script,newline='\n');files.append(LOCAL/'run_response.sh')
    a.write(OUT/'contract.json',dict(jobs=jobs,remote=REMOTE,source='same integrated particle source as local experiment; native cm^-3 s^-1 rate = particle source/(native Si box m2 * 1e6). Added PMI, retains native SRH and PhuMob.',manual='Local T-2022.03 UG pp1214-1216,1230,1240,1266-1268; simplified PMI output is constant w.r.t. all independent variables, derivatives exactly zero.',gates=dict(kcl=1e-8,bias_V=1e-10,zero_Id_relative=1e-8,matched_sites=8,prediction=1e-3,two_amplitude=1e-3,even_over_odd=.01,signal_over_drift=100),finite_volume_replacement=False))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json',s.REPO/'build-release/m79_research/sdevice_ug_local_2022.txt'])
    with tarfile.open(LOCAL/'input.tgz','w:gz') as t:
        for name in ('bundle','pairsource.C','run_zero.sh','run_response.sh'):t.add(LOCAL/name,arcname=name)
    print('Prepared four zero controls then 16 signed responses',flush=True)

def unpack(label):
    a.verify(OUT/'freeze.json');dest=LOCAL/(label+'_raw');dest.mkdir(exist_ok=False)
    with tarfile.open(LOCAL/(label+'_results.tgz' if label=='zero' else 'results.tgz')) as tar:
        for item in tar.getmembers():assert (dest/item.name).resolve().is_relative_to(dest.resolve()) and not item.issym() and not item.islnk()
        tar.extractall(dest,filter='data')
    for src in (LOCAL/'bundle').rglob('*'):
        if src.is_file():assert a.sha(src)==a.sha(dest/'bundle'/src.relative_to(LOCAL/'bundle'))
    rows=[];refs=a.rows(s.prior.NOUT/'native_points.csv')
    for job in a.read(OUT/'contract.json')['jobs']:
        if label=='zero' and job['label']!='zero':continue
        root=dest/'bundle'/job['key']/job['label'];code=int((root/'exit_code.txt').read_text());log=(root/'console.log').read_text(errors='replace')
        points=s.prior.native.original.b.exporter.pltrows(root/'native_des.plt');assert len(points)==1
        p=points[0];currents={n:p[n+' TotalCurrent'] for n in ('drain','source','gate','substrate')};kcl=abs(math.fsum(currents.values()))/abs(currents['drain']);bias=max(abs(p['gate OuterVoltage']),abs(p['drain OuterVoltage']-job['vd']))
        seen=[x.split(',') for x in (root/'pair_source_seen.csv').read_text().splitlines()];indices={int(x[0]) for x in seen};assert indices==set(range(8)) and all(int(x[5])==2 for x in seen)
        ref=next(r for r in refs if r['model']=='phumob' and r['case']==job['key'] and int(r['index'])==0)
        zeroerr=abs(currents['drain']/float(ref['Id_A_per_um'])-1) if job['label']=='zero' else 0.
        row=dict(**job,exit_code=code,Id_A_per_um=currents['drain'],kcl_over_Id=kcl,bias_error_V=bias,zero_Id_relative=zeroerr,matched_sites=len(indices),extended_precision_128=True,qualified=code==0 and 'Good Bye' in log and 'T-2022.03-SP2' in log and kcl<=1e-8 and bias<=1e-10 and zeroerr<=1e-8)
        row.update({n+'_A_per_um':v for n,v in currents.items()});rows.append(row)
    a.write_csv(OUT/(label+'_points.csv'),rows);d.matrix.freeze(OUT/(label+'_evidence.json'),[OUT/'freeze.json',OUT/(label+'_points.csv')]+[p for p in dest.rglob('*') if p.is_file()]);print(label,sum(r['qualified'] for r in rows),'/',len(rows),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','zero','response'));x=p.parse_args().action;prepare() if x=='prepare' else unpack(x)
