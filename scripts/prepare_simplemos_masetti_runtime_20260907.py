"""Freeze native runtime data, separate from TDR plot semantics."""
import argparse
from pathlib import Path
import shutil
import tarfile
import prepare_simplemos_masetti_local_20260907 as prior
a=prior.a;d=prior.d;REPO=prior.REPO
LOCAL=REPO/'build-release/simplemos_masetti_runtime_20260907'
OUT=d.ROOT/'masetti_runtime_calibration_20260907'
REMOTE='/tmp/vela_simplemos_masetti_runtime_20260907'

TCL=r'''
proc tcl_cp_constructor {} {
 upvar #1 dump_iteration dump_iteration
 set dump_iteration 0
}
proc tcl_cp_Compute_Dataset_Names {} { return [list "VelaRuntimeAudit"] }
proc tcl_cp_Compute_Function_Names {} { return [list "VelaRuntimeAudit"] }
proc tcl_cp_Compute_Plot_Values {} {
 upvar #1 tcl_cp_adr tcl_cp_adr
 upvar #1 dump_iteration dump_iteration
 set mesh [$tcl_cp_adr Mesh];set data [$tcl_cp_adr Data]
 set prefix [format "runtime_%03d" $dump_iteration];incr dump_iteration
 set meta [open "${prefix}_meta.txt" w]
 puts $meta "vertices=[$mesh size_vertex] edges=[$mesh size_edge] elements=[$mesh size_element] element_vertices=[$mesh size_element_vertex]"
 foreach {tag location size} [list vertex $::des_data_vertex [$mesh size_vertex] edge $::des_data_edge [$mesh size_edge] element $::des_data_element [$mesh size_element] element_vertex $::des_data_element_vertex [$mesh size_element_vertex]] {
  foreach name {eMobility hMobility} {
   set ptr [$data ReadScalar $location $name]
   set f [open "${prefix}_${tag}_${name}.csv" w];puts $f "index,value"
   for {set i 0} {$i<$size} {incr i} {puts $f "$i,[format %.17g [tcl_cp_get_double $ptr $i]]"}
   close $f
  }
 }
 foreach name {eDensity hDensity Potential eQuasiFermi hQuasiFermi DonorConcentration AcceptorConcentration SRHRecombination} {
  set ptr [$data ReadScalar $::des_data_vertex $name]
  set f [open "${prefix}_vertex_${name}.csv" w];puts $f "index,value"
  for {set i 0} {$i<[$mesh size_vertex]} {incr i} {puts $f "$i,[format %.17g [tcl_cp_get_double $ptr $i]]"}
  close $f
 }
 foreach {tag location size} [list edge $::des_data_edge [$mesh size_edge] element $::des_data_element [$mesh size_element]] {
  foreach name {eCurrent hCurrent} {
   set ptr [$data ReadVector $location $name]
   set f [open "${prefix}_${tag}_${name}.csv" w];puts $f "index,x,y"
   for {set i 0} {$i<$size} {incr i} {puts $f "$i,[format %.17g [tcl_cp_get_double2 $ptr 0 $i]],[format %.17g [tcl_cp_get_double2 $ptr 1 $i]]"}
   close $f
  }
 }
 set f [open "${prefix}_vertices.csv" w];puts $f "index,x_um,y_um"
 for {set i 0} {$i<[$mesh size_vertex]} {incr i} {
  set v [$mesh vertex $i];puts $f "[$v index],[format %.17g [$v coord 0]],[format %.17g [$v coord 1]]"
 };close $f
 set f [open "${prefix}_edges.csv" w];puts $f "index,start,end"
 for {set i 0} {$i<[$mesh size_edge]} {incr i} {
  set e [$mesh edge $i];set v0 [$e start];set v1 [$e end];puts $f "[$e index],[$v0 index],[$v1 index]"
 };close $f
 set measures [$data ReadMeasure];set coeffs [$data ReadCoefficient]
 set f [open "${prefix}_element_vertices.csv" w];puts $f "element,region,material,local,vertex,element_vertex,measure_um2"
 set ef [open "${prefix}_element_edges.csv" w];puts $ef "element,local,edge,coefficient"
 for {set i 0} {$i<[$mesh size_element]} {incr i} {
  set el [$mesh element $i];set ei [$el index];set bulk [$el bulk]
  for {set j 0} {$j<[$el size_vertex]} {incr j} {
   set v [$el vertex $j]
   puts $f "$ei,[$bulk name],[$bulk material],$j,[$v index],[$v element_vertex_index $el],[format %.17g [tcl_cp_get_double2 $measures $ei $j]]"
  }
  for {set j 0} {$j<[$el size_edge]} {incr j} {
   set e [$el edge $j];puts $ef "$ei,$j,[$e index],[format %.17g [tcl_cp_get_double2 $coeffs $ei $j]]"
  }
 };close $f;close $ef;close $meta
 return [list 0.0]
}
'''

def prepare():
 a.verify(prior.OUT/'validation_evidence.json');files=[Path(__file__).resolve(),prior.OUT/'validation_evidence.json',prior.MANUAL];jobs=[]
 for c in a.read(prior.OUT/'vela_contract.json')['cases']:
  if c['vg']!=1.:continue
  source=prior.prev.LOCAL/'native_raw/bundle/masetti'/c['case'];root=LOCAL/'bundle'/c['key'];root.mkdir(parents=True,exist_ok=False)
  for name in ('input_fps.tdr','result_020_des.sav','result_020_circuit_des.sav'):
   shutil.copyfile(source/name,root/name);files += [source/name,root/name]
  text=(source/'native_des.cmd').read_text().split('Plot {',1)[0]
  text+='''Plot { eDensity hDensity Potential eQuasiFermi hQuasiFermi eMobility/Element hMobility/Element eCurrent/Vector/Element hCurrent/Vector/Element SRHRecombination }
CurrentPlot { Tcl(tcl="source runtime.tcl") }
'''+prior.prev.p.MATH+'''
Solve { Load(FilePrefix="result_020") Coupled { Poisson Electron Hole } Plot(FilePrefix="runtime_export") }
'''
  (root/'native_des.cmd').write_text(text,newline='\n');(root/'runtime.tcl').write_text(TCL,newline='\n')
  files += [root/'native_des.cmd',root/'runtime.tcl'];jobs.append(dict(key=c['key'],device=c['device'],vd=c['vd'],vg=1.,native_Id_A_per_um=c['native_Id_A_per_um']))
 shell='''#!/bin/bash
set -u
cd "$(dirname "$0")"
for path in bundle/*; do
 (cd "$path"; sdevice native_des.cmd > console.log 2>&1; printf '%s\\n' "$?" > exit_code.txt)
done
tar -czf results.tgz bundle
printf 'complete\\n' > complete.txt
'''
 (LOCAL/'run.sh').write_text(shell,newline='\n');files.append(LOCAL/'run.sh')
 a.write(OUT/'native_contract.json',dict(status='frozen_before_execution',jobs=jobs,remote=REMOTE,
  scope='Read-only CurrentPlot Tcl runtime scalar/vector data and explicit internal topology; four strict Masetti endpoints at Vg=1.',model_and_acceptance_changes=False,
  forbidden_interpretation='ReadFlux is a scalar Laplacian-like quantity, not carrier edge current; it is not called.',
  gates=dict(current_reclosure_relative=1e-8,kcl_over_Id=1e-8,coordinate_um=1e-12,mobility_identity_relative=1e-8,edge_terminal_relative=1e-6),
  semantics='No runtime mobility or current location is declared the solver edge operator before numerical conservation and same-state checks.',
  eventual_response_gates=dict(prediction_relative=.001,two_amplitude_relative=.001,even_over_odd=.01,zero_drift_dex=1e-5,minimum_signal_over_zero_drift=100),
  source_calibration='A distributed source is frozen only after its physical definition and support are explicit; no production model replacement.'))
 files.append(OUT/'native_contract.json');d.matrix.freeze(OUT/'native_freeze.json',files)
 with tarfile.open(LOCAL/'input.tgz','w:gz') as t:t.add(LOCAL/'bundle',arcname='bundle');t.add(LOCAL/'run.sh',arcname='run.sh')
 print('Frozen four runtime exports',flush=True)

def unpack():
 a.verify(OUT/'native_freeze.json');dest=LOCAL/'native_raw';assert not dest.exists()
 with tarfile.open(LOCAL/'results.tgz') as t:
  for m in t.getmembers():assert (dest/m.name).resolve().is_relative_to(dest.resolve()) and not m.issym() and not m.islnk()
  t.extractall(dest,filter='data')
 for source in (LOCAL/'bundle').rglob('*'):
  if source.is_file():assert a.sha(source)==a.sha(dest/'bundle'/source.relative_to(LOCAL/'bundle'))
 print('Returned input hashes verified',flush=True)

if __name__=='__main__':
 q=argparse.ArgumentParser();q.add_argument('action',choices=('prepare','unpack'));globals()[q.parse_args().action]()
