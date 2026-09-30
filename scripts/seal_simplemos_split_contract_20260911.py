"""Seal final split-state qualification without rewriting failed evidence."""
from pathlib import Path
import ast,re
import validate_simplemos_split_contract_v2_20260911 as v
R=v.ROOT;L=v.LOCAL.parent;O=v.OUT
# Verify historical build inputs against preserved source preimages, not today's files.
substitutions={
 str(O.parent/'build_evidence.json'):{'scripts/validate_simplemos_split_contract_20260911.py':L/'validate_before_wide.py','include/vela/numerics/SplitDDState.h':L/'before_wide.h','include/vela/equation/SplitDDOperator.h':L/'operator_before_frame_guard.h'},
 str(O/'build_evidence.json'):{'include/vela/numerics/SplitDDState.h':L/'state_before_frame_guard.h','include/vela/equation/SplitDDOperator.h':L/'operator_before_frame_guard.h'},
 str(O/'frame_guard/build.json'):{'include/vela/equation/SplitDDOperator.h':L/'operator_before_port_sign.h'}}
proof=[]
for manifest,replacements in substitutions.items():
 for original,digest in v.read(manifest)['input_hashes'].items():
  resolved=replacements.get(original,R/original);assert v.prior.sha(resolved)==digest,(manifest,original,resolved)
  if original in replacements:proof.append(dict(manifest=str(Path(manifest).relative_to(R)),original=original,archived=str(resolved.relative_to(R)),sha256=digest))
checked=[O/'run_evidence.json',O/'analysis_inputs.json',O/'analysis_v3/analysis_evidence.json',O/'hotspot/evidence.json',O/'response_extension/evidence.json',O/'merit_correct/evidence.json',O/'frame_guard/replay.json',O/'final_kernel/build.json',O/'final_kernel/replay.json',O/'port_micro/evidence.json']
for p in checked:v.verify(p)
assert len(v.rows(O/'final_kernel/replay.csv'))==6
assert all(r['nonport_outputs_identical']=='True' for r in v.rows(O/'final_kernel/replay.csv'))
assert sum(r['descent']=='True' for r in v.rows(O/'merit_correct/merit.csv'))==12
assert 'All tests passed (48 assertions in 3 test cases)' in (L/'final_unit.log').read_text()
assert '100% tests passed out of 23' in (L/'final_ctest.log').read_text()
report=R/'docs/validation/simplemos_unified_split_state_contract_2026-09-11.md';status=R/'docs/validation/simplemos_branch_status.md'
for target in re.findall(r'\]\(([^)]+)\)',report.read_text()):
 if target.startswith('http'):continue
 p=(report.parent/target).resolve()
 assert p.exists() or p==O/'completion_evidence.json',target
scripts=[R/'scripts'/n for n in ('validate_simplemos_split_contract_20260911.py','validate_simplemos_split_contract_v2_20260911.py','analyze_simplemos_split_contract_20260911.py','analyze_simplemos_split_contract_v2_20260911.py','analyze_simplemos_split_contract_v3_20260911.py','replay_simplemos_split_contract_guard_20260911.py','replay_simplemos_split_contract_final_20260911.py','localize_simplemos_split_contract_baseline_20260911.py','check_simplemos_split_contract_ports_20260911.py','check_simplemos_split_contract_merit_20260911.py','check_simplemos_split_contract_port_micro_20260911.py')]
for p in scripts:ast.parse(p.read_text(encoding='utf-8-sig'))
v.write(O/'source_preimages.json',proof)
v.write(O/'qualification_summary.json',dict(scope='experimental unified state and fixed-state qualification; no production Newton migration or DC',split_partition_cases=6,main_checkpoint_jobs=54,additional_port_jobs=8,hotspot_micro_directions=18,nonzero_edge_responses=198,nonzero_port_directions=6,small_step_jv_blocks_passed=18,full_step_jv_blocks_passed=17,full_step_jv_blocks_total=18,whole_dd_merit_descent=12,old_production_residual_cases_passed=5,old_production_residual_cases_total=6,old_production_residual_failure='n23 Vd=1 Vg=0 hole node 1087; frozen 1e-10 gate retained',formal_dual_initialization_cases='10/16 unchanged',new_test_assertions=48,ctest='23/23',sentaurus_runs=0))
files=checked+[O.parent/'build_evidence.json',O/'build_evidence.json',O/'frame_guard/build.json',O/'source_preimages.json',O/'qualification_summary.json',O/'acceptance.json',report,status,R/'CMakeLists.txt',R/'tests/test_split_dd_state.cpp',R/'include/vela/numerics/SplitDDState.h',R/'include/vela/equation/SplitDDOperator.h',R/'build-release/test_split_dd_state.exe',L/'final_build.log',L/'final_unit.log',L/'final_ctest.log',L/'analysis.log',L/'analysis_v2.log',L/'analysis_v3.log']+scripts+[R/p['archived'] for p in proof]+[Path(__file__).resolve()]
v.freeze(O/'completion_evidence.json',list(dict.fromkeys(files)));v.verify(O/'completion_evidence.json');print('Final evidence verified:',len(set(files)),'files')
