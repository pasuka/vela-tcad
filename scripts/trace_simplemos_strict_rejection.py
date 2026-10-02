"""Diagnose three strict Newton rejections without changing acceptance or numerics."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import math
import run_simplemos_convergence_acceptance_isolation as previous
import run_simplemos_convergence_audit as audit

REPO,ROOT=previous.REPO,previous.ROOT
LOCAL=REPO/'build-release/simplemos_strict_rejection_trace'
OUT=ROOT/'strict_rejection_trace'
RUNNER=LOCAL/'trace_runner.exe'
CONTRACT=OUT/'contract.json'
SCRIPT=Path(__file__).resolve()

TRACE=r'''
    nlohmann::json trace = nlohmann::json::array();
    for (const auto& t : result.trace) {
        nlohmann::json trials = nlohmann::json::array();
        for (const auto& h : t.lineSearchHistory) {
            trials.push_back({{"attempt", h.attempt}, {"damping", h.damping},
                {"residual", h.residualNorm}, {"target", h.targetResidualNorm},
                {"finite", h.finite}, {"caller_valid", h.acceptedByCaller},
                {"accepted", h.accepted}, {"reason", h.rejectionReason}});
        }
        trace.push_back({{"iteration", t.iter}, {"residual", t.residualNorm},
            {"blocks", blockResidualsJson(t.blockResiduals)},
            {"event", t.event}, {"accepted", t.lineSearchAccepted},
            {"carrier_rows", carrierRowConvergenceJson(t.carrierRowConvergence)},
            {"trials", trials}});
    }
'''


def env():
    e=os.environ.copy();e['PATH']='D:/msys64/ucrt64/bin'+os.pathsep+e['PATH'];return e


def build():
    if RUNNER.exists():raise FileExistsError(RUNNER)
    LOCAL.mkdir(parents=True,exist_ok=True)
    text=(REPO/'src/tools/vela_example_runner.cpp').read_text()
    start=text.index('nlohmann::json runNewtonSolveFromState(')
    end=text.index('std::vector<vela::Real> readNodeScalarCsv',start)
    part=text[start:end]
    assert part.count('    return {')==1
    part=part.replace('    return {',TRACE+'\n    return {\n        {"newton_trace", trace},',1)
    (LOCAL/'runner.cpp').write_text(text[:start]+part+text[end:])
    args=audit.read(REPO/'build-release/simplemos_convergence_audit/build_command.json')
    args=[v for v in args if not v.endswith('/NewtonSolver.cpp') and not v.endswith('\\NewtonSolver.cpp') and not ('simplemos_convergence_audit' in v and v.startswith('-I'))]
    args=[str(LOCAL/'runner.cpp') if v.endswith('runner.cpp') else v for v in args]
    args[-1]=str(RUNNER)
    audit.write(LOCAL/'build_command.json',args)
    result=subprocess.run(args,env=env(),capture_output=True,text=True)
    (LOCAL/'build.log').write_text(result.stdout+result.stderr)
    if result.returncode:raise RuntimeError(result.stderr[-3000:])
    print('Built runner-only trace overlay; original core library retained',flush=True)


def prepare():
    if CONTRACT.exists():raise FileExistsError(CONTRACT)
    audit.verify(previous.FREEZE)
    cases=[r for r in audit.rows(previous.OUT/'case_ledger.csv') if r['qualified']=='False']
    inputs=[SCRIPT,RUNNER,REPO/'build-release/libvela_core.a',REPO/'src/solver/NewtonSolver.cpp',
            previous.FREEZE,previous.CONTRACT,previous.OUT/'case_ledger.csv']
    for c in cases:
        source=previous.LOCAL/c['case'];dest=LOCAL/c['case']
        deck=audit.read(source/'config.json')
        mesh=audit.read(Path(deck['mesh_file']));nodes=[int(n['id']) for n in mesh['nodes']]
        deck['output_state_file']=str(dest/'state.csv')
        deck['solver']['local_update_diagnostics']={'enabled':True,'nodes':nodes,'csv_file':str(dest/'updates.csv'),'first_iterations':200,'every_iterations':1}
        cr=deck['solver']['carrier_row_convergence']
        cr.update(diagnostic_csv=str(dest/'violations.csv'),trace_csv=str(dest/'rows.csv'),trace_nodes=nodes,trace_first_iterations=200,trace_every_iterations=1)
        audit.write(dest/'config.json',deck)
        inputs += [source/'config.json',source/'state.csv',source/'status.json',dest/'config.json',Path(deck['state_file']),Path(deck['mesh_file'])]
    audit.write(CONTRACT,{'status':'frozen_before_execution','cases':cases,'additional_nonlinear_reclosures':3,
        'change':'Output-only runner serialization and existing local-update/row diagnostic CSV switches. Same original states, same solver library, same local eps=1e-6 and external KCL=1e-8. No cap, Jacobian, damping, scaling, tolerance or recovery changes.',
        'identity_gate':'Final CSV SHA256, iteration count and failure reason must equal prior isolation run before interpreting traces.',
        'trace_warning':'Carrier trace CSV scale/ratio is a local flux/source display formula, not necessarily the authoritative eligibility/scale used by carrier_row_convergence; use status violations for final qualification.',
        'inputs':{audit.rel(p):audit.sha(p) for p in sorted(set(inputs))}})


def verify():
    for p,h in audit.read(CONTRACT)['inputs'].items():
        if audit.sha(REPO/p)!=h:raise ValueError(f'Input changed: {p}')


def run():
    verify();results=[];violations=[];trials=[]
    for c in audit.read(CONTRACT)['cases']:
        dest=LOCAL/c['case'];source=previous.LOCAL/c['case']
        if not (dest/'status.json').exists():
            proc=subprocess.run([str(RUNNER),'--config',str(dest/'config.json'),'--log','off'],env=env(),capture_output=True,text=True)
            (dest/'stdout.txt').write_text(proc.stdout);(dest/'stderr.txt').write_text(proc.stderr)
            status=json.loads(proc.stdout.strip().splitlines()[-1]);status['exit_code']=proc.returncode
            audit.write(dest/'status.json',status)
        status=audit.read(dest/'status.json')
        same=audit.sha(dest/'state.csv')==audit.sha(source/'state.csv')
        if not same or status['iterations']!=int(c['iterations']) or status['failure_reason']!=c['failure']:
            raise ValueError('Diagnostic-only identity failed')
        updates=audit.rows(dest/'updates.csv');last=max(int(r['iteration']) for r in updates)
        final=[r for r in updates if int(r['iteration'])==last]
        traces=status['newton_trace'];last_trace=traces[-1]
        for t in last_trace['trials']:trials.append({'case':c['case'],**t})
        mesh=audit.read(Path(audit.read(dest/'config.json')['mesh_file']))
        coords={int(n['id']):(n['x'],n['y']) for n in mesh['nodes']}
        caps=sum(float(r['raw_linear_step_V'])!=float(r['capped_step_V']) for r in final)
        result={'case':c['case'],'state_identity':same,'iterations':status['iterations'],'failed_iteration':last,
                'failure_reason':status['failure_reason'],'local_violations':status['carrier_row_convergence']['violation_count'],
                'capped_carrier_rows_on_failed_step':caps,'max_raw_step_V':max(abs(float(r['raw_linear_step_V'])) for r in final),
                'raw_linear_residual_L2':float(final[0]['raw_linear_residual_l2']),
                'capped_linear_residual_L2':float(final[0]['capped_linear_residual_l2']),
                'global_residual':last_trace['residual'],'trial_count':len(last_trace['trials']),
                'best_trial_over_current':min(t['residual'] for t in last_trace['trials'])/last_trace['residual']}
        for v in status['carrier_row_convergence']['violations']:
            row=next(r for r in final if int(r['node_id'])==v['node_id'] and r['carrier']==v['carrier'])
            x,y=coords[v['node_id']]
            violations.append({'case':c['case'],'node_id':v['node_id'],'carrier':v['carrier'],'x_um':x,'y_um':y,
                'ratio':v['ratio'],'residual_A_per_um':v['residual']*1e-8,'scale_A_per_um':v['scale']*1e-8,
                **{k:row[k] for k in ('qf_reference_V','state_increment_V','state_physical_qf_V','raw_linear_step_V','capped_step_V','applied_step_V','raw_linear_residual','selected_trial_residual','jacobian_diagonal','row_weight')},
                'increment_ulp_V':math.ulp(float(row['state_increment_V']))})
        results.append(result);print(json.dumps(result),flush=True)
    audit.write_csv(OUT/'case_ledger.csv',results);audit.write_csv(OUT/'violating_rows.csv',violations);audit.write_csv(OUT/'rejected_trials.csv',trials)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('build','prepare','run','verify'))
    {'build':build,'prepare':prepare,'run':run,'verify':verify}[p.parse_args().action]()
