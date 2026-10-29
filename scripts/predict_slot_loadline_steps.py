"""Read-only qualified-evidence replay; never invokes a solver or resumes a run."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path

import run_slot_ldmos_bv_campaign as campaign
import slot_loadline_step_control as steps


def replay(checkpoint, executed_driver, output):
    data = campaign.load(checkpoint)
    files = dict(data['identity']['files'])
    driver = str(Path(campaign.__file__).resolve())
    expected = files.pop(driver)
    if campaign.ll.sha256(executed_driver) != expected:
        raise campaign.ll.EvidenceError('archived executed driver hash mismatch')
    campaign.ll.verify_hashes(files)
    runner = next(Path(p) for p in files if Path(p).name == 'vela_example_runner.exe')
    records = data['records']
    accepted = [a for a in data['attempts'] if a['status']=='accepted' and 'candidate' in a]
    # Latest completed attempt and any earlier recovered attempts, not the
    # environment-interrupted proposal at the tail.
    chosen = [a for a in accepted if a['candidate'].get('recovered')][-1:] + accepted[-1:]
    checked = set()
    cases = []
    for attempt in chosen:
        folder = Path(attempt['folder'])
        k = next(i for i,r in enumerate(records) if Path(r['point']['iv_file']).is_relative_to(folder))
        if k < 3:
            continue
        history = records[:k]
        actual = steps.samples(records[:k+1])[-1]
        for i in set([p['record_index'] for p in steps.samples(history)] + [k]):
            if i in checked:
                continue
            r = records[i]
            backend = campaign.ll.SolverBackend(r['base'],Path.cwd(),runner)
            backend.validate(r['point'])
            backend.validate_confirmation(r['point'],r['outer_V'],r['confirmation'])
            checked.add(i)
        campaign.ll.verify_hashes(attempt['candidate']['evidence'])
        scalar = campaign.load(folder/'scalar/checkpoint.json')
        cost = steps.measure_cost(folder,scalar['elapsed_s'])
        options = []
        for cap in (.1,.2,.5):
            policy = steps.Policy(outer_max_V=cap,accepted_inner_max_V=cap,trial_inner_max_V=.1)
            old = steps.propose(history,policy,attempt['span_V'],200000.-history[-1]['outer_V']-
                                history[-1]['base']['potential_origin_V'])
            observed_span = actual['outer_V']-old['parent']['outer_V']
            feedback = steps.feedback(old,actual['inner_V'],actual['current_A_per_um'],
                                      observed_span,policy,cost,attempt['candidate'].get('recovered',False))
            next_step = steps.propose(records[:k+1],policy,feedback['next_span_V'],200000.-actual['outer_V'])
            ceiling = steps.propose(records[:k+1],policy,cap,200000.-actual['outer_V'])
            options.append(dict(policy=asdict(policy),retrospective_prediction=old,feedback=feedback,
                                next_step=next_step,port_only_ceiling=ceiling))
        cases.append(dict(attempt=folder.name, record_index=k, actual=actual, options=options))
    result = dict(schema='vela.slot.port_step_replay.v1', solver_launched=False,
        checkpoint=str(checkpoint.resolve()), checkpoint_sha256=campaign.ll.sha256(checkpoint),
        executed_driver_sha256=expected, checked_record_indices=sorted(checked),
        replay_source_hashes=campaign.ll.hashes([Path(__file__),Path(steps.__file__),Path(campaign.__file__)]),
        validation='device, all-active rows, global closure, hashes, independent circuit confirmation',
        elapsed_note='historical elapsed is scalar-controller time; excludes outer orchestration and final confirmation',
        limitation='counterfactual port estimates only; larger steps and local convergence are unverified', cases=cases)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x',encoding='utf-8') as handle:
        json.dump(result,handle,indent=2,allow_nan=False)
    print(json.dumps(dict(output=str(output),checked_records=len(checked),cases=[dict(
        attempt=c['attempt'],outer=c['actual']['outer_V'],
        predictions=[dict(cap=o['policy']['outer_max_V'],next=o['next_step']['span_V'],
                          g=o['next_step']['fit'].get('slope_A_per_um_V'),
                          error=o['feedback']['normalized_port_prediction_error']) for o in c['options']])
        for c in cases]),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkpoint',type=Path,required=True)
    parser.add_argument('--executed-driver',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    replay(args.checkpoint,args.executed_driver,args.output)
