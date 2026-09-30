"""Read live baseline results across the immutable original and resumed trees."""
import json
from pathlib import Path
import status_simplemos_hfs_curves_20260914 as old


def main():
    rows = []
    for root in (old.LOCAL/'vela', old.LOCAL/'resume_20260915/vela'):
        for path in root.rglob('result.json'):
            try:
                rows.append(json.loads(path.read_text(encoding='utf-8')))
            except (FileNotFoundError, json.JSONDecodeError):
                continue
    good = {(r['case'], r['index'], r['arm']): r for r in rows if r['qualified']}
    pairs = [abs(r['current_A_per_um']/good[(case, index, 'native')]['current_A_per_um']-1)
             for (case, index, arm), r in good.items()
             if arm == 'continuation' and (case, index, 'native') in good]
    result = dict(qualified_states=len(good), target_states=408,
        finished_attempts=len(rows), failed_attempts=sum(not r['qualified'] for r in rows),
        administrative_interruptions_retained=4, paired_currents=len(pairs),
        max_paired_Id_relative=max(pairs, default=None),
        interpretation='Live current pairing only; full dual field qualification pending final analysis.',
        groups=[dict(device=device, vd=vd, arm=arm, qualified=sum(
            r['device']==device and r['vd']==vd and r['arm']==arm for r in good.values()))
            for device in ('n19', 'n23') for vd in (.05, 1.) for arm in ('continuation', 'native')],
        failures=[r for r in rows if not r['qualified']])
    print(json.dumps(result), flush=True)
    return result


if __name__ == '__main__':
    main()
