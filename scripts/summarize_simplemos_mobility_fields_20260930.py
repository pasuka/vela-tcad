"""Summarize an exact 16 x 51 fixed-state mobility audit without changing gates."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

METRICS=('projected_cell_to_split_edge','native_seed_cell_vs_native_cell',
         'vela_cell_vs_native_seed_cell','vela_state_cell_vs_native_state_cell',
         'split_edge_vs_native_cell_reconstruction','incommensurate_node_average_to_edge')

def summarize(results):
    expected={(f'n{n}_vd_{vd}',i) for n in range(17,25) for vd in ('0p05','1') for i in range(51)}
    keys=[(r['case'],r['index']) for r in results]
    if len(keys)!=len(set(keys)) or set(keys)!=expected:
        raise ValueError('Expected unique complete 816-point coverage')
    out={'points':len(results),'curves':16,'new_dc_solves':0,'gate_changes':False,'carriers':{},
         'native_load_identity_max':{key:max(r['native_load_identity'][key] for r in results)
                                     for key in results[0]['native_load_identity']}}
    for car in ('electron','hole'):
        block={}
        for metric in METRICS:
            r=max(results,key=lambda r:r['carriers'][car][metric]['max'])
            block[metric]={'max':r['carriers'][car][metric]['max'],
                          'case':r['case'],'index':r['index'],'vg':r['vg']}
        failures=[{'case':r['case'],'index':r['index'],'vg':r['vg']}
                  for r in results if not r['carriers'][car]['same_export_cell_gate_1e_7_passed']]
        block['same_export_cell_gate_1e_7']={'passed_points':len(results)-len(failures),'failed_points':len(failures),'failures':failures}
        worst=max(results,key=lambda r:r['carriers'][car]['native_seed_cell_vs_native_cell']['max'])
        block['same_export_worst_cell']=worst['carriers'][car]['same_export_worst_cell']
        block['node_display_hypotheses_global_max']={mode:max(r['carriers'][car]['native_node_display_hypotheses'][mode]['max'] for r in results)
                                                   for mode in ('equal','cell_area','box_vertex')}
        out['carriers'][car]=block
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--summary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();source=json.loads(args.summary.read_text(encoding='utf-8-sig'))
    assert source['points']==816 and source['new_dc_solves']==0 and source['gate_changes'] is False
    result=summarize(source['results']);result['source_sha256']=hashlib.sha256(args.summary.read_bytes()).hexdigest()
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with (args.output/'points.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['case','index','vg','carrier','same_export_cell_gate_1e_7_passed']+list(METRICS));writer.writeheader()
        for r in source['results']:
            for car,b in r['carriers'].items():
                writer.writerow(dict(case=r['case'],index=r['index'],vg=r['vg'],carrier=car,same_export_cell_gate_1e_7_passed=b['same_export_cell_gate_1e_7_passed'],**{k:b[k]['max'] for k in METRICS}))
    print(json.dumps({c:{k:v for k,v in b.items() if k!='same_export_cell_gate_1e_7'} for c,b in result['carriers'].items()},indent=2))

if __name__=='__main__':main()
