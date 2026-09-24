"""Post-timing field audit for thermal assembly reuse (off versus on).

Uses the saved 62-point full-curve pair and the same silicon physics probe.
This is not a new native-field qualification. Native electrical/thermal scoring
is performed by the staged runner. Run after timing so probes do not compete.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

import h5py
import numpy as np
from electrothermal_state import read_bound_record
from audit_templates_ldmos_joint_local_fields import carrier_metric


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full',type=Path,required=True)
    parser.add_argument('--probe',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    assert read(args.full/'summary.json')['status']=='completed'
    args.output.mkdir(parents=True,exist_ok=False)
    contract=read(Path(__file__).resolve().parents[1]/'reference_tcad/templates_ldmos_sentaurus2022/thermal/d0_joint_acceptance_20260914_v2.json')['local_fields']
    report=dict(status='running',scope=__doc__,probe_sha256=hashlib.sha256(args.probe.read_bytes()).hexdigest(),points=[])
    def save():
        (args.output/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    save()
    try:
        for gate in (4,8):
            folders=[args.full/f'r0_vg{gate}_{e}'/'results' for e in (0,1)]
            ledgers=[read(f/'ledger.json') for f in folders]
            assert all(l['status']=='complete' and len(l['exact_points'])==31 for l in ledgers)
            cfg=read_bound_record(folders[0].parent/'input.json')
            mesh=read(cfg['mesh_file']); silicon=set();oxide=set()
            for cell in mesh['triangles']:
                (silicon if cell['region_id']==0 else oxide).update(cell['node_ids'])
            ids=sorted(silicon);areas=np.array(cfg['silicon_area_m2'])[ids]
            interface=np.array([i in oxide for i in ids])
            for index in range(31):
                fields=[];results=[]
                for enabled,ledger in enumerate(ledgers):
                    point=ledger['exact_points'][index]
                    assert abs(point['bias_V']-index*40/30)<1e-9
                    result=read_bound_record(point['result']);results.append(result)
                    x=result['referenced_state_interleaved'];origin=result.get('potential_origin_V',0.)
                    states=[dict(id=i,potential_V=x[4*i]+origin,electron_qf_V=x[4*i+1],hole_qf_V=x[4*i+2],temperature_K=x[4*i+3],
                        electron_qf_reference_V=result['electron_qf_reference_V'][i]+origin,hole_qf_reference_V=result['hole_qf_reference_V'][i]+origin,
                        donors_m3=cfg['donors_m3'][i],acceptors_m3=cfg['acceptors_m3'][i]) for i in ids]
                    with tempfile.TemporaryDirectory(prefix='probe_',dir=args.output) as temp:
                        inp=Path(temp)/'input.json';inp.write_text(json.dumps(dict(states_SI=states,auger_with_generation=cfg['auger_with_generation'])),encoding='utf-8')
                        process=subprocess.run([str(args.probe.resolve()),str(inp.resolve())],capture_output=True,text=True,check=True)
                        local=json.loads(process.stdout)['results']
                    assert [r['id'] for r in local]==ids
                    values={k:np.array([r[k]['value'] for r in local]) for k in contract['carriers']+contract['band_edges']}
                    assert all(np.isfinite(v).all() for v in values.values())
                    fields.append(values)
                    with h5py.File(args.output/f'vg{gate}_p{index:02d}_{enabled}.h5','w') as f:
                        f.attrs['source_result']=point['result'];f.create_dataset('node_id',data=ids)
                        for k,v in values.items():f.create_dataset(k,data=v)
                row=dict(gate_V=gate,index=index,bias_V=index*40/30,carriers={},bands={},contacts={},pass_gate=True)
                for key in contract['carriers']:
                    row['carriers'][key]={}
                    for group,mask in (('all_silicon',np.ones(len(ids),dtype=bool)),('silicon_oxide_interface',interface)):
                        metric=carrier_metric(fields[0][key],fields[1][key],areas,mask,contract['reference_density_floor_fraction_of_carrier_peak'])
                        row['carriers'][key][group]=metric
                        row['pass_gate'] &= metric['weighted_relative_rms']<=contract['volume_weighted_relative_rms_limit']
                for key in contract['band_edges']:
                    row['bands'][key]=float(np.max(np.abs(fields[0][key]-fields[1][key])))
                    row['pass_gate'] &= row['bands'][key]<=contract['maximum_absolute_band_edge_error_eV']
                assert [c['contact'] for c in results[0]['contacts']]==[c['contact'] for c in results[1]['contacts']]
                for a,b in zip(results[0]['contacts'],results[1]['contacts']):
                    row['contacts'][a['contact']]=abs(a['total_outflow_A_per_m']-b['total_outflow_A_per_m'])
                for key in ('boundary_heat_W_per_m','lattice_source_W_per_m'):
                    row[key+'_absolute_difference']=abs(results[0][key]-results[1][key])
                row['temperature_maximum_absolute_difference_K']=max(abs(a-b) for a,b in zip(results[0]['temperature_K'],results[1]['temperature_K']))
                report['points'].append(row);save()
                if not row['pass_gate']:raise ValueError('Derived field comparison failed')
        report['status']='completed';save()
    except BaseException as error:
        report.update(status='failed',error=repr(error));save();raise


if __name__=='__main__':
    main()
