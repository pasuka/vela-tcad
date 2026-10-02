"""Audit supplemental M60 field-export curves without replacing references."""
import argparse
import math
from pathlib import Path
import simplemos_hfs_cloud_20260926 as h
from sentaurus_import import parse_quoted_list,parse_values_block


def run(base):
    references={(r['case'],int(r['index'])):float(r['tight_default_Id_A_per_um']) for r in h.rows(base/'data/m60_joined.csv')}
    rows=[]
    for c in h.read(base/'data/inputs/contract.json')['cases']:
        case=c['case'];folder=base/'native/bundle'/c['device'];prefix='m60fields_'+case
        if (base/'native'/(prefix+'.exitcode')).read_text().strip()!='0':
            raise ValueError('Unsuccessful native run '+case)
        path=folder/('IdVg_'+prefix+'_des.plt');text=path.read_text()
        names=parse_quoted_list(text,'datasets');points=parse_values_block(text,len(names))
        if len(points)!=51 or len(list(folder.glob(prefix+'_state_*_des.tdr')))!=51:
            raise ValueError('Incomplete native field curve '+case)
        for index,p in enumerate(points):
            point=dict(zip(names,p))
            if abs(point['gate OuterVoltage']-index*.05)>1e-10 or abs(point['drain OuterVoltage']-c['vd'])>1e-10:
                raise ValueError('Wrong bias grid '+case)
            old=references[case,index];new=point['drain TotalCurrent']
            delta=100*(new/old-1)
            if not math.isfinite(delta):
                raise ValueError('Nonfinite current '+case)
            rows.append(dict(case=case,index=index,vg=index*.05,frozen_M60_Id_A_per_um=old,
                supplemental_Id_A_per_um=new,change_percent=delta,plt_sha256=h.sha(path)))
    if len(rows)!=816 or {(r['case'],r['index']) for r in rows}!=set(references):
        raise ValueError('Native matrix keys do not match frozen M60')
    h.csvout(base/'native_curve_comparison.csv',rows)
    worst=max(rows,key=lambda r:abs(r['change_percent']))
    result=dict(points=len(rows),complete=True,max_abs_change_percent=abs(worst['change_percent']),
        worst=worst,original_reference_replaced=False,scope='Default-observer current only; state field precision separately assessed')
    h.write(base/'native_curve_audit.json',result)
    print(result)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);run(p.parse_args().base)
