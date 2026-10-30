"""Validate generated GRD/DAT via licensed TDX readback, independent of TDR encoding."""
import argparse
import hashlib,json
from pathlib import Path
from sentaurus_sde_cases import load_manifest,DEFAULT_MANIFEST
from sde_tdx_validation import validate_pair


def audit(original,readback,output):
    rows=[]
    for case in load_manifest(DEFAULT_MANIFEST)['cases']:
        name=case['id'];grid=readback/name/'roundtrip.grd'
        try:
            bundle=json.loads((original/name/'1/mesh_field_bundle.json').read_text())
            result=validate_pair(bundle,grid.with_suffix(''))
            rows.append({'id':name,**result,'input_sha256':{suffix:hashlib.sha256((original/name/'1'/('device.'+suffix)).read_bytes()).hexdigest() for suffix in ('grd','dat')},'readback_sha256':{
                p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (grid,grid.with_suffix('.dat'))}})
        except Exception as error:rows.append({'id':name,'status':'failed','error':str(error)})
        print(name,rows[-1]['status'],rows[-1].get('error',''),flush=True)
    report={'schema':'vela.sde_tdx_readback.v1','status':'passed' if all(r['status']=='passed' for r in rows) else 'failed',
            'scope':'generated GRD/DAT -> TDX TDR -> TDX GRD/DAT -> independent validation-only regional mapping',
            'direct_ogawa_import':'unsupported; existing HDF5 TDR reader unchanged','cases':rows}
    output.write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original-root',type=Path,required=True);parser.add_argument('--readback-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    return 0 if audit(args.original_root,args.readback_root,args.output)['status']=='passed' else 2


if __name__=='__main__':raise SystemExit(main())
