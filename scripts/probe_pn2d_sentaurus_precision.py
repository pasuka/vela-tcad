#!/usr/bin/env python3
"""Isolated native arithmetic-precision diagnostic; never a formal sweep result."""
import argparse
import copy
import subprocess
import tarfile
from pathlib import Path
from pn2d_variants import DEFAULT_ROOT, SPEC, read_json, write_json, sdevice_text, sha, run_cmd, tag


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=DEFAULT_ROOT)
    parser.add_argument('--id',default='D1',choices=['D1','D2'])
    parser.add_argument('--mesh',default='M2',choices=['M2','M3'])
    parser.add_argument('--precision',type=int,default=128,choices=[80,128])
    parser.add_argument('--branch',default='reverse',choices=['forward','reverse'])
    parser.add_argument('--rhs-min',type=float)
    parser.add_argument('--ssh-target',default='sentaurus')
    args=parser.parse_args();root=args.root.resolve()
    spec=copy.deepcopy(read_json(SPEC))
    cfg={**spec['baseline'],**next(e['overrides'] for e in spec['experiments'] if e['id']==args.id)}
    bias=.02 if args.branch=='forward' else -.1
    spec['branches'][args.branch]={'stop_V':bias,'step_V':bias}
    spec['saved_biases_V']=[0.,bias]
    name=f'{args.id}_{args.mesh}_EP{args.precision}'
    if args.branch!='reverse': name+='_'+args.branch
    if args.rhs_min is not None: name+=f'_rhs{args.rhs_min:g}'
    directory=root/'diagnostics/native_precision'/name;directory.mkdir(parents=True,exist_ok=True)
    numerics={**spec['meshes'][args.mesh],'sentaurus_precision_bits':args.precision}
    if args.branch=='forward': numerics['sdevice_goal_initial_step']=1.0
    if args.rhs_min is not None: numerics['sdevice_rhs_min']=args.rhs_min
    deck=sdevice_text(spec,cfg,args.branch,numerics)
    (directory/'probe.cmd').write_text(deck)
    source=root/args.id/args.mesh
    (directory/'models.par').write_bytes((source/'models.par').read_bytes())
    remote=read_json(root/'remote_root.json')['directory']
    target=f'{remote}/diagnostics/native_precision/{name}'
    ssh=['C:/Windows/System32/OpenSSH/ssh.exe','-o','BatchMode=yes','-o','ConnectTimeout=15',args.ssh_target]
    scp=['C:/Windows/System32/OpenSSH/scp.exe','-o','BatchMode=yes']
    banner=subprocess.run(ssh+['sdevice -h'],capture_output=True,text=True,timeout=45)
    if spec['sentaurus_version'] not in banner.stdout+banner.stderr: raise ValueError('Native version mismatch')
    subprocess.run(ssh+[f'mkdir -p {target} && cp {remote}/{args.id}/{args.mesh}/pn2d_msh.tdr {target}/pn2d_msh.tdr'],check=True)
    subprocess.run(scp+[str(directory/'probe.cmd'),str(directory/'models.par'),f'{args.ssh_target}:{target}/'],check=True)
    identity={'scope':'diagnostic_only','id':args.id,'mesh':args.mesh,'precision_bits':args.precision,'digits':8,
              'branch':args.branch,'bias_V':bias,'rhs_min':args.rhs_min,
              'version':spec['sentaurus_version'],'remote':target,'source_mesh_sha256':sha(source/'pn2d_msh.tdr'),
              'deck_sha256':sha(directory/'probe.cmd'),'models_sha256':sha(directory/'models.par')}
    write_json(directory/'run.json',{**identity,'status':'running'})
    rc=run_cmd(ssh+[f'cd {target} && sdevice probe.cmd > probe.out 2>&1'],root,directory/'remote.log',timeout=3600)
    subprocess.run(ssh+[f'cd {target} && tar --ignore-failed-read -czf transfer.tgz -- *.cmd *.par *.out *.log *.plt {args.branch}_p0p00_des.tdr {args.branch}_{tag(bias)}_des.tdr'],check=True)
    subprocess.run(scp+[f'{args.ssh_target}:{target}/transfer.tgz',str(directory/'transfer.tgz')],check=True)
    with tarfile.open(directory/'transfer.tgz') as archive: archive.extractall(directory,filter='data')
    write_json(directory/'run.json',{**identity,'returncode':rc,'status':'complete' if rc==0 else 'fail'})
    return rc


if __name__=='__main__':
    raise SystemExit(main())
