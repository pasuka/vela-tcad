"""Run the sealed eight-source corpus through an extracted offline SDE package."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

REPO=Path(__file__).resolve().parents[1]


def qualify(root,output,components,repeat=2,case_ids=None):
    sys.path.insert(0,str(REPO/'scripts'))
    from sentaurus_sde_cases import load_manifest,DEFAULT_MANIFEST,check_contract,compare_native
    from vela_sdevice.mesh_import import import_grid
    data=load_manifest(DEFAULT_MANIFEST)
    python=root/('bin/python/python.exe' if os.name=='nt' else 'bin/python/bin/python3')
    assert python.is_file() and repeat>=1
    output.mkdir(parents=True,exist_ok=True)
    allowed=('SystemRoot','WINDIR','TEMP','TMP','HOME','LANG')
    env={k:v for k,v in os.environ.items() if k in allowed}
    env['PATH']=str(components[2].parent)+os.pathsep+(str(Path(os.environ['SystemRoot'])/'System32') if os.name=='nt' else '/usr/bin:/bin')
    env['PYTHONPATH']=str(output/'poison')
    (output/'poison').mkdir(exist_ok=True)
    for name in ('jsonschema.py','gmsh.py'):
        (output/'poison'/name).write_text('raise RuntimeError("host import contamination")\n')
    rows=[]
    for case in data['cases']:
        if case_ids and case['id'] not in case_ids:continue
        runs=[]
        for iteration in range(repeat):
            directory=output/case['id']/str(iteration)
            cmd=[str(python),'-B','-I','-S',str(root/'run_sde.py'),'--sde',str(REPO/case['source']),
                 '--output-dir',str(directory),'--gmsh-python',str(components[0]),
                 '--gmsh-library',str(components[1]),'--gmsh-executable',str(components[2])]
            for key,value in case.get('parameters',{}).items():cmd+=['--parameter',key+'='+value]
            result=subprocess.run(cmd,cwd=output,env=env,capture_output=True,text=True)
            (output/(case['id']+'-'+str(iteration)+'.log')).write_text(result.stdout+result.stderr)
            if result.returncode:raise RuntimeError(case['id']+': '+result.stderr[-1000:])
            quality=check_contract(directory,case['contract'])
            with tempfile.TemporaryDirectory(prefix='p6-independent-import-') as temp:
                import_grid(directory/'device.grd',Path(temp));compare_native(directory,Path(temp))
            runs.append({'sha256':{name:hashlib.sha256((directory/name).read_bytes()).hexdigest()
                                  for name in ('device.grd','device.dat','device.msh')},'contract':'passed',
                         'nodes':quality['node_count'],'triangles':quality['cell_count']})
        if case['repeat_policy']=='byte_identical' and any(r['sha256']!=runs[0]['sha256'] for r in runs):
            raise RuntimeError(case['id']+': byte-identical repeat contract failed')
        rows.append({'id':case['id'],'source_sha256_lf':case['source_sha256_lf'],'runs':runs})
        print(case['id'],'passed',flush=True)
    if not rows:raise ValueError('no selected corpus cases')
    report={'schema':'vela.sde_package_qualification.v1','status':'passed','cases':rows,
            'isolated_python':True,'poisoned_pythonpath_ignored':True,'external_gmsh':True}
    (output/'qualification.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package-root',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--gmsh-python',type=Path,required=True)
    parser.add_argument('--gmsh-library',type=Path,required=True)
    parser.add_argument('--gmsh-executable',type=Path,required=True)
    parser.add_argument('--repeat',type=int,default=2)
    parser.add_argument('--case',action='append')
    args=parser.parse_args()
    qualify(args.package_root.resolve(),args.output_dir.resolve(),
            [p.resolve() for p in (args.gmsh_python,args.gmsh_library,args.gmsh_executable)],args.repeat,args.case)


if __name__=='__main__':main()
