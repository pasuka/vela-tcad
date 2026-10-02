"""Archive an authorized task-local PMI build, without running SDevice."""
import argparse,shutil,tarfile
from pathlib import Path
import prepare_simplemos_hfs_observer_v2_20260912 as p
a,d=p.a,p.d
L=p.L/'native_build';O=p.O/'native_build'
REMOTE='/tmp/vela_simplemos_hfs_observer_build_20260912'


def prepare():
    a.verify(p.O/'evidence.json');assert a.read(p.O/'summary.json')['all_qualified']
    assert not (O/'input_evidence.json').exists()
    dest=L/'input';(dest/'gcc_compat').mkdir(parents=True,exist_ok=False)
    files=[Path(__file__).resolve(),p.O/'evidence.json']
    for name in ('canali_diagnostic.h','pmi_vela_hfs_observer.C'):
        shutil.copyfile(p.L/name,dest/name);files.extend([p.L/name,dest/name])
    src=p.e.L/'native_raw/pmi/gcc_compat/g++';shutil.copyfile(src,dest/'gcc_compat/g++');files.extend([src,dest/'gcc_compat/g++'])
    text='''#!/bin/bash
set -u
cd "$(dirname "$0")" || exit 90
test ! -f compile.exit_code.txt || exit 91
chmod +x gcc_compat/g++
export PATH="$PWD/gcc_compat:$PATH"
/usr/bin/g++ --version > compiler.txt
/atctools/Synopsys/tcad/T-2022.03/bin/cmi -v pmi_vela_hfs_observer.C > compile.log 2>&1
code=$?
printf '%s\\n' "$code" > compile.exit_code.txt
sha256sum canali_diagnostic.h pmi_vela_hfs_observer.C gcc_compat/g++ > input_hashes.txt
if [ "$code" -eq 0 ]; then sha256sum pmi_vela_hfs_observer.so.linux64 > binary_hash.txt; fi
tar czf ../results.tgz .
exit "$code"
'''
    (dest/'compile.sh').write_text(text,newline='\n');files.append(dest/'compile.sh')
    a.write(O/'contract.json',dict(remote_root=REMOTE,scope='CMI compile only. No mesh/state upload, no SDevice run, no production changes or model release.',native_identity_required_before_observer_analysis=True))
    d.matrix.freeze(O/'input_evidence.json',files+[O/'contract.json'])
    with tarfile.open(L/'input.tgz','w:gz') as t:t.add(dest,arcname='input')
    print(L/'input.tgz',a.sha(L/'input.tgz'),flush=True)


def unpack(expected):
    a.verify(O/'input_evidence.json');assert expected and a.sha(L/'results.tgz')==expected
    root=L/'raw';assert not root.exists()
    with tarfile.open(L/'results.tgz','r:gz') as t:
        for m in t.getmembers():
            assert (root/m.name).resolve().is_relative_to(root.resolve()) and not m.issym() and not m.islnk()
        t.extractall(root,filter='data')
    for source in (L/'input').rglob('*'):
        if source.is_file():assert a.sha(source)==a.sha(root/source.relative_to(L/'input'))
    code=int((root/'compile.exit_code.txt').read_text())
    summary=dict(exit_code=code,compiled=code==0,native_loaded=False,native_identity_qualified=False)
    if code==0:
        summary['binary_sha256']=a.sha(root/'pmi_vela_hfs_observer.so.linux64')
        assert summary['binary_sha256']==(root/'binary_hash.txt').read_text().split()[0]
    a.write(O/'summary.json',summary)
    d.matrix.freeze(O/'evidence.json',[Path(__file__).resolve(),O/'input_evidence.json',L/'results.tgz',O/'summary.json']+[f for f in root.rglob('*') if f.is_file()])
    print(summary,flush=True);assert code==0


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=('prepare','unpack'));ap.add_argument('--sha');args=ap.parse_args()
    if args.action=='prepare':prepare()
    else:unpack(args.sha)
