"""Compile and run bounded Catch2 geometry/conservative-flux properties."""
from pathlib import Path
import subprocess
import build_simplemos_joint_geometry_20260907 as b

def main():
    args=b.a.read(b.LOCAL/'build_command.json')
    args=[x for x in args if not x.endswith('.cpp')];args[-1]=str(b.LOCAL/'test_geometry.exe')
    index=next(i for i,x in enumerate(args) if x.endswith('libvela_core.a'))
    args.insert(index,str(b.p.REPO/'tests/diagnostics/test_simplemos_native_geometry.cpp'))
    index=args.index('-o');args[index:index]=['D:/msys64/ucrt64/lib/libCatch2Main.a','D:/msys64/ucrt64/lib/libCatch2.a']
    b.a.write(b.LOCAL/'test_build_command.json',args)
    r=subprocess.run(args,env=b.env(b.LOCAL),capture_output=True,text=True)
    (b.LOCAL/'test_build.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    r=subprocess.run([str(b.LOCAL/'test_geometry.exe')],env=b.env(b.LOCAL),capture_output=True,text=True)
    (b.LOCAL/'test_result.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stdout+r.stderr
    print(r.stdout,flush=True)

if __name__=='__main__':main()
