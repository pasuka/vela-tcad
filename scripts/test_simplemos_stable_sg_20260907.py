"""Independent multiprecision derivative property regression."""
import subprocess
import build_simplemos_joint_stable_20260907 as b

def main():
    original=b.a.read(b.old.LOCAL/'build_command.json')
    flags=[x for x in original[1:] if x.startswith(('-I','-D','-O','-std='))]
    args=[original[0]]+flags+[str(b.p.REPO/'tests/diagnostics/test_simplemos_stable_sg_derivative.cpp')]
    args += [x for x in original if x.endswith('.a')]+['D:/msys64/ucrt64/lib/libCatch2Main.a','D:/msys64/ucrt64/lib/libCatch2.a','-o',str(b.LOCAL/'test_stable_sg.exe')]
    b.a.write(b.LOCAL/'test_build_command.json',args)
    r=subprocess.run(args,env=b.env(b.LOCAL),capture_output=True,text=True)
    (b.LOCAL/'test_build.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-3000:]
    r=subprocess.run([str(b.LOCAL/'test_stable_sg.exe')],env=b.env(b.LOCAL),capture_output=True,text=True)
    (b.LOCAL/'test_result.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stdout+r.stderr
    print(r.stdout,flush=True)

if __name__=='__main__':main()
