"""Compile the opt-in diagnostic kernel and check physical invariants.

This header is deliberately absent from production solver targets.
"""
from pathlib import Path
import os
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
HARNESS = r'''
#include "bvmethods_private_element_source.h"
#include <iostream>
#include <string>
using namespace bvmethods_diagnostic;
int main(int argc,char** argv) {
 const std::string test=argv[1];
 std::array<std::array<double,2>,3> points{{{0,0},{1e-6,0},{0,1e-6}}};
 std::array<double,3> n{1e15,1e15,1e15},psi{0,-.1,-.2},ref{},inc{0,-.1,-.2};
 std::array<std::array<double,2>,2> normals{{{1,0},{0,1}}};
 const auto close=[](long double a,long double b,long double tol=1e-12L) {return std::abs(a-b)<=tol*std::max(1e-30L,std::abs(b));};
 auto run=[&]<class T>(bool e=true,int count=0,bool contact=false) {return evaluate<T>(points,n,psi,ref,inc,100,e,normals,count,contact);};
 bool ok=false;
 if(test=="equilibrium") {inc={0,0,0};auto r=run.operator()<double>();ok=r.generation==0 && r.current[0]==0 && r.current[1]==0;}
 if(test=="uniform_drift") {auto r=run.operator()<double>();ok=close(r.current[0],1.602192e-19*1e15*.1/1e-6) && close(r.current[1],2*r.current[0]);}
 if(test=="carrier_sign") {auto e=run.operator()<double>();auto h=run.operator()<double>(false);ok=close(e.current[0],h.current[0]) && close(e.current[1],h.current[1]);}
 if(test=="counter_field") {psi={0,.1,.2};auto r=run.operator()<double>();ok=r.force==0 && r.generation==0;}
 if(test=="boundary") {auto r=run.operator()<double>(true,1,false);auto c=run.operator()<double>(true,1,true);ok=r.current[0]==0 && close(r.current[1],c.current[1]) && c.current[0]>0;}
 if(test=="corner") {auto r=run.operator()<double>(true,2,false);ok=r.current[0]==0 && r.current[1]==0 && r.generation==0;}
 if(test=="split_qf") {ref={1,1,1};inc={0,1e-20,2e-20};auto r=run.operator()<long double>();ok=r.qfDrop[0]>0 && r.current[0]!=0 && close(r.qfDrop[0],1e-20L);}
 if(test=="invalid") {n[1]=0;try{run.operator()<double>();}catch(const std::invalid_argument&){ok=true;}}
 if(test=="threshold") {long double f=400000;auto a=alpha(f,false);auto b=(1582000*std::exp(-2036000/f)+671000*std::exp(-1693000/f))/2;ok=close(a,b,1e-9L) && alpha(0.L,false)==0;}
 if(test=="permutation") {auto a=run.operator()<long double>();std::swap(points[1],points[2]);std::swap(n[1],n[2]);std::swap(psi[1],psi[2]);std::swap(inc[1],inc[2]);auto b=run.operator()<long double>();ok=close(a.current[0],b.current[0]) && close(a.current[1],b.current[1]) && close(a.generation,b.generation);}
 if(!ok){std::cerr<<"failed "<<test<<'\n';return 1;}return 0;
}
'''


class PrivateElementSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        build = ROOT / "build-release"
        build.mkdir(exist_ok=True)
        cls.directory = tempfile.TemporaryDirectory(prefix="private-element-test-", dir=build)
        d = Path(cls.directory.name)
        source = d / "test.cpp"
        source.write_text(HARNESS)
        compiler = os.environ.get("CXX") or shutil.which("g++")
        if compiler is None and os.name == "nt":
            compiler = "D:/msys64/ucrt64/bin/g++.exe"
        if not compiler or not Path(compiler).exists():
            raise RuntimeError("C++20 compiler required for diagnostic invariants")
        cls.binary = d / ("test.exe" if os.name == "nt" else "test")
        cls.environment = os.environ.copy()
        if os.name == "nt":
            cls.environment["PATH"] = str(Path(compiler).parent) + ";" + cls.environment.get("PATH", "")
        result = subprocess.run([compiler, "-std=c++20", "-O2", "-I", str(ROOT / "scripts"),
                                 str(source), "-o", str(cls.binary)], env=cls.environment,
                                capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(result.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.directory.cleanup()


def make_test(case):
    def test(self):
        subprocess.run([str(self.binary), case], env=self.environment,
                       check=True, capture_output=True)
    return test


for case in ["equilibrium", "uniform_drift", "carrier_sign", "counter_field",
             "boundary", "corner", "split_qf", "invalid", "threshold", "permutation"]:
    setattr(PrivateElementSourceTests, "test_" + case, make_test(case))

if __name__ == "__main__":
    unittest.main()
