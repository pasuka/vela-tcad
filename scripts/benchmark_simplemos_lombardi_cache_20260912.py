"""Isolated timing of state-independent Lombardi constants; no solver edits."""
from pathlib import Path
import subprocess,os
import simplemos_enormal_curves_20260912 as e
a,d=e.a,e.d
L=e.R/'build-release/lombardi_cache_preflight_20260912';O=e.R/'reference_tcad/simplemos_sentaurus2022/lombardi_cache_preflight_20260912'
SOURCE=r'''#include "vela/equation/SplitDDOperator.h"
#include <chrono>
#include <iostream>
#include <iomanip>
using vela::split_dd::Wide;
struct Prepared {
    vela::LombardiParameters p;Wide damp,C;
    Prepared(Wide distance,Wide impurity,const vela::LombardiParameters& params):p(params) {
        damp=p.criticalLength>1.?Wide(1):Wide(exp(-distance/Wide(p.criticalLength)));
        C=Wide(p.C)*pow(Wide(300./300.),Wide(-p.k))*pow((impurity+Wide(p.N2))/Wide(p.N0),Wide(p.lambda));
    }
    Wide value(Wide field) const {
        if(field<=0)return Wide(0);
        Wide ac=field/(Wide(p.B)+C*pow(field,Wide(2)/3));
        Wide sr=pow(field/Wide(100),Wide(p.A))/Wide(p.delta)+field*field*field/Wide(p.eta);
        return Wide(damp*(Wide(p.acousticFactor)*ac+Wide(p.roughnessFactor)*sr));
    }
};
int main() {
    vela::MobilityModelConfig cfg;
    struct Item { Wide distance,impurity,field;vela::LombardiParameters p;Prepared prep; };
    std::vector<Item> items;
    for(auto p:{cfg.electronLombardi,cfg.holeLombardi})
      for(double distance:{0.,1e-10,1e-8,2e-7,1e-6})
       for(double impurity:{1e19,1e23,2e23,1e26})
        for(double field:{0.,1e-24,1e-6,1.,1e5,1e8}) {
            Prepared prep(Wide(distance),Wide(impurity),p);
            items.push_back({Wide(distance),Wide(impurity),Wide(field),p,prep});
        }
    int exact=0;Wide worst=0;
    for(const auto& v:items) {
        const Wide base=vela::elementLombardiInverse(v.field,v.distance,v.impurity,300.,v.p),next=v.prep.value(v.field);
        if(base==next)++exact;
        if(base!=0)worst=std::max(worst,Wide(abs(next/base-1)));
    }
    std::cout<<std::setprecision(17);
    std::cout<<"{\"checks\":"<<items.size()<<",\"exact\":"<<exact<<",\"max_relative\":"<<worst;
    std::array<double,2> time{};std::array<Wide,2> sums{{0,0}};
    for(int mode=0;mode<2;++mode) {
        const auto t=std::chrono::steady_clock::now();
        for(int repeat=0;repeat<4;++repeat)for(const auto& v:items)
            sums[mode]+=mode?v.prep.value(v.field):vela::elementLombardiInverse(v.field,v.distance,v.impurity,300.,v.p);
        time[mode]=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();
    }
    std::cout<<",\"uncached_seconds\":"<<time[0]<<",\"cached_seconds\":"<<time[1]<<",\"kernel_speedup\":"<<time[0]/time[1]<<",\"sums_identical\":"<<(sums[0]==sums[1]?"true":"false")<<"}\n";
    return exact==int(items.size())&&sums[0]==sums[1]?0:1;
}
'''

def main():
    assert not (O/'evidence.json').exists();L.mkdir(parents=True,exist_ok=True)
    (L/'benchmark.cpp').write_text(SOURCE,encoding='utf-8',newline='\n')
    env=os.environ.copy();env['PATH']='D:/msys64/ucrt64/bin'+os.pathsep+env['PATH']
    r=subprocess.run(['D:/msys64/ucrt64/bin/g++.exe','-std=c++20','-O2','-I'+str(e.R/'include'),'-ID:/msys64/ucrt64/include/eigen3',str(L/'benchmark.cpp'),'-o',str(L/'benchmark.exe')],env=env,capture_output=True,text=True)
    (L/'compile.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr
    r=subprocess.run([str(L/'benchmark.exe')],env=env,capture_output=True,text=True);(L/'benchmark.log').write_text(r.stdout+r.stderr)
    import json
    result=json.loads(r.stdout);result.update(scope='Kernel-only isolated cache of static distance/doping coefficients at 300 K. Current frozen solver and curves unchanged. No production migration, matrix or self-consistent qualification. Wall time observed while four solver jobs run; not a whole-solver speedup estimate.')
    a.write(O/'summary.json',result)
    files=[Path(__file__).resolve(),O/'summary.json',e.R/'include/vela/physics/ElementLombardi.h',e.R/'include/vela/equation/SplitDDOperator.h',e.R/'include/vela/numerics/SplitDDState.h',e.R/'include/vela/physics/MobilityModel.h']+[p for p in L.iterdir() if p.is_file()]
    d.matrix.freeze(O/'evidence.json',files);print(result,flush=True);assert r.returncode==0

if __name__=='__main__':main()
