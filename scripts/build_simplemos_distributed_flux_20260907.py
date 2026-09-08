"""Isolated, state-independent distributed edge source with direct port term."""
from pathlib import Path
import subprocess
import prepare_simplemos_masetti_runtime_20260907 as p
import validate_simplemos_masetti_channel_20260907 as channel
import validate_simplemos_local_conservative_flux as old
a=p.a;LOCAL=p.LOCAL/'distributed_response';RUNNER=LOCAL/'runner.exe'
FILE_ENV='VELA_DIAGNOSTIC_DISTRIBUTED_ELECTRON_FLUX_FILE'
SCALE_ENV='VELA_DIAGNOSTIC_DISTRIBUTED_ELECTRON_FLUX_SCALE'

HELPER=r'''
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <map>
#include <stdexcept>
namespace {
double diagnosticConstantElectronEdgeFlux(std::size_t edge) {
    static const std::map<std::size_t,double> source = [] {
        std::map<std::size_t,double> values;
        const char* path=std::getenv("VELA_DIAGNOSTIC_DISTRIBUTED_ELECTRON_FLUX_FILE");
        const char* scale=std::getenv("VELA_DIAGNOSTIC_DISTRIBUTED_ELECTRON_FLUX_SCALE");
        if (!path && !scale) return values;
        if (!path || !scale) throw std::runtime_error("Incomplete distributed edge source");
        const double multiplier=std::stod(scale);
        if (!std::isfinite(multiplier)) throw std::runtime_error("Nonfinite source scale");
        std::ifstream stream(path);
        if (!stream) throw std::runtime_error("Cannot read distributed edge source");
        long long id;double value;
        while (stream >> id >> value) {
            if (id<0 || !std::isfinite(value) || values.count(static_cast<std::size_t>(id)))
                throw std::runtime_error("Invalid or duplicate source edge");
            values.emplace(static_cast<std::size_t>(id),multiplier*value);
        }
        if (!stream.eof()) throw std::runtime_error("Malformed distributed source");
        return values;
    }();
    auto item=source.find(edge);
    return item==source.end() ? 0.0 : item->second;
}
}
'''

def env(root,path=None,scale=0.):
 e=channel.env(root)
 e.pop(FILE_ENV,None);e.pop(SCALE_ENV,None)
 if path is not None:e[FILE_ENV]=str(path);e[SCALE_ENV]=format(scale,'.17g')
 return e

def main():
 a.verify(p.prior.OUT/'validation_evidence.json');LOCAL.mkdir(parents=True,exist_ok=False)
 assembler=channel.LOCAL/'CoupledDDAssembler.cpp';contact=old.LOCAL/'ContactCurrent.cpp'
 for src in (assembler,contact):
  text=src.read_text();assert text.count(old.HELPER)==1
  (LOCAL/src.name).write_text(text.replace(old.HELPER,HELPER),newline='\n')
 args=a.read(channel.LOCAL/'build_command.json')
 args=[str(LOCAL/'CoupledDDAssembler.cpp') if x==str(assembler) else str(LOCAL/'ContactCurrent.cpp') if x==str(contact) else x for x in args];args[-1]=str(RUNNER)
 a.write(LOCAL/'build_command.json',args)
 r=subprocess.run(args,env=env(LOCAL),capture_output=True,text=True)
 (LOCAL/'build.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stderr[-2500:]
 print('Built isolated distributed source runner; production source unchanged',flush=True)

if __name__=='__main__':main()
