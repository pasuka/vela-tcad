// Isolated production scalar evaluation; inputs cm^-3, outputs cm^2/(V s).
// This diagnostic does not change the production element-box model guard.
#include "vela/physics/MobilityModel.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
int main(int argc,char** argv) {
    if(argc<3||argc>4)return 2;
    const bool derivatives=argc==4;
    std::ifstream in(argv[1]);std::ofstream out(argv[2]);
    if(!in||!out)return 3;
    out<<std::setprecision(17)<<"id,mu_e,mu_h,P_e,P_h,G_e,G_h";
    const std::vector<double> steps{1e-4,3e-5,1e-5,3e-6,1e-6,3e-7};
    if(derivatives)for(int c=0;c<2;++c)for(int p=0;p<2;++p)for(int k=0;k<6;++k)
        out<<",dlog_"<<c<<'_'<<p<<'_'<<k;
    out<<'\n';std::string line;std::getline(in,line);
    while(std::getline(in,line)) {
        std::istringstream row(line);std::vector<std::string> f;std::string token;
        while(std::getline(row,token,','))f.push_back(token);
        if(f.size()!=5)return 4;
        vela::PhuMobScalarState state{std::stod(f[1])*1e6,std::stod(f[2])*1e6,
                                      std::stod(f[3])*1e6,std::stod(f[4])*1e6,300.};
        const auto e=vela::evaluatePhuMobScalar(vela::CarrierType::Electron,state);
        const auto h=vela::evaluatePhuMobScalar(vela::CarrierType::Hole,state);
        out<<f[0]<<','<<e.mobility*1e4<<','<<h.mobility*1e4<<','<<e.screeningParameter
           <<','<<h.screeningParameter<<','<<e.screeningG<<','<<h.screeningG;
        if(derivatives)for(int c=0;c<2;++c)for(int p=0;p<2;++p)for(double step:steps) {
            auto plus=state,minus=state;
            if(p==0){plus.electrons*=std::exp(step);minus.electrons*=std::exp(-step);}
            else {plus.holes*=std::exp(step);minus.holes*=std::exp(-step);}
            const auto carrier=c==0?vela::CarrierType::Electron:vela::CarrierType::Hole;
            const double fp=vela::evaluatePhuMobScalar(carrier,plus).mobility*1e4;
            const double fm=vela::evaluatePhuMobScalar(carrier,minus).mobility*1e4;
            out<<','<<(fp-fm)/(2*step);
        }
        out<<'\n';
    }
}
