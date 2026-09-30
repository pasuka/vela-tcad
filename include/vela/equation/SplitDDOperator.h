#pragma once
#include "vela/numerics/SplitDDState.h"
#include "vela/physics/ElementLombardi.h"
#include <map>
#include <array>

namespace vela::split_dd {
// Qualification operator for the exported plain-PhuMob/Boltzmann/OldSlotboom
// geometry bundle. All state-dependent quantities consume State::potential;
// no double projection, cached nodal populations, or frozen mobilities are used.
struct PhuMob {
    std::array<Wide,5> electron,hole;
    std::array<Wide,19> c;
    std::array<Wide,2> minimumP,minimumG;
    explicit PhuMob(const nlohmann::json& j) {
        for(int i=0;i<5;++i){electron[i]=j["electron"][i].get<double>();hole[i]=j["hole"][i].get<double>();}
        for(int i=0;i<19;++i)c[i]=j["common"][i].get<double>();
        for(int b=0;b<2;++b) {
            Wide lo=-80,hi=80;
            for(int k=0;k<160;++k) {
                Wide mid=(lo+hi)/2,p=exp(mid),mass=c[4+b];
                Wide first=pow(Wide(1)/mass,c[13]),second=pow(mass,c[14]);
                Wide positive=c[10]*c[15]*first*pow(c[11]+p*first,-c[15]-1);
                Wide negative=c[12]*c[16]*pow(p*second,-c[16])/p;
                if(positive>negative)hi=mid;else lo=mid;
            }
            minimumP[b]=exp((lo+hi)/2);minimumG[b]=abs(rawG(minimumP[b],b));
        }
    }
    Wide rawG(const Wide& p,int b) const {
        return 1-c[10]*pow(c[11]+p*pow(Wide(1)/c[4+b],c[13]),-c[15])+
            c[12]*pow(p*pow(c[4+b],c[14]),-c[16]);
    }
    Wide mobility(int b,const Wide& donors,const Wide& acceptors,const Wide& n,const Wide& p) const {
        auto cluster=[](Wide x,Wide ref,Wide factor)->Wide {Wide z=x/ref;return x*(1+z*z/(1+factor*z*z));};
        Wide nd=cluster(donors*c[17],c[0]*c[17],c[2]);
        Wide na=cluster(acceptors*c[17],c[1]*c[17],c[3]);
        Wide en=n*c[17],hp=p*c[17],other=b?en:hp,sc=nd+na+other;
        const auto& m=b?hole:electron;Wide maximum=m[0]*c[18],minimum=m[1]*c[18];
        if(sc<=0)return maximum/c[18];
        Wide free=en+hp,mass=c[4+b],ratio=mass/c[5-b];
        Wide screening=1/(c[6]*pow(sc,Wide(2)/3)/Wide(3.97e13)+c[7]*free/(Wide(1.36e20)*mass));
        Wide power=pow(screening,Wide(0.6478));
        Wide f=(Wide(0.7643)*power+Wide(2.2999)+Wide(6.5502)*ratio)/(power+Wide(2.3670)-Wide(0.8552)*ratio);
        Wide g=minimumG[b];if(screening>=minimumP[b])g=std::max(g,rawG(screening,b));
        Wide effective=b?na+g*nd+c[9]*en/f:nd+g*na+c[8]*hp/f;
        Wide muN=maximum*maximum/(maximum-minimum),muC=maximum*minimum/(maximum-minimum);
        Wide scatter=muN*(sc/effective)*pow(m[3]*c[17]/sc,m[4])+muC*free/effective;
        return (1/(1/maximum+1/scatter))/c[18];
    }
};

struct Evaluation {
    std::vector<Wide> psi,phin,phip,n,p,mun,mup,srh,residual,electronFlux,holeFlux;
    std::vector<Wide> edgeMun,edgeMup;
    std::map<std::string,Wide> currents;
};
class Operator {
    nlohmann::json data_;
    PhuMob mobility_;
    std::array<LombardiParameters,2> surfaceParameters_;
    std::array<FieldMobilityParameters,2> hfsParameters_;
    static Wide bernoulli(Wide x) {return x==0?Wide(1):x/expm1(x);}
public:
    explicit Operator(nlohmann::json data):data_(std::move(data)),mobility_(data_.at("phumob")) {
        if(data_.at("schema")!="vela.split-dd-operator.v1" || data_.at("temperature_K")!=300.0)
            throw std::invalid_argument("Split operator requires its qualified 300 K schema");
        if(data_.contains("hfs"))for(int b=0;b<2;++b){hfsParameters_[b].saturationVelocity=data_["hfs"]["parameters"][b][0];hfsParameters_[b].beta=data_["hfs"]["parameters"][b][1];}
        if(data_.contains("enormal"))for(int b=0;b<2;++b) {
            const auto& j=data_["enormal"]["parameters"][b];auto& p=surfaceParameters_[b];
            p.B=j[0];p.C=j[1];p.N0=j[2];p.N2=j[3];p.lambda=j[4];p.k=j[5];p.delta=j[6];p.A=j[7];p.eta=j[8];p.criticalLength=j[9];p.acousticFactor=j[10];p.roughnessFactor=j[11];
        }
    }
    Evaluation evaluate(const State& state) const {
        const std::size_t N=state.size();
        if(N!=data_["nodes"].size() || state.meshFingerprint()!=data_.at("mesh_fingerprint").get<std::string>())
            throw std::invalid_argument("Split operator and checkpoint mismatch");
        // Encoded Dirichlet rows belong to this particular coordinate frame.
        // A matching mesh alone cannot qualify a different scale or QF origin.
        if(state.potentialScale()!=data_.at("potential_scale").get<double>())
            throw std::invalid_argument("Split operator potential scale mismatch");
        for(std::size_t i=0;i<N;++i)
            if(state.reference(1,i)!=data_["nodes"][i][7].get<double>() ||
               state.reference(2,i)!=data_["nodes"][i][8].get<double>())
                throw std::invalid_argument("Split operator QF reference mismatch");
        Evaluation o;
        for(auto* v:{&o.psi,&o.phin,&o.phip,&o.n,&o.p,&o.mun,&o.mup,&o.srh})v->resize(N);
        o.residual.resize(3*N);std::vector<bool> conducting(N,false);
        for(std::size_t i=0;i<N;++i) {
            o.psi[i]=state.potential(0,i);o.phin[i]=state.potential(1,i);o.phip[i]=state.potential(2,i);
            const auto& a=data_["nodes"][i];Wide ni=a[9].get<double>(),vt=a[10].get<double>();
            if(ni<=0)continue;
            Wide zn=(o.psi[i]-o.phin[i])/vt,zp=(o.phip[i]-o.psi[i])/vt;
            if(abs(zn)>=500 || abs(zp)>=500)throw std::invalid_argument("Split operator clipping branch is unqualified");
            o.n[i]=ni*exp(zn);o.p[i]=ni*exp(zp);
            const auto& material=data_["node_physics"][i];
            o.mun[i]=mobility_.mobility(0,Wide(material[0].get<double>()),Wide(material[1].get<double>()),o.n[i],o.p[i]);
            o.mup[i]=mobility_.mobility(1,Wide(material[0].get<double>()),Wide(material[1].get<double>()),o.n[i],o.p[i]);
            Wide numerator=ni*ni*expm1((o.phip[i]-o.phin[i])/vt);
            o.srh[i]=numerator/(Wide(material[3].get<double>())*(o.n[i]+ni)+Wide(material[2].get<double>())*(o.p[i]+ni));
            Wide source=o.srh[i]*Wide(material[4].get<double>());
            o.residual[N+i]+=source;o.residual[2*N+i]+=source;
        }
        for(const auto& e:data_["poisson_edges"]) {
            int i=e[0],j=e[1];Wide f=Wide(e[2].get<double>())*(o.psi[i]-o.psi[j]);
            o.residual[i]+=f;o.residual[j]-=f;
        }
        for(std::size_t i=0;i<N;++i) {
            const auto& a=data_["nodes"][i];
            o.residual[i]=(o.residual[i]+Wide(a[15].get<double>())*Wide(a[16].get<double>())*
                (o.n[i]*Wide(a[12].get<double>())-o.p[i]*Wide(a[13].get<double>())-Wide(a[11].get<double>())*Wide(a[14].get<double>()))-Wide(a[18].get<double>()))/Wide(a[17].get<double>());
        }
        std::vector<Wide> cellMun,cellMup;
        if(data_.contains("enormal")) {
            const auto& surface=data_["enormal"];Wide muFactor=surface["mobility_to_SI"].get<double>(),densityFactor=surface["concentration_to_SI"].get<double>();
            for(const auto& cell:surface["cells"]) {
                int i=cell["nodes"][0],j=cell["nodes"][1],k=cell["nodes"][2];
                Wide F=abs((o.psi[j]-o.psi[i])*Wide(cell["stencil"][0].get<double>())+(o.psi[k]-o.psi[i])*Wide(cell["stencil"][1].get<double>()));
                std::array<Wide,2> drive{{0,0}};
                if(data_.contains("hfs"))for(int b=0;b<2;++b) {
                    const auto& z=cell["hfs_contact"].get<bool>()?o.psi:(b?o.phip:o.phin);
                    const auto& s=cell["hfs_stencil"];
                    Wide gx=(z[j]-z[i])*Wide(s[0].get<double>())+(z[k]-z[i])*Wide(s[2].get<double>());
                    Wide gy=(z[j]-z[i])*Wide(s[1].get<double>())+(z[k]-z[i])*Wide(s[3].get<double>());
                    drive[b]=sqrt(gx*gx+gy*gy);
                }
                std::array<Wide,2> mu{{0,0}};
                for(int v=0;v<3;++v) {
                    int n=cell["nodes"][v];Wide distance=cell["distances_m"][v].get<double>();
                    Wide impurity=(Wide(data_["node_physics"][n][0].get<double>())+Wide(data_["node_physics"][n][1].get<double>()))*densityFactor;
                    for(int b=0;b<2;++b) {
                        Wide bulk=(b?o.mup[n]:o.mun[n])*muFactor;
                        Wide inv=elementLombardiInverse(F,distance,impurity,300.,surfaceParameters_[b]);
                        Wide low=1/(1/bulk+inv)/muFactor;
                        if(data_.contains("hfs"))low=elementCanali(low,drive[b],hfsParameters_[b],data_["hfs"]["cutoff"].get<double>());
                        mu[b]+=Wide(cell["weights"][v].get<double>())*low;
                    }
                }
                cellMun.push_back(mu[0]);cellMup.push_back(mu[1]);
            }
        }
        for(const auto& e:data_["transport_edges"]) {
            int i=e["i"],j=e["j"];Wide en=0,hp=0,mn=0,mp=0;
            if(!e["weights"].empty()) {

        if(data_.contains("enormal"))for(const auto& w:e["cell_weights"]) {int k=w[0];mn+=Wide(w[1].get<double>())*cellMun[k];mp+=Wide(w[1].get<double>())*cellMup[k];}
                else for(const auto& w:e["weights"]) {int k=w[0];mn+=Wide(w[1].get<double>())*o.mun[k];mp+=Wide(w[1].get<double>())*o.mup[k];}
                Wide vt=data_["nodes"][i][10].get<double>(),ni0=data_["nodes"][i][9].get<double>(),ni1=data_["nodes"][j][9].get<double>();
                Wide eta=(o.psi[j]-o.psi[i])/vt+log(ni1/ni0);
                Wide geo=e["coefficient"].get<double>();
                en=mn*geo*bernoulli(eta)*o.n[j]*expm1((o.phin[j]-o.phin[i])/vt);
                // Hole ni-gradient sign is opposite to the electron drift term.
                Wide etaH=(o.psi[j]-o.psi[i])/vt+log(ni0/ni1);
                hp=mp*geo*bernoulli(-etaH)*o.p[j]*expm1((o.phip[i]-o.phip[j])/vt);
                conducting[i]=conducting[j]=true;
            }
            o.electronFlux.push_back(en);o.holeFlux.push_back(hp);
            o.edgeMun.push_back(mn);o.edgeMup.push_back(mp);
            o.residual[N+i]+=en;o.residual[N+j]-=en;
            o.residual[2*N+i]+=hp;o.residual[2*N+j]-=hp;
            // ContactCurrent defines total = electronCurrent - holeCurrent,
            // with each particle flux opposite to its continuity flux.
            for(const auto& [name,sign]:e["ports"].items())o.currents[name]+=Wide(sign.get<double>())*(hp-en)*Wide(data_["current_factor"].get<double>());
        }
        for(std::size_t i=0;i<N;++i)if(!conducting[i]) {
            // Canonical physical-zero inactive-carrier gauge. The old double
            // operator rounds reference/scale before addition; these gauge
            // rows are compared separately from physical carrier equations.
            o.residual[N+i]=o.phin[i]/Wide(data_["potential_scale"].get<double>());
            o.residual[2*N+i]=o.phip[i]/Wide(data_["potential_scale"].get<double>());
            if(data_.value("gauge_policy",std::string{})=="encoded_reference") {
                const double scale=data_["potential_scale"];
                o.residual[N+i]=state.coordinate(N+i)+Wide(data_["nodes"][i][7].get<double>()/scale);
                o.residual[2*N+i]=state.coordinate(2*N+i)+Wide(data_["nodes"][i][8].get<double>()/scale);
            }
        }
        for(int b=0;b<3;++b)for(const auto& bc:data_["boundary"][b]) {
            int i=bc[0];o.residual[b*N+i]=state.coordinate(b*N+i)-Wide(bc[1].get<double>());
        }
        return o;
    }
};
}
