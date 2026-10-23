#include "vela/physics/NonlocalPathPhononKernel.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <boost/math/special_functions/airy.hpp>

using namespace vela;
using namespace vela::experimental;
using Catch::Approx;
namespace {
using Big = boost::multiprecision::cpp_dec_float_50;
constexpr Real q = 1.602192e-19, kb = 1.380662e-23;
NonlocalPathBandToBandConfig config() {
    NonlocalPathBandToBandConfig c;
    c.regions = {1}; c.constants = TunnelingConstants::SentaurusT202203;
    c.ab = IndirectABParameters{4e20, 1.9e9, 1.1082092307692308 * q};
    c.phononEnergy_J = .037 * q; c.maximumPathLength_m = 1;
    return c;
}
TunnelingBandProfile uniform(Real field = 1e8, Real temperature = 300, Index pieces = 1) {
    const Real gap = config().ab->referenceGap_J, length = gap / (q * field);
    TunnelingBandProfile p;
    p.initialValenceGradient_J_per_m = q * field;
    for (Index i = 0; i <= pieces; ++i) {
        const Real u = Real(i) / pieces;
        p.knots.push_back({u * length, -gap * u, gap * (1 - u), temperature});
    }
    return p;
}
TunnelingReservoirs reservoirs() { return {-.5*q, .5*q, 100*q, 100*q}; }
TunnelingReservoirs saturated() { return {-100*q, 100*q, 100*q, 100*q}; }
NonlocalPathBandToBandConfig microscopic() {
    auto c = config(); c.ab.reset();
    c.microscopic = IndirectMicroscopicParameters{.10208356687911158 * 9.109534e-31,
        .12827262005550075 * 9.109534e-31, 6.456560487646502e-19};
    return c;
}
struct Oracle { Big action, iv, ic, logRate; };
// Independent 50-digit numerical quadrature and bisection, NOT production
// primitives, root formula or parameter conversion. Table constants are explicit.
Oracle oracle(const TunnelingBandProfile& p, const TunnelingReservoirs& r,
              const NonlocalPathBandToBandConfig& c) {
    const Big h("6.626176e-34"), boltzmann("1.380662e-23"), pi = acos(Big(-1));
    const Big hb = h / (2*pi), mc(c.microscopic->conductionMass_kg), mv(c.microscopic->valenceMass_kg);
    const Big mr = 1/(1/mc+1/mv), energy = Big(p.knots.front().valence_J)-Big(c.valenceOffset_J);
    const Big total(p.knots.back().distance_m);
    auto bands = [&](const Big& x) {
        Index i=1; while(i+1<p.knots.size() && Big(p.knots[i].distance_m)<x)++i;
        const auto& a=p.knots[i-1]; const auto& b=p.knots[i];
        const Big u=(x-Big(a.distance_m))/(Big(b.distance_m)-Big(a.distance_m));
        const Big ev=Big(a.valence_J)+(Big(b.valence_J)-Big(a.valence_J))*u;
        const Big ec=Big(a.conduction_J)+(Big(b.conduction_J)-Big(a.conduction_J))*u;
        return std::pair<Big,Big>{energy-ev+Big(c.valenceOffset_J),ec+Big(c.conductionOffset_J)-energy};
    };
    Big lower=0,upper=total;
    for(int i=0;i<180;++i) {
        const Big x=(lower+upper)/2;const auto [a,b]=bands(x);
        if(mv*a<mc*b)lower=x;else upper=x;
    }
    const Big crossover=(lower+upper)/2;
    Big action=0,iv=0,ic=0,prefactor=0;
    auto quad = [&](auto function) {
        return boost::math::quadrature::gauss_kronrod<Big,31>::integrate(function,Big(0),Big(1),10,Big("1e-27"));
    };
    for(Index i=1;i<p.knots.size();++i) {
        const Big a(p.knots[i-1].distance_m),b(p.knots[i].distance_m);
        for(bool valence:{true,false}) {
            const Big lo=valence?a:std::max(a,crossover),hi=valence?std::min(b,crossover):b;
            if(hi<=lo)continue;
            // Square-coordinate substitution removes the true endpoint root.
            auto integrand = [&](const Big& t,bool reciprocal)->Big {
                const Big x=valence?Big(lo+(hi-lo)*t*t):Big(hi-(hi-lo)*t*t);
                const auto [v,w]=bands(x);
                const Big wave=sqrt(2*(valence?mv:mc)*(valence?v:w))/hb;
                return reciprocal?Big(2*(hi-lo)*t/wave):Big(2*(hi-lo)*t*wave);
            };
            const Big s=quad([&](const Big& t)->Big{return integrand(t,false);});
            const Big inv=quad([&](const Big& t)->Big{return integrand(t,true);});
            action+=s;if(valence)iv+=inv;else ic+=inv;
        }
        prefactor+=quad([&](const Big& t)->Big {
            const auto [v,w]=bands(a+(b-a)*t);const Big gap=v+w;
            const Big N=1/(exp(Big(c.phononEnergy_J)/(boltzmann*Big(p.knots[0].temperature_K)))-1);
            return Big((b-a)*Big(c.microscopic->coupling_J2_m_per_kg)*(1+2*N)/(64*pi*pi*Big(c.phononEnergy_J)*gap)*
                sqrt(mv*mc/(h*total*sqrt(2*mr*gap))));
        });
    }
    const Big temperature(p.knots[0].temperature_K);
    const Big fp=1/(1+exp((energy-Big(r.holeFermi_J))/(boltzmann*temperature)));
    const Big fn=1/(1+exp((energy-Big(r.electronFermi_J))/(boltzmann*temperature)));
    const Big wv=1-exp(-2*mv*Big(r.valenceAvailableEnergy_J)*iv/(hb*hb));
    const Big wc=1-exp(-2*mc*Big(r.conductionAvailableEnergy_J)*ic/(hb*hb));
    const Big logRate=log(Big(p.initialValenceGradient_J_per_m)*prefactor/iv/ic*wv*wc*abs(fp-fn))-2*action;
    return {action,iv,ic,logRate};
}
}
TEST_CASE("NLP N3 uniform field recovers independently normalized Kane Keldysh law", "[nonlocalpath][phonon]") {
    for(Real field:{1e7,3e7,1e8,5e8}) for(Real ratio:{0.,.5,2.}) {
        auto c=config();c.ab->massRatio=ratio;
        const auto r=evaluatePhononTunneling(uniform(field),saturated(),c);
        const Real reference=c.ab->A_per_m3_s*std::pow(field/100,2.5)*std::exp(-c.ab->B_V_per_m/field);
        REQUIRE(r.generation_per_m3_s/reference==Approx(1).epsilon(1e-8));
        REQUIRE(2*r.action==Approx(c.ab->B_V_per_m/field).epsilon(2e-13));
        REQUIRE(r.valenceWindowFactor==1);REQUIRE(r.conductionWindowFactor==1);
    }
}
TEST_CASE("NLP N3 nonuniform P1 action reciprocal and prefactor agree with 50 digit quadrature", "[nonlocalpath][phonon]") {
    auto p=uniform(8e7,400,4);auto c=microscopic();auto r=reservoirs();
    p.knots[1].conduction_J*=1.12;p.knots[2].valence_J*=.92;p.knots[3].conduction_J*=.88;
    for(Real window:{100.,.03,1e-12}) {
        r.valenceAvailableEnergy_J=window*q;r.conductionAvailableEnergy_J=.7*window*q;
        const auto actual=evaluatePhononTunneling(p,r,c);const auto reference=oracle(p,r,c);
        std::ostringstream record;record<<std::setprecision(17)<<"N3_ORACLE window_eV="<<window
            <<" action_relative_error="<<actual.action/reference.action.convert_to<Real>()-1
            <<" Iv_relative_error="<<actual.valenceReciprocalIntegral_m2/reference.iv.convert_to<Real>()-1
            <<" Ic_relative_error="<<actual.conductionReciprocalIntegral_m2/reference.ic.convert_to<Real>()-1
            <<" source_relative_error="<<std::expm1(actual.logAbsoluteGeneration-reference.logRate.convert_to<Real>());
        std::cout<<record.str()<<'\n';
        REQUIRE(actual.action/reference.action.convert_to<Real>()==Approx(1).epsilon(1e-8));
        REQUIRE(actual.valenceReciprocalIntegral_m2/reference.iv.convert_to<Real>()==Approx(1).epsilon(1e-8));
        REQUIRE(actual.conductionReciprocalIntegral_m2/reference.ic.convert_to<Real>()==Approx(1).epsilon(1e-8));
        REQUIRE(std::expm1(actual.logAbsoluteGeneration-reference.logRate.convert_to<Real>())==Approx(0).margin(1e-8));
    }
}
TEST_CASE("NLP N3 equilibrium recombination and tiny Fermi separation preserve sign", "[nonlocalpath][phonon]") {
    auto r=reservoirs();const auto p=uniform();const auto c=config();
    r.electronFermi_J=r.holeFermi_J=0;
    const auto equilibrium=evaluatePhononTunneling(p,r,c);
    REQUIRE(equilibrium.generation_per_m3_s==0);REQUIRE(equilibrium.sign==0);REQUIRE_FALSE(equilibrium.underflow);
    r.electronFermi_J=-.01*q;r.holeFermi_J=.01*q;
    const auto positive=evaluatePhononTunneling(p,r,c);
    std::swap(r.electronFermi_J,r.holeFermi_J);
    const auto negative=evaluatePhononTunneling(p,r,c);
    REQUIRE(positive.generation_per_m3_s>0);REQUIRE(negative.generation_per_m3_s==-positive.generation_per_m3_s);
    r.electronFermi_J=0;r.holeFermi_J=1e-30*kb*300;
    const auto tiny=evaluatePhononTunneling(p,r,c),full=evaluatePhononTunneling(p,saturated(),c);
    REQUIRE(tiny.sign==1);REQUIRE(tiny.generation_per_m3_s/full.generation_per_m3_s==Approx(2.5e-31).epsilon(2e-12).margin(0));
}
TEST_CASE("NLP N3 finite transverse windows vanish independently and use expm1 tails", "[nonlocalpath][phonon]") {
    auto r=saturated();const auto p=uniform();const auto c=config();
    const auto full=evaluatePhononTunneling(p,r,c);
    for(bool valence:{true,false}) {
        auto small=r;(valence?small.valenceAvailableEnergy_J:small.conductionAvailableEnergy_J)=1e-15*q;
        const auto finite=evaluatePhononTunneling(p,small,c);
        REQUIRE(finite.generation_per_m3_s>0);REQUIRE(finite.generation_per_m3_s<full.generation_per_m3_s);
        (valence?small.valenceAvailableEnergy_J:small.conductionAvailableEnergy_J)=0;
        REQUIRE(evaluatePhononTunneling(p,small,c).generation_per_m3_s==0);
    }
}
TEST_CASE("NLP N3 temperature and phonon parameter controls retain reference 300K calibration", "[nonlocalpath][phonon]") {
    auto c=config();const auto r=saturated();
    const auto base=evaluatePhononTunneling(uniform(),r,c);
    const auto warm=evaluatePhononTunneling(uniform(1e8,400),r,c);
    const auto thermal=[&](Real T){return 1+2/std::expm1(c.phononEnergy_J/(kb*T));};
    REQUIRE(warm.generation_per_m3_s/base.generation_per_m3_s==Approx(thermal(400)/thermal(300)).epsilon(2e-13));
    for(auto [phonon,ratio]:{std::pair{.0185,1.0459893033702452},std::pair{.074,.9007866902069033}}) {
        c.phononEnergy_J=phonon*q;
        REQUIRE(evaluatePhononTunneling(uniform(),r,c).generation_per_m3_s/base.generation_per_m3_s==Approx(1).epsilon(2e-13));
        REQUIRE(evaluatePhononTunneling(uniform(1e8,400),r,c).generation_per_m3_s/warm.generation_per_m3_s==Approx(ratio).epsilon(1e-8));
    }
}
TEST_CASE("NLP N3 retains logarithms for subnormal and underflowed sources", "[nonlocalpath][phonon]") {
    auto r=saturated();r.holeFermi_J=-800*kb*300;r.electronFermi_J=-801*kb*300;
    const auto subnormal=evaluatePhononTunneling(uniform(),r,config());
    REQUIRE(subnormal.generation_per_m3_s>0);REQUIRE(subnormal.generation_per_m3_s<std::numeric_limits<Real>::min());
    const auto reference=oracle(uniform(),r,microscopic());
    REQUIRE(std::abs(subnormal.logAbsoluteGeneration-reference.logRate.convert_to<Real>())<1e-8);
    // cpp_dec_float's direct double conversion flushes this subnormal to zero.
    // Round in units of the binary64 subnormal quantum before narrowing.
    const Big oracleUlps=exp(reference.logRate)/Big(std::numeric_limits<Real>::denorm_min());
    const Real roundedUlps=floor(oracleUlps+Big("0.5")).convert_to<Real>();
    REQUIRE(subnormal.generation_per_m3_s/std::numeric_limits<Real>::denorm_min()==roundedUlps);
    const auto underflow=evaluatePhononTunneling(uniform(1e4),saturated(),config());
    REQUIRE(underflow.generation_per_m3_s==0);REQUIRE(underflow.underflow);REQUIRE(underflow.sign==1);
    REQUIRE(std::isfinite(underflow.logAbsoluteGeneration));REQUIRE(underflow.logAbsoluteGeneration< -10000);
    const auto c=config();
    REQUIRE(std::abs(underflow.logAbsoluteGeneration-(std::log(c.ab->A_per_m3_s)+2.5*std::log(1e4/100)-c.ab->B_V_per_m/1e4))<1e-8);
}
TEST_CASE("NLP N3 both constant profiles preserve the physical electric field normalization", "[nonlocalpath][phonon]") {
    for(auto constants:{TunnelingConstants::Codata2018,TunnelingConstants::SentaurusT202203}) {
        auto c=config();c.constants=constants;auto p=uniform();
        const Real charge=tunnelingElementaryCharge(constants);
        p.initialValenceGradient_J_per_m=charge*1e8;
        for(auto& k:p.knots)k.distance_m*=q/charge;
        const auto r=evaluatePhononTunneling(p,saturated(),c);
        const Real expected=c.ab->A_per_m3_s*std::pow(1e8/100,2.5)*std::exp(-c.ab->B_V_per_m/1e8);
        REQUIRE(r.generation_per_m3_s/expected==Approx(1).epsilon(1e-8));
    }
}
TEST_CASE("NLP N3 gauge and mesh partition preserve action and signed source", "[nonlocalpath][phonon]") {
    const auto c=config();const auto r=reservoirs();const auto base=evaluatePhononTunneling(uniform(),r,c);
    for(Index pieces:{2,8,32,128}) {
        auto p=uniform(1e8,300,pieces);auto shifted=r;
        for(auto& k:p.knots){k.valence_J+=3*q;k.conduction_J+=3*q;}
        shifted.electronFermi_J+=3*q;shifted.holeFermi_J+=3*q;
        const auto actual=evaluatePhononTunneling(p,shifted,c);
        REQUIRE(actual.generation_per_m3_s/base.generation_per_m3_s==Approx(1).epsilon(2e-12));
    }
}
TEST_CASE("NLP N3 runtime bandgap and parameter responses do not refit the reference", "[nonlocalpath][phonon]") {
    auto c=config();const auto p=uniform();const auto r=saturated();const auto base=evaluatePhononTunneling(p,r,c);
    c.ab->A_per_m3_s*=2;
    REQUIRE(evaluatePhononTunneling(p,r,c).generation_per_m3_s/base.generation_per_m3_s==Approx(2).epsilon(2e-13));
    c=config();c.ab->B_V_per_m*=1.1;
    REQUIRE(evaluatePhononTunneling(p,r,c).generation_per_m3_s<base.generation_per_m3_s);
    c=config();auto reduced=p;
    for(auto& k:reduced.knots)k.conduction_J*=.95;
    REQUIRE(evaluatePhononTunneling(reduced,r,c).generation_per_m3_s>base.generation_per_m3_s);
}
TEST_CASE("NLP N3 invalid unsupported and exhausted profiles fail explicitly", "[nonlocalpath][phonon]") {
    const auto p=uniform(1e8,300,4);const auto c=config();const auto r=reservoirs();
    auto bad=p;bad.knots[2].temperature_K=400;REQUIRE_THROWS(evaluatePhononTunneling(bad,r,c));
    bad=p;bad.knots.back().conduction_J=.1*q;REQUIRE_THROWS(evaluatePhononTunneling(bad,r,c));
    bad=p;bad.knots[1].valence_J=.1*q;REQUIRE_THROWS(evaluatePhononTunneling(bad,r,c));
    bad=p;bad.knots[1].conduction_J=.01*q;bad.knots[2].valence_J=-.01*q;
    REQUIRE_THROWS(evaluatePhononTunneling(bad,r,c)); // Multiple Keldysh crossovers.
    auto limit=c;limit.maximumIntegrationIntervals=2;REQUIRE_THROWS(evaluatePhononTunneling(p,r,limit));
    limit=c;limit.relativeIntegrationTolerance=1e-30;REQUIRE_THROWS(evaluatePhononTunneling(p,r,limit));
    auto reservoir=r;reservoir.valenceAvailableEnergy_J=-q;REQUIRE_THROWS(evaluatePhononTunneling(p,reservoir,c));
    bad=p;bad.knots[2].distance_m=bad.knots[1].distance_m;REQUIRE_THROWS(evaluatePhononTunneling(bad,r,c));
}
TEST_CASE("NLP N3 actual-state adapter rejects corners and stale paths", "[nonlocalpath][phonon]") {
    const auto c=config();const Real gap=c.ab->referenceGap_J,field=1e8,L=gap/(q*field);
    TunnelingGeometryMesh mesh;
    mesh.positions_m={Point2(0,0),Point2(3*L,0),Point2(0,3*L)};
    mesh.cells={{{0,1,2},0,1,true}};
    mesh.boundaries={{{0,1}},{{1,2}},{{0,2}}};
    TunnelingBandState s;s.nodes.resize(3);
    for(Index i=0;i<3;++i){auto& n=s.nodes[i];n.region=1;n.valence_J=-q*field*mesh.positions_m[i].x();
        n.conduction_J=n.valence_J+gap;n.conductionDOS_per_m3=n.valenceDOS_per_m3=1e25;}
    auto path=NonlocalPathGeometry(mesh).traceFromPoint(s,c,0,Point2(.1*L,.1*L));
    auto p=tunnelingProfileFromGeometry(path,s,c);auto r=saturated();
    REQUIRE(evaluatePhononTunneling(p,r,c).generation_per_m3_s>0);
    path.path.termination=TunnelingTermination::AmbiguousCorner;REQUIRE_THROWS(tunnelingProfileFromGeometry(path,s,c));
    path.path.termination=TunnelingTermination::Complete;++path.path.identity.trial;REQUIRE_THROWS(tunnelingProfileFromGeometry(path,s,c));
}

TEST_CASE("NLP independent unequal-slope triangular barriers preserve action under partition",
          "[nonlocalpath][phonon][independent_reference]") {
    // Integrate the two triangular legs directly. No A/B conversion or
    // production root/primitive is used by this reference.
    const auto c=microscopic();
    const long double mv=c.microscopic->valenceMass_kg,mc=c.microscopic->conductionMass_kg;
    const long double hb=6.626176e-34L/(2*acos(-1.L));
    for(double ae:{.7,1.2}) for(double be:{.9,1.4}) {
        const long double A=ae*q,B=be*q,L=17e-9L,u=mc*B/(mv*A+mc*B);
        const long double action=2*L/(3*hb)*
            (sqrt(2*mv*A)*pow(u,1.5L)+sqrt(2*mc*B)*pow(1-u,1.5L));
        const long double iv=2*L*hb*sqrt(u)/sqrt(2*mv*A);
        const long double ic=2*L*hb*sqrt(1-u)/sqrt(2*mc*B);
        for(int count:{1,7,31}) {
            TunnelingBandProfile p;p.initialValenceGradient_J_per_m=A/L;
            // Uneven partitions also straddle the crossover.
            for(int i=0;i<=count;++i) {
                const double z=std::pow(double(i)/count,1.7);
                p.knots.push_back({double(L)*z,-double(A)*z,double(B)*(1-z),300});
            }
            const auto a=evaluatePhononTunneling(p,saturated(),c);
            REQUIRE(a.action==Approx(double(action)).epsilon(1e-12));
            REQUIRE(a.crossover_m==Approx(double(L*u)).epsilon(1e-12).margin(0));
            REQUIRE(a.valenceReciprocalIntegral_m2==Approx(double(iv)).epsilon(1e-12).margin(0));
            REQUIRE(a.conductionReciprocalIntegral_m2==Approx(double(ic)).epsilon(1e-12).margin(0));
        }
    }
}

TEST_CASE("NLP fixed-path window and occupation directional derivatives compose",
          "[nonlocalpath][phonon][independent_reference]") {
    const auto c=microscopic();const auto p=uniform(8e7);
    TunnelingReservoirs r{-.01*q,.02*q,.02*q,.03*q};
    const auto a=evaluatePhononTunneling(p,r,c);
    const long double hb=6.626176e-34L/(2*acos(-1.L)),kt=kb*300.L;
    const long double cv=2*c.microscopic->valenceMass_kg*a.valenceReciprocalIntegral_m2/(hb*hb);
    const long double cc=2*c.microscopic->conductionMass_kg*a.conductionReciprocalIntegral_m2/(hb*hb);
    auto occ=[&](long double f){return 1/(1+exp(-f/kt));};
    const long double fp=occ(r.holeFermi_J),fn=occ(r.electronFermi_J);
    for(int mask:{1,2,3,4,8,15}) {
        const double dn=(mask&1)?.001*q:0,dp=(mask&2)?-.0007*q:0;
        const double dv=(mask&4)?.0001*q:0,dc=(mask&8)?-.0002*q:0;
        const long double tangent=(fp*(1-fp)*dp-fn*(1-fn)*dn)/(kt*(fp-fn))
            +cv*dv/expm1(cv*r.valenceAvailableEnergy_J)+cc*dc/expm1(cc*r.conductionAvailableEnergy_J);
        for(double h:{.01,.003,.001}) {
            auto difference=[&](double step) {
                auto plus=r,minus=r;
                plus.electronFermi_J+=step*dn;minus.electronFermi_J-=step*dn;
                plus.holeFermi_J+=step*dp;minus.holeFermi_J-=step*dp;
                plus.valenceAvailableEnergy_J+=step*dv;minus.valenceAvailableEnergy_J-=step*dv;
                plus.conductionAvailableEnergy_J+=step*dc;minus.conductionAvailableEnergy_J-=step*dc;
                const auto ap=evaluatePhononTunneling(p,plus,c),am=evaluatePhononTunneling(p,minus,c);
                return (ap.logAbsoluteGeneration-am.logAbsoluteGeneration)/(2*step);
            };
            // Cancel central-difference truncation rather than loosen the
            // derivative gate when the coarsest control is nonlinear.
            const double fd=(4*difference(h/2)-difference(h))/3;
            REQUIRE(fd==Approx(double(tangent)).epsilon(1e-8).margin(0));
        }
    }
}

TEST_CASE("NLP closed transverse window has a finite one-sided source derivative",
          "[nonlocalpath][phonon][independent_reference]") {
    const auto c=microscopic();const auto p=uniform();const auto r=saturated();
    const auto full=evaluatePhononTunneling(p,r,c);
    const long double hb=6.626176e-34L/(2*acos(-1.L));
    for(bool valence:{true,false}) {
        auto closed=r;
        (valence?closed.valenceAvailableEnergy_J:closed.conductionAvailableEnergy_J)=0;
        REQUIRE(evaluatePhononTunneling(p,closed,c).generation_per_m3_s==0);
        const long double coefficient=2*(valence?c.microscopic->valenceMass_kg:c.microscopic->conductionMass_kg)*
            (valence?full.valenceReciprocalIntegral_m2:full.conductionReciprocalIntegral_m2)/(hb*hb);
        // Set the dimensionless window exponent explicitly; never step into W<0.
        for(double exponent:{1e-9,3e-10,1e-10}) {
            auto open=closed;const double width=exponent/coefficient;
            (valence?open.valenceAvailableEnergy_J:open.conductionAvailableEnergy_J)=width;
            const auto rate=evaluatePhononTunneling(p,open,c).generation_per_m3_s;
            REQUIRE(rate/(width*full.generation_per_m3_s*double(coefficient))==Approx(1).epsilon(1e-8));
        }
    }
}

TEST_CASE("NLP fixed-geometry band scaling differentiates action spectrum and windows",
          "[nonlocalpath][phonon][independent_reference]") {
    const auto c=microscopic();const auto p=uniform(8e7,300,7);
    const TunnelingReservoirs r{-.01*q,.02*q,.02*q,.03*q};
    const auto base=evaluatePhononTunneling(p,r,c);
    const long double hb=6.626176e-34L/(2*acos(-1.L));
    const long double sv=2*c.microscopic->valenceMass_kg*
        base.valenceReciprocalIntegral_m2*r.valenceAvailableEnergy_J/(hb*hb);
    const long double sc=2*c.microscopic->conductionMass_kg*
        base.conductionReciprocalIntegral_m2*r.conductionAvailableEnergy_J/(hb*hb);
    // At fixed geometry: action ~ sqrt(lambda), reciprocal integrals ~
    // lambda^-1/2, and the gradient-inclusive prefactor ~ lambda^3/4.
    // Ev(start)=0 keeps the occupation arguments fixed under this scaling.
    const long double expected=.75L-base.action-.5L*(sv/expm1(sv)+sc/expm1(sc));
    auto scaled=[&](double t) {
        auto profile=p;const double factor=std::exp(t);
        profile.initialValenceGradient_J_per_m*=factor;
        for(auto& knot:profile.knots){knot.valence_J*=factor;knot.conduction_J*=factor;}
        return evaluatePhononTunneling(profile,r,c);
    };
    for(double h:{.01,.003,.001}) {
        auto difference=[&](double step) {
            const auto plus=scaled(step),minus=scaled(-step);
            return (plus.logAbsoluteGeneration-minus.logAbsoluteGeneration)/(2*step);
        };
        REQUIRE((4*difference(h/2)-difference(h))/3==
                Approx(double(expected)).epsilon(1e-8).margin(0));
        const auto plus=scaled(h);
        REQUIRE(plus.action/base.action==Approx(std::exp(h/2)).epsilon(1e-12));
        REQUIRE(plus.valenceReciprocalIntegral_m2/base.valenceReciprocalIntegral_m2==
                Approx(std::exp(-h/2)).epsilon(1e-12));
        REQUIRE(plus.conductionReciprocalIntegral_m2/base.conductionReciprocalIntegral_m2==
                Approx(std::exp(-h/2)).epsilon(1e-12));
    }
}

TEST_CASE("Airy spectral reference has positive integral and an independent derivative",
          "[nonlocalpath][phonon][independent_reference]") {
    auto spectrum=[](long double y) {
        const auto a=boost::math::airy_ai(y),b=boost::math::airy_ai_prime(y);
        return b*b-y*a*a;
    };
    for(long double y:{0.L,.5L,2.L,8.L}) {
        const auto a=boost::math::airy_ai(y);
        const auto integral=boost::math::quadrature::gauss_kronrod<long double,61>::integrate(
            [](long double t){auto z=boost::math::airy_ai(t);return z*z;},y,40.L,10,1e-14L);
        REQUIRE(spectrum(y)>0);
        REQUIRE(double(spectrum(y)/integral)==Approx(1).epsilon(1e-11));
        const long double h=1e-3L;
        auto fd=[&](long double step){return (spectrum(y+step)-spectrum(y-step))/(2*step);};
        const auto richardson=(4*fd(h/2)-fd(h))/3;
        REQUIRE(double(richardson/(-a*a))==Approx(1).epsilon(1e-8));
    }
    // The leading asymptotic form is not exact; verify its overlap improves.
    long double previous=1;
    for(long double y:{8.L,16.L,32.L}) {
        const long double asym=exp(-4*pow(y,1.5L)/3)/(8*acos(-1.L)*y);
        const long double error=abs(asym/spectrum(y)-1);
        REQUIRE(error<previous);previous=error;
        REQUIRE(error<.1L);
    }
}
