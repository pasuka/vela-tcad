#include <catch2/catch_test_macros.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "simplemos_stable_sg_psi_derivative.hpp"
#include "vela/discretization/ScharfetterGummel.h"
using HP=boost::multiprecision::cpp_dec_float_100;

namespace {
HP B(HP x) { if(x==0)return HP(1);return x/expm1(x); }
HP flux(HP a,HP b,HP q0,HP q1,bool electron) {
    const HP vt(".02585"),ni("1e10"),coef("2.4");const HP eta=(b-a)/vt;
    if(electron)return coef*ni*(exp((a-q0)/vt)*B(-eta)-exp((b-q1)/vt)*B(eta));
    return coef*ni*(exp((q0-a)/vt)*B(eta)-exp((q1-b)/vt)*B(-eta));
}
}

TEST_CASE("stable SG psi derivatives agree with independent 100-digit density-flux differences", "[simplemos][stable_sg]") {
    for(bool electron:{false,true})for(double potential:{-.43,.13})for(double delta:{0.,1e-18,1e-10,.03}) {
        const double a=potential,b=potential+.011,vt=.02585;const double eta=(b-a)/vt;
        const double f=electron?vela::sgElectronBoltzmannContinuityFlux(1e10,1e10,a,b,0.,delta,vt,2.4,{false,true})
                               :vela::sgHoleBoltzmannContinuityFlux(1e10,1e10,a,b,0.,delta,vt,2.4,{false,true});
        const HP x=electron?-HP(eta):HP(eta);const HP em1=expm1(x);const double db=static_cast<double>((em1-x*(em1+1))/(em1*em1));
        const auto derivative=simplemos_stable_sg::psiDerivative(f,vela::bernoulli(static_cast<double>(x)),db,vt,electron);
        const HP step("1e-20");
        for(int end=0;end<2;++end) {
            const HP plus=flux(HP(a)+(end==0?step:HP(0)),HP(b)+(end==1?step:HP(0)),HP(0),HP(delta),electron);
            const HP minus=flux(HP(a)-(end==0?step:HP(0)),HP(b)-(end==1?step:HP(0)),HP(0),HP(delta),electron);
            const HP ref=(plus-minus)/(2*step);
            if(delta==0.)REQUIRE(derivative[end]==0.);
            else REQUIRE(static_cast<double>(abs((HP(derivative[end])-ref)/ref))<2e-12);
        }
    }
}
