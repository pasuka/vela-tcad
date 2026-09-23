#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/equation/IalElementAvalanche.h"
#include "vela/physics/CarrierStatistics.h"
#include <algorithm>
#include <cmath>
#include <limits>

using namespace vela;
namespace {
constexpr Real vt=.025851999786435;
using State=std::array<IalElementVertexState,3>;
struct Fixture {
    IalMobility em{IalMobility::siliconDefaults(true),true};
    IalMobility hm{IalMobility::siliconDefaults(false),false};
    std::array<const IalMobility*,3> electrons{&em,&em,&em},holes{&hm,&hm,&hm};
    IalElementGeometry g;
    IalElementAvalancheOptions options;
    State reference;
    Fixture() {
        g.coordinates_m={{{0.,0.},{1e-8,0.},{0.,1e-8}}};
        g.interfaceDistance_m={1e-9,4e-9,2e-9};
        g.vertexMeasure_m2={1.2e-17,2.3e-17,1.5e-17};
        g.boundaryTangent={.6,.8};
        options.mobilityVertex=1;
        options.sourceMeasure_m2={1e-17,1.5e-17,2.5e-17};
        // Independent low/high curves exercise the mixing-weight derivative.
        options.electron={2e7,1e6,7e7,3e7};
        for (int i=0;i<3;++i) {
            reference[i].potential_V=std::array<Real,3>{.12,.22,.16}[i];
            reference[i].electronQf_V=std::array<Real,3>{.02,.35,.25}[i];
            reference[i].holeQf_V=std::array<Real,3>{.015,.34,.23}[i];
            reference[i].donors_m3=1e23*(i+1);
            reference[i].acceptors_m3=2e22;
        }
    }
    void populate(State& state,bool fermi) const {
        for (int i=0;i<3;++i) {
            auto& s=state[i];const auto& r=reference[i];
            const Real en=.5+i+(s.potential_V-r.potential_V-s.electronQf_V+r.electronQf_V)/vt;
            const Real ep=2.5-.7*i+(s.holeQf_V-r.holeQf_V-s.potential_V+r.potential_V)/vt;
            s.electrons_m3=1e23*(fermi?fermiDiracHalf(en):std::exp(en));
            s.holes_m3=8e22*(fermi?fermiDiracHalf(ep):std::exp(ep));
            s.electronResponse_m3_per_V=1e23/vt*(fermi?fermiDiracHalfDerivative(en):std::exp(en));
            s.holeResponse_m3_per_V=8e22/vt*(fermi?fermiDiracHalfDerivative(ep):std::exp(ep));
        }
    }
    auto evaluate(const State& s) const {
        return evaluateIalElementAvalanche(g,s,electrons,holes,options);
    }
};
void perturb(State& s,int col,Real delta) {
    const std::array<Real IalElementVertexState::*,3> fields{
        &IalElementVertexState::potential_V,&IalElementVertexState::electronQf_V,&IalElementVertexState::holeQf_V};
    s[col/3].*fields[col%3]+=delta;
}
}

TEST_CASE("Element density-gradient avalanche includes all nine live mobility columns",
          "[avalanche][ialmob][element-density][jacobian]") {
    for (bool fermi:{false,true}) for (int boundary:{0,1,2,3}) {
        Fixture f;f.g.partialBoundaryLayer=(boundary&1);f.g.touchesEffectiveElectrode=(boundary&2);
        if(f.g.touchesEffectiveElectrode)
            for(auto& s:f.reference)s.potential_V=s.electronQf_V;
        auto state=f.reference;f.populate(state,fermi);const auto actual=f.evaluate(state);
        for(int col=0;col<9;++col) for(Real step:{1e-6,3e-7}) {
            INFO("fermi="<<fermi<<" boundary="<<boundary<<" column="<<col<<" step="<<step);
            auto plus=state,minus=state;perturb(plus,col,step);perturb(minus,col,-step);
            f.populate(plus,fermi);f.populate(minus,fermi);
            const auto a=f.evaluate(plus),b=f.evaluate(minus);
            for(bool hole:{false,true}) {
                const auto& r=hole?actual.hole:actual.electron;
                const auto& rp=hole?a.hole:a.electron;const auto& rm=hole?b.hole:b.electron;
                REQUIRE(r.generation_m3_per_s.value>0.);
                const Real fd=(rp.generation_m3_per_s.value-rm.generation_m3_per_s.value)/(2.*step);
                CHECK(r.generation_m3_per_s.derivative[col]/r.generation_m3_per_s.value
                    ==Catch::Approx(fd/r.generation_m3_per_s.value).epsilon(3e-5).margin(3e-7));
                for(int vertex=0;vertex<3;++vertex) {
                    const Real sourceFd=(rp.source_per_m_s[vertex].value-rm.source_per_m_s[vertex].value)/(2.*step);
                    CHECK(r.source_per_m_s[vertex].derivative[col]/r.source_per_m_s[vertex].value
                        ==Catch::Approx(sourceFd/r.source_per_m_s[vertex].value).epsilon(3e-5).margin(3e-7));
                }
            }
        }
        // A frozen-mobility Jacobian must differ materially in this fixture.
        for(const auto& r:{actual.electron,actual.hole}) {
            Real contribution=0.,gauge=0.;
            for(int k=0;k<9;++k) {
                contribution=std::max(contribution,std::abs(r.mobility_m2_per_Vs.derivative[k]/r.mobility_m2_per_Vs.value));
                gauge+=r.generation_m3_per_s.derivative[k]/r.generation_m3_per_s.value;
            }
            CHECK(contribution>1e-2);CHECK(std::abs(gauge)<1e-10);
            Real sum=0.;for(const auto& x:r.source_per_m_s)sum+=x.value;
            CHECK(sum==Catch::Approx(r.generation_m3_per_s.value*5e-17).epsilon(2e-14));
        }
    }
}

TEST_CASE("Element avalanche density uses extrema and declares a one-sided tie branch",
          "[avalanche][ialmob][element-density]") {
    Fixture f;auto state=f.reference;f.populate(state,false);
    for(int i=0;i<3;++i) {
        state[i].electrons_m3=2e22*std::pow(2.,i);
        state[i].electronResponse_m3_per_V=state[i].electrons_m3/vt;
    }
    auto r=f.evaluate(state).electron.density_m3;
    CHECK(r.value==Catch::Approx(3.2e22));
    CHECK(r.derivative[3]==0.);CHECK(r.derivative[4]==0.);
    CHECK(r.derivative[0]==Catch::Approx(1.28*state[0].electronResponse_m3_per_V));
    CHECK(r.derivative[6]==Catch::Approx(.08*state[2].electronResponse_m3_per_V));
    state[1].electrons_m3=state[0].electrons_m3;
    state[1].electronResponse_m3_per_V=state[0].electronResponse_m3_per_V;
    r=f.evaluate(state).electron.density_m3;
    // Decreasing the first tied minimum preserves the documented selected branch.
    auto minus=state;constexpr Real step=1e-8;
    minus[0].electrons_m3*=std::exp(-step/vt);
    minus[0].electronResponse_m3_per_V=minus[0].electrons_m3/vt;
    const Real fd=(r.value-f.evaluate(minus).electron.density_m3.value)/step;
    CHECK(r.derivative[0]==Catch::Approx(fd).epsilon(1e-6));
}

TEST_CASE("Element avalanche contact alpha does not create current at constant quasi Fermi potential",
          "[avalanche][ialmob][element-density]") {
    Fixture f;f.g.touchesEffectiveElectrode=true;auto state=f.reference;f.populate(state,false);
    for(auto& s:state){s.electronQf_V=.1;s.holeQf_V=.2;}
    const auto r=f.evaluate(state);
    for(const auto& c:{r.electron,r.hole}) {
        CHECK(c.alpha_per_m.value>0.);CHECK(c.currentField_V_per_m.value==0.);
        CHECK(c.generation_m3_per_s.value==0.);
        for(Real x:c.generation_m3_per_s.derivative)CHECK(x==0.);
    }
}

TEST_CASE("Element avalanche cosine curve is C1 at both transition boundaries",
          "[avalanche][ialmob][element-density]") {
    Fixture f;auto state=f.reference;f.populate(state,false);
    for(Real ratio:{.95,1.05}) {
        const Real force=ratio*f.options.electron.switchField_V_per_m;
        state[0].electronQf_V=0.;state[1].electronQf_V=force*1e-8;state[2].electronQf_V=0.;
        const auto center=f.evaluate(state).electron.alpha_per_m;
        for(Real direction:{-1.,1.}) {
            auto moved=state;moved[1].electronQf_V+=direction*1e-10;
            const auto far=f.evaluate(moved).electron.alpha_per_m;
            moved=state;moved[1].electronQf_V+=direction*1e-12;
            const auto near=f.evaluate(moved).electron.alpha_per_m;
            // C1 is a limiting property: test contraction toward each join,
            // including the large curvature of the independent test curves.
            CHECK(std::abs(near.derivative[4]-center.derivative[4])
                <=.02*std::abs(far.derivative[4]-center.derivative[4])+1e-6);
            CHECK(near.value==Catch::Approx(center.value).epsilon(1e-8));
            CHECK(near.derivative[4]==Catch::Approx(center.derivative[4]).epsilon(1e-6));
        }
    }
}

TEST_CASE("Element avalanche rejects unsupported temperature and invalid support",
          "[avalanche][ialmob][element-density]") {
    Fixture f;auto state=f.reference;f.populate(state,false);
    SECTION("temperature"){state[0].temperature_K=301.;}
    SECTION("selected vertex"){f.options.mobilityVertex=3;}
    SECTION("negative measure"){f.options.sourceMeasure_m2[1]=-1.;}
    SECTION("empty support"){f.options.sourceMeasure_m2={0.,0.,0.};}
    SECTION("bad coefficient"){f.options.hole.lowA_per_m=std::numeric_limits<Real>::quiet_NaN();}
    SECTION("bad density"){state[1].electrons_m3=-1.;}
    CHECK_THROWS_AS(f.evaluate(state),std::invalid_argument);
}
