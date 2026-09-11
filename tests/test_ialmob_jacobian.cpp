#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/equation/IalElementMobility.h"
#include "vela/physics/CarrierStatistics.h"
#include "vela/physics/IalInterfaceGeometry.h"
#include <algorithm>
#include <cmath>

using namespace vela;

TEST_CASE("IALMob geometry separates segment distance from vertex orientation", "[mobility][ialmob][geometry]")
{
    std::vector<std::array<Real,2>> xy{{0.,0.},{1e-6,0.},{1e-6,1e-6},{.6e-6,.2e-6}};
    std::vector<IalInterfaceSegment> edges{{0,1,{.5e-6,.5e-6}},{1,2,{.5e-6,.5e-6}}};
    const auto r=buildIalInterfaceGeometry(xy,edges);
    CHECK(r[0].orientationFamily==100);
    CHECK(r[1].orientationFamily==110);
    CHECK(r[2].orientationFamily==100);
    CHECK(r[3].orientationFamily==110);
    CHECK(r[3].distance_m==Catch::Approx(.2e-6));
    CHECK(r[3].nearestInterfaceVertex==1);
    CHECK_FALSE(r[3].onInterface);
    for (auto& p:xy) for (Real& v:p) v=17.*v+2e-6;
    for (auto& e:edges) for (Real& v:e.semiconductorPoint_m) v=17.*v+2e-6;
    const auto transformed=buildIalInterfaceGeometry(xy,edges,{0.,1.,0.},{1.,0.,0.});
    for (int i=0;i<4;++i) {
        CHECK(transformed[i].orientationFamily==r[i].orientationFamily);
        CHECK(transformed[i].distance_m==Catch::Approx(17.*r[i].distance_m).margin(1e-20));
    }
    CHECK_THROWS_AS(buildIalInterfaceGeometry(xy,edges,{1.,0.,0.},{1.,0.,0.}),std::invalid_argument);
    edges.push_back(edges[0]);
    CHECK_THROWS_AS(buildIalInterfaceGeometry(xy,edges),std::invalid_argument);
}

TEST_CASE("IALMob six SI partial derivatives match independent state perturbations", "[mobility][ialmob][jacobian]")
{
    const std::array<Real IalMobilityState::*,6> fields{
        &IalMobilityState::donors_m3,&IalMobilityState::acceptors_m3,
        &IalMobilityState::electrons_m3,&IalMobilityState::holes_m3,
        &IalMobilityState::normalField_V_per_m,&IalMobilityState::interfaceDistance_m};
    for (bool electron:{false,true}) {
        auto params=IalMobility::siliconDefaults(electron);
        params.alphaSr=5e-22; params.nu=0.; params.lambdaSr=0.;
        params.lCrit=params.lCritC=1e-6;
        const IalMobility model(params,electron);
        for (Real density:{1e18,1e23,1e26})
        for (Real field:{1e4,1e6,1e8}) {
            IalMobilityState state{2e23,8e22,density,density*.23,field,2e-8};
            const auto actual=model.evaluateWithDerivatives(state);
            for (int k=0;k<6;++k) {
                INFO("carrier="<<electron<<" density="<<density<<" field="<<field<<" variable="<<k);
                const Real scale=state.*fields[k];
                for (Real fraction:{1e-4,3e-5}) {
                    const Real step=scale*fraction;
                    auto plus=state,minus=state;
                    plus.*fields[k]+=step; minus.*fields[k]-=step;
                    const Real fd=(model.evaluate(plus).mobility_m2_per_Vs-
                        model.evaluate(minus).mobility_m2_per_Vs)/(2.*step);
                    const Real normalized=scale/actual.result.mobility_m2_per_Vs;
                    CHECK(actual.derivative_SI[k]*normalized ==
                        Catch::Approx(fd*normalized).epsilon(2e-5).margin(2e-8));
                }
            }
        }
    }
}

TEST_CASE("IALMob selected zero-input branches have finite derivatives", "[mobility][ialmob][jacobian]")
{
    for (bool electron:{false,true}) {
        const IalMobility model(IalMobility::siliconDefaults(electron),electron);
        for (Real nd:{0.,1e23})
        for (Real na:{0.,1e23})
        for (Real n:{0.,1e22})
        for (Real p:{0.,1e22})
        for (Real field:{0.,1e6}) {
            const auto result=model.evaluateWithDerivatives({nd,na,n,p,field,0.});
            CHECK(result.result.mobility_m2_per_Vs>0.);
            for (Real d:result.derivative_SI) CHECK(std::isfinite(d));
        }
    }
}

namespace {
auto geometry() {
    IalElementGeometry g;
    g.coordinates_m={{{0.,0.},{1e-7,0.},{0.,1e-7}}};
    // Deliberately non-unit distance gradient and asymmetric local measures.
    g.interfaceDistance_m={0.,2e-8,1e-8};
    g.vertexMeasure_m2={1e-15,3e-15,1e-15};
    g.boundaryTangent={.6,.8};
    return g;
}
void populate(std::array<IalElementVertexState,3>& state, bool fermi) {
    constexpr Real vt=.025852;
    for (auto& s:state) {
        const Real en=(s.potential_V-s.electronQf_V-.1)/vt;
        const Real ep=(s.holeQf_V-s.potential_V+.06)/vt;
        s.electrons_m3=1e23*(fermi?fermiDiracHalf(en):std::exp(en));
        s.holes_m3=8e22*(fermi?fermiDiracHalf(ep):std::exp(ep));
        s.electronResponse_m3_per_V=1e23/vt*(fermi?fermiDiracHalfDerivative(en):std::exp(en));
        s.holeResponse_m3_per_V=8e22/vt*(fermi?fermiDiracHalfDerivative(ep):std::exp(ep));
    }
}
}

TEST_CASE("IALMob element nine potential columns match live Fermi and Boltzmann differences", "[mobility][ialmob][jacobian]")
{
    auto ep=IalMobility::siliconDefaults(true), hp=IalMobility::siliconDefaults(false);
    ep.lCrit=ep.lCritC=hp.lCrit=hp.lCritC=1e-6;
    ep.alphaSr=5e-22; ep.nu=0.;
    const IalMobility e(ep,true),h(hp,false);
    auto ep110=ep; ep110.C*=2.;
    auto hp110=hp; hp110.delta*=.7;
    const IalMobility e110(ep110,true),h110(hp110,false);
    const std::array<const IalMobility*,3> em{&e,&e110,&e},hm{&h,&h,&h110};
    const std::array<Real IalElementVertexState::*,3> variables{
        &IalElementVertexState::potential_V,&IalElementVertexState::electronQf_V,
        &IalElementVertexState::holeQf_V};
    for (bool fermi:{false,true})
    for (bool highField:{false,true})
    for (bool boundary:{false,true})
    for (bool contact:{false,true}) {
        auto g=geometry(); g.partialBoundaryLayer=boundary;g.touchesEffectiveElectrode=contact;
        IalElementMobilityOptions options; options.highField=highField;options.referenceDensity_m3=1e23;
        std::array<IalElementVertexState,3> state;
        for (int i=0;i<3;++i) {
            state[i].potential_V=.12+.017*i;
            state[i].electronQf_V=.005*i;
            state[i].holeQf_V=.04-.003*i;
            state[i].donors_m3=(i+1)*1e23;
            state[i].acceptors_m3=4e22/(i+1);
        }
        populate(state,fermi);
        const auto r=evaluateIalElementMobility(g,state,em,hm,options);
        for (int k=0;k<9;++k) {
            INFO("Fermi="<<fermi<<" HFS="<<highField<<" partial="<<boundary<<" contact="<<contact<<" column="<<k);
            for (Real step:{2e-6,5e-7}) {
                auto plus=state,minus=state;
                plus[k/3].*variables[k%3]+=step;minus[k/3].*variables[k%3]-=step;
                populate(plus,fermi);populate(minus,fermi);
                const auto a=evaluateIalElementMobility(g,plus,em,hm,options);
                const auto b=evaluateIalElementMobility(g,minus,em,hm,options);
                CHECK(r.electron.derivative[k]/r.electron.value ==
                    Catch::Approx((a.electron.value-b.electron.value)/(2.*step*r.electron.value)).epsilon(2e-5).margin(2e-7));
                CHECK(r.hole.derivative[k]/r.hole.value ==
                    Catch::Approx((a.hole.value-b.hole.value)/(2.*step*r.hole.value)).epsilon(2e-5).margin(2e-7));
            }
        }
        // Gauge invariance involves all three potentials at every vertex.
        for (const auto& mu:{r.electron,r.hole}) {
            Real sum=0.;for (Real d:mu.derivative) sum+=d;
            CHECK(std::abs(sum)/mu.value<1e-11);
        }
        CHECK(r.electron.value<=r.electronLowField.value);
        CHECK(r.hole.value<=r.holeLowField.value);
    }
}

TEST_CASE("IALMob element rejects malformed support and preserves zero-drive limit", "[mobility][ialmob][jacobian]")
{
    const IalMobility e(IalMobility::siliconDefaults(true),true), h(IalMobility::siliconDefaults(false),false);
    const std::array<const IalMobility*,3> em{&e,&e,&e},hm{&h,&h,&h};
    auto g=geometry();std::array<IalElementVertexState,3> state{};
    const auto r=evaluateIalElementMobility(g,state,em,hm);
    CHECK(r.electron.value==Catch::Approx(.1417));
    CHECK(r.hole.value==Catch::Approx(.04705));
    for (Real d:r.electron.derivative) CHECK(d==0.);
    g.vertexMeasure_m2={0.,0.,0.};
    CHECK_THROWS_AS(evaluateIalElementMobility(g,state,em,hm),std::invalid_argument);
    g=geometry();g.partialBoundaryLayer=true;g.boundaryTangent={2.,0.};
    CHECK_THROWS_AS(evaluateIalElementMobility(g,state,em,hm),std::invalid_argument);
    g=geometry();g.coordinates_m[2]=g.coordinates_m[1];
    CHECK_THROWS_AS(evaluateIalElementMobility(g,state,em,hm),std::invalid_argument);
}
