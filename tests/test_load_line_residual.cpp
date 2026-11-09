#include "vela/simulation/LoadLineResidual.h"
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <array>
#include <limits>

using namespace vela;

TEST_CASE("Load-line coordinates share the same near-root physical RHS",
          "[coupled_load_line][load_line_residual]")
{
    // Frozen high-voltage rejection: scaling the three large terms separately
    // changes the RHS enough to fail the original componentwise linear gate.
    constexpr Real inner=7.170377120558898, current=0.000998722869893414;
    constexpr Real target=9994.399075811647, resistance=1e7;
    constexpr Real vt=0.025851995266484913, i0=2.98741982845799e-5;
    const auto physical=LoadLineResidual::evaluate(1,inner,resistance,current,target);
    const Real oldScaled=inner/vt+(resistance*i0/vt)*(current/i0)-target/vt;
    constexpr Real denominator=6.591291849509355e-7;
    REQUIRE(std::abs(oldScaled*vt-physical.voltage_V)/denominator>1e-10);
    const Real mapped=physical.inVoltageUnits(vt);
    REQUIRE(std::abs(mapped*vt-physical.voltage_V)/denominator<=1e-10);
    REQUIRE(physical.physicalRhs()==-physical.voltage_V);
    REQUIRE(physical.physicalRhs()==-2.4305336410179734e-7);
}

TEST_CASE("Load-line residual mapping preserves sign units and source orientation",
          "[coupled_load_line][load_line_residual]")
{
    for (Real coefficient : {-1.,0.,1.,2.}) {
        const auto r=LoadLineResidual::evaluate(coefficient,3.,4.,2.,5.);
        REQUIRE(r.voltage_V==3*coefficient+3);
        for (Real scale : {1e-3,0.025851995266484913,1.,1e3}) {
            REQUIRE(r.inVoltageUnits(scale)*scale==Catch::Approx(r.voltage_V).margin(2e-15));
            REQUIRE(-r.inVoltageUnits(scale)*scale==Catch::Approx(r.physicalRhs()).margin(2e-15));
        }
    }
    REQUIRE(LoadLineResidual::evaluate(1,2,4,2,10).voltage_V==0);
}

TEST_CASE("Load-line residual rejects invalid coordinate scales",
          "[coupled_load_line][load_line_residual]")
{
    const auto r=LoadLineResidual::evaluate(1,2,4,2,10);
    for(Real scale : {0.,-1.,std::numeric_limits<Real>::infinity(),std::numeric_limits<Real>::quiet_NaN()})
        REQUIRE_THROWS_AS(r.inVoltageUnits(scale),std::invalid_argument);
}
