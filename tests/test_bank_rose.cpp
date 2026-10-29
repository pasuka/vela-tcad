#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "vela/numerics/BankRose.h"
#include "vela/core/PopulationDomainError.h"
#include <cmath>
#include <limits>

using namespace vela;
namespace {
VectorXd one(Real value) { return VectorXd::Constant(1, value); }
VectorXd squareResidual(const VectorXd& x) { return one(x[0]*x[0]-1); }
}

TEST_CASE("BankRose recovers an overshooting Newton step and approaches full steps", "[bank_rose]")
{
    BankRoseConfig c; c.notDamped = 0;
    BankRoseDamping br(c);
    auto x = one(.1);
    bool rejectedFull = false, returnedToFull = false;
    for (int i=1; i<=16; ++i) {
        const auto r = squareResidual(x);
        if (r.norm() < 1e-12) break;
        const auto step = one(-r[0]/(2*x[0]));
        auto result = br.search(i, x, step, r, squareResidual, {}, {}, {}, true);
        REQUIRE(result.accepted);
        REQUIRE(result.residualNorm <= (1-c.sufficientDecrease*result.damping)*r.norm());
        if (i == 1) {
            REQUIRE(result.attempts > 1);
            REQUIRE_FALSE(result.history.front().accepted);
            REQUIRE(result.history.front().damping == 1);
            REQUIRE(result.history[1].bankRoseK == 1);
            rejectedFull = true;
        }
        const auto& last = result.history.back();
        REQUIRE(result.damping == Catch::Approx(1/(1+last.bankRoseK*r.norm())));
        REQUIRE(br.nextK() == Catch::Approx(last.bankRoseK*c.decay));
        if (i > 1 && result.damping > .999) returnedToFull = true;
        x = result.x;
    }
    REQUIRE(rejectedFull);
    REQUIRE(returnedToFull);
    REQUIRE(x[0] == Catch::Approx(1).margin(1e-12));
}

TEST_CASE("BankRose NotDamped activates at iteration N plus one", "[bank_rose]")
{
    BankRoseConfig c; c.notDamped = 25;
    BankRoseDamping br(c);
    const auto x=one(.1), r=squareResidual(x), d=one(4.95);
    const auto full = br.search(25,x,d,r,squareResidual,{}, {}, {}, true);
    REQUIRE(full.accepted);
    REQUIRE(full.damping == 1);
    REQUIRE(full.residualNorm > r.norm());
    REQUIRE_FALSE(full.history[0].bankRoseActive);
    REQUIRE(br.nextK() == 0);
    const auto damped = br.search(26,x,d,r,squareResidual,{}, {}, {}, true);
    REQUIRE(damped.accepted);
    REQUIRE(damped.damping < 1);
    REQUIRE(damped.history.back().bankRoseActive);
}

TEST_CASE("BankRose cannot accept invalid or nonfinite trials even before activation", "[bank_rose]")
{
    BankRoseConfig c; c.notDamped=1; c.maxTrials=4;
    BankRoseDamping br(c);
    auto invalid = [](const VectorXd&,const VectorXd&){return false;};
    auto result = br.search(1,one(2),one(-1),one(1),[](const VectorXd&){return one(0);},invalid);
    REQUIRE_FALSE(result.accepted);
    REQUIRE(result.attempts == 1);
    REQUIRE(result.x[0] == 2);
    result=br.search(2,one(2),one(-1),one(1),[](const VectorXd&){return one(std::numeric_limits<Real>::quiet_NaN());});
    REQUIRE_FALSE(result.accepted);
    REQUIRE(result.attempts == 4);
    REQUIRE(result.x[0] == 2);
    REQUIRE(result.residual[0] == 1);
    REQUIRE(br.nextK() == 0);
}

TEST_CASE("BankRose failure is bounded and restores state and K", "[bank_rose]")
{
    BankRoseConfig c; c.notDamped=0; c.initialK=2; c.minDamping=.2;
    BankRoseDamping br(c);
    const auto result=br.search(1,one(0),one(1),one(1),[](const VectorXd&){return one(2);}, {}, {}, {}, true);
    REQUIRE_FALSE(result.accepted);
    REQUIRE(result.failureReason == "bank_rose_min_damping");
    REQUIRE(result.attempts == 1);
    REQUIRE(result.x[0] == 0);
    REQUIRE(result.bestRejectedCandidate);
    REQUIRE(br.nextK() == 2);
    REQUIRE_FALSE(br.search(1,one(0),one(0),one(0),squareResidual).accepted);
}

TEST_CASE("BankRose residual scaling is compensated by inverse K units", "[bank_rose]")
{
    BankRoseConfig a; a.notDamped=0; a.initialK=2;
    BankRoseConfig b=a; b.initialK/=1e8; b.seedK/=1e8;
    BankRoseDamping first(a), second(b);
    auto f=first.search(1,one(.1),one(4.95),one(-.99),squareResidual, {}, {}, {}, true);
    auto g=second.search(1,one(.1),one(4.95),one(-.99e8),
        [](const VectorXd& x)->VectorXd{return 1e8*squareResidual(x);}, {}, {}, {}, true);
    REQUIRE(f.accepted); REQUIRE(g.accepted);
    REQUIRE(f.attempts == g.attempts);
    REQUIRE(f.damping == Catch::Approx(g.damping));
    REQUIRE(f.x[0] == Catch::Approx(g.x[0]));
}

TEST_CASE("BankRose evaluates actual candidate callback and strict decrease", "[bank_rose]")
{
    BankRoseConfig c; c.notDamped=0; c.maxTrials=2;
    BankRoseDamping br(c);
    // A direction projection that discards the step must not masquerade as progress.
    auto result=br.search(1,one(2),one(-1),one(1),squareResidual, {}, {},
        [](const VectorXd& x,const VectorXd&,Real){return x;},true);
    REQUIRE_FALSE(result.accepted);
    REQUIRE(result.attempts == 2);
    REQUIRE(result.x[0] == 2);
    c.sufficientDecrease = .5;
    BankRoseDamping strict(c);
    result = strict.search(1,one(0),one(1),one(1),[](const VectorXd&){return one(.9);});
    REQUIRE_FALSE(result.accepted); // Unlike legacy fallback, any decrease is not enough.
}

TEST_CASE("BankRose rejects invalid parameters and handles K overflow", "[bank_rose]")
{
    BankRoseConfig c; c.notDamped=-1; REQUIRE_THROWS_AS(BankRoseDamping(c),std::invalid_argument);
    c={}; c.sufficientDecrease=1; REQUIRE_THROWS_AS(BankRoseDamping(c),std::invalid_argument);
    c={}; c.decay=0; REQUIRE_THROWS_AS(BankRoseDamping(c),std::invalid_argument);
    c={}; c.maxTrials=0; REQUIRE_THROWS_AS(BankRoseDamping(c),std::invalid_argument);
    c={}; c.seedK=std::numeric_limits<Real>::infinity(); REQUIRE_THROWS_AS(BankRoseDamping(c),std::invalid_argument);
    c={}; c.notDamped=0; c.initialK=1e308;
    BankRoseDamping br(c);
    const auto result=br.search(1,one(0),one(1),one(1e-308),[](const VectorXd&){return one(1e-308);}, {},
        [](const VectorXd& r){return std::abs(r[0]);});
    REQUIRE_FALSE(result.accepted);
    REQUIRE(result.failureReason == "bank_rose_k_overflow");
    REQUIRE(br.nextK() == 1e308);
}

TEST_CASE("BankRose backs away from a population domain boundary without clamping", "[bank_rose]")
{
    BankRoseConfig c; c.notDamped=0;
    BankRoseDamping br(c);
    int admissibleCalls=0;
    auto residual=[](const VectorXd& x) {
        if (!(x[0]>0)) throw PopulationDomainError("nonpositive population");
        return one(std::log(x[0]));
    };
    // Exact Newton for log(x)=0 from x=10 overshoots the positive domain.
    const auto x=one(10), r=residual(x), d=one(-10*std::log(10.));
    const auto result=br.search(1,x,d,r,residual,
        [&](const VectorXd& y,const VectorXd&){++admissibleCalls; REQUIRE(y[0]>0); return true;}, {}, {}, true);
    REQUIRE(result.accepted);
    REQUIRE(result.attempts==2);
    REQUIRE(result.history[0].rejectionReason=="bank_rose_population_domain");
    REQUIRE_FALSE(result.history[0].finite);
    REQUIRE(result.history[1].bankRoseK==1);
    REQUIRE(admissibleCalls==1);
    REQUIRE(result.x[0]==Catch::Approx(x[0]+result.damping*d[0]));
    REQUIRE(result.residualNorm<r.norm());
    REQUIRE(x[0]==10);
}

TEST_CASE("BankRose population domain failure preserves bounded search and activation", "[bank_rose]")
{
    const auto invalid=[](const VectorXd&)->VectorXd {throw PopulationDomainError("population overflow");};
    BankRoseConfig c; c.notDamped=1; c.maxTrials=3;
    BankRoseDamping br(c);
    for (int iteration:{1,2}) {
        const auto result=br.search(iteration,one(2),one(-1),one(1),invalid, {}, {}, {}, true);
        REQUIRE_FALSE(result.accepted);
        REQUIRE(result.attempts==(iteration==1?1:3));
        REQUIRE(result.failureReason=="bank_rose_population_domain");
        REQUIRE(result.x[0]==2);
        REQUIRE(result.residual[0]==1);
        REQUIRE_FALSE(result.bestRejectedCandidate);
        REQUIRE(br.nextK()==0);
    }
    c.notDamped=0; c.minDamping=.2;
    BankRoseDamping floor(c);
    const auto result=floor.search(1,one(2),one(-1),one(1),invalid);
    REQUIRE(result.attempts==2);
    REQUIRE(result.failureReason=="bank_rose_min_damping");
}

TEST_CASE("BankRose recovers typed candidate errors but propagates other errors", "[bank_rose]")
{
    BankRoseConfig c; c.notDamped=0;
    BankRoseDamping br(c);
    const auto result=br.search(1,one(2),one(-1),one(1),
        [](const VectorXd& x){return one(x[0]-1);}, {}, {},
        [](const VectorXd& x,const VectorXd& d,Real alpha)->VectorXd {
            if(alpha==1) throw PopulationDomainError("inverse population mapping");
            return x+alpha*d;
        },true);
    REQUIRE(result.accepted);
    REQUIRE(result.attempts==2);
    REQUIRE(result.history[0].rejectionReason=="bank_rose_population_domain");
    REQUIRE_THROWS_AS(br.search(1,one(2),one(-1),one(1),
        [](const VectorXd&)->VectorXd {throw std::invalid_argument("invalid model configuration");}),std::invalid_argument);
    REQUIRE_THROWS_AS(br.search(1,one(2),one(-1),one(1),squareResidual, {}, {},
        [](const VectorXd&,const VectorXd&,Real)->VectorXd {throw std::runtime_error("assembly defect");}),std::runtime_error);
    REQUIRE_THROWS_AS(br.search(1,one(2),one(-1),one(1),squareResidual,
        [](const VectorXd&,const VectorXd&)->bool {throw PopulationDomainError("guard defect");}),PopulationDomainError);
}
