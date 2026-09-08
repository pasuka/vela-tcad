#include <catch2/catch_test_macros.hpp>
#include "simplemos_stable_merit.hpp"
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
using simplemos_stable_merit::compare;

static vela::VectorXd vec(std::initializer_list<double> xs) {
    vela::VectorXd v(xs.size()); int i=0; for(double x:xs)v(i++)=x; return v;
}
TEST_CASE("stable merit retains hidden minority descent and rejects ascent", "[simplemos][merit]") {
    const auto x=vec({1e-10,1e-25});
    const auto down=vec({1e-10,0.999e-25}),up=vec({1e-10,1.001e-25});
    REQUIRE(x.norm()==down.norm()); REQUIRE(x.norm()==up.norm());
    REQUIRE(compare(x,down).accepted); REQUIRE_FALSE(compare(x,up).accepted);
    REQUIRE_FALSE(compare(x,x).accepted);
    REQUIRE(compare(vec({0,0}),vec({0,0})).accepted);
    REQUIRE_FALSE(compare(vec({0,0}),vec({0,1e-200})).accepted);
}
TEST_CASE("exact fallback resolves cancellation and subnormal descent", "[simplemos][merit]") {
    const auto x=vec({1,1e-150,0}),y=vec({0,0.5e-150,1});
    REQUIRE(compare(x,y).exactFallback); REQUIRE(compare(x,y).accepted);
    REQUIRE_FALSE(compare(y,x).accepted);
    REQUIRE(compare(vec({1e-310}),vec({0})).accepted);
    REQUIRE(compare(vec({1e308}),vec({0.5e308})).accepted);
    REQUIRE_FALSE(compare(vec({1}),vec({std::numeric_limits<double>::infinity()})).accepted);
    REQUIRE_FALSE(compare(vec({1}),vec({std::numeric_limits<double>::quiet_NaN()})).accepted);
}
static vela::VectorXd residual(const std::string& path) {
    std::ifstream in(path); REQUIRE(in.good());std::string line;std::getline(in,line);std::vector<double> values;
    while(std::getline(in,line)) {std::istringstream row(line);std::string field;for(int i=0;i<=5;++i)std::getline(row,field,',');values.push_back(std::stod(field));}
    vela::VectorXd out(values.size());for(int i=0;i<out.size();++i)out(i)=values[i];return out;
}
TEST_CASE("26 frozen trial signs match independently computed high precision norm", "[simplemos][merit][fixture]") {
    const char* path=std::getenv("VELA_MERIT_FIXTURE_MANIFEST");REQUIRE(path!=nullptr);
    std::ifstream in(path); REQUIRE(in.good());std::string line;int count=0,descending=0;
    while(std::getline(in,line)) {
        std::istringstream row(line);std::string base,trial,expected;
        std::getline(row,base,'\t');std::getline(row,trial,'\t');std::getline(row,expected,'\t');
        const auto result=compare(residual(base),residual(trial));
        REQUIRE(result.accepted==(expected=="1"));++count;descending+=result.accepted;
    }
    REQUIRE(count==26);REQUIRE(descending==7);
}
