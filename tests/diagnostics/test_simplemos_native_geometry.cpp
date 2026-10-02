#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include "simplemos_native_geometry.hpp"
#include "vela/discretization/ScharfetterGummel.h"
#include <filesystem>
#include <fstream>

namespace {
struct Ratios {
    std::filesystem::path path;
    Ratios(const std::string& rows) {
        path=std::filesystem::temp_directory_path()/"vela_joint_geometry_unit_ratios.txt";
        std::ofstream out(path);out<<rows;out.close();
        _putenv_s("VELA_TEST_GEOMETRY_PATH",path.string().c_str());
        _putenv_s("VELA_TEST_GEOMETRY_ALPHA","1");
    }
    ~Ratios() {
        _putenv_s("VELA_TEST_GEOMETRY_PATH","");_putenv_s("VELA_TEST_GEOMETRY_ALPHA","");
        std::filesystem::remove(path);
    }
    void apply(std::vector<double>& x) {
        simplemos_native_geometry::apply(x,"VELA_TEST_GEOMETRY_PATH","VELA_TEST_GEOMETRY_ALPHA");
    }
};
}

TEST_CASE("geometry replacement preserves semiconductor thermal equilibrium", "[simplemos][native_geometry]") {
    Ratios ratios("0 0.2\n1 0\n2 1.3\n");
    std::vector<double> coefficients{2.,7.,3.};ratios.apply(coefficients);
    for(double c:coefficients) {
        REQUIRE(vela::sgElectronContinuityFluxFromQuasiFermiStable(1e10,-.12,.17,.03,.03,.02585,c)==0.);
    }
}

TEST_CASE("modified two-edge transport obeys series conductance and local conservation", "[simplemos][native_geometry]") {
    Ratios ratios("0 0.2\n1 0.6\n");
    std::vector<double> g{2.,5.};ratios.apply(g);
    const double n0=10.,n2=2.;
    const double n1=(g[0]*n0+g[1]*n2)/(g[0]+g[1]);
    const double left=vela::sgElectronContinuityFlux(n0,n1,0.,.02585,g[0]);
    const double right=vela::sgElectronContinuityFlux(n1,n2,0.,.02585,g[1]);
    REQUIRE(left==Catch::Approx(right).epsilon(1e-14));
    REQUIRE(left==Catch::Approx((n0-n2)/(1./.4+1./3.)).epsilon(1e-14));
    REQUIRE(vela::sgElectronContinuityFlux(n1,n0,0.,.02585,g[0])==Catch::Approx(-left));
}

TEST_CASE("common Poisson charge support preserves local neutrality", "[simplemos][native_geometry]") {
    Ratios ratios("0 0.2\n1 0\n2 1.3\n");
    std::vector<double> electron{2.,7.,3.},hole=electron,dopant=electron;
    ratios.apply(electron);ratios.apply(hole);ratios.apply(dopant);
    const double n=9.,p=1.,nd=8.;
    for(std::size_t i=0;i<electron.size();++i) {
        const double charge=n*electron[i]-p*hole[i]-nd*dopant[i];
        REQUIRE(std::abs(charge)<=1e-14*std::max(1.,electron[i]));
    }
}
