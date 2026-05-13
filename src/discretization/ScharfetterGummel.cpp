#include "vela/discretization/ScharfetterGummel.h"
#include "vela/discretization/Bernoulli.h"

#include <cmath>

namespace vela {

SGEdgeWeights sgEdgeWeights(Real dpsi, Real Vt)
{
    const Real u = dpsi / Vt;
    return SGEdgeWeights{bernoulli(u), bernoulli(-u)};
}

Real sgElectronContinuityFlux(Real n0, Real n1, Real dpsi, Real Vt, Real coef)
{
    const SGEdgeWeights w = sgEdgeWeights(dpsi, Vt);
    return coef * (w.b_minus * n0 - w.b_plus * n1);
}

Real sgHoleContinuityFlux(Real p0, Real p1, Real dpsi, Real Vt, Real coef)
{
    const SGEdgeWeights w = sgEdgeWeights(dpsi, Vt);
    return coef * (w.b_plus * p0 - w.b_minus * p1);
}

Real sgElectronContinuityFluxFromQuasiFermi(Real ni0,
                                            Real psi1,
                                            Real phin0,
                                            Real phin1,
                                            Real dpsi,
                                            Real Vt,
                                            Real coef)
{
    return sgElectronContinuityFluxFromQuasiFermiFactors(
        ni0,
        std::exp(psi1 / Vt),
        std::exp(-phin0 / Vt),
        std::exp(-phin1 / Vt),
        dpsi,
        Vt,
        coef);
}

Real sgElectronContinuityFluxFromQuasiFermiFactors(Real ni0,
                                                   Real expPsi1,
                                                   Real expNegPhin0,
                                                   Real expNegPhin1,
                                                   Real dpsi,
                                                   Real Vt,
                                                   Real coef)
{
    const Real Bu = bernoulli(dpsi / Vt);
    return coef * Bu * ni0 * expPsi1 * (expNegPhin0 - expNegPhin1);
}

Real sgHoleContinuityFluxFromQuasiFermi(Real ni0,
                                        Real psi0,
                                        Real phip0,
                                        Real phip1,
                                        Real dpsi,
                                        Real Vt,
                                        Real coef)
{
    return sgHoleContinuityFluxFromQuasiFermiFactors(
        ni0,
        std::exp(-psi0 / Vt),
        std::exp(phip0 / Vt),
        std::exp(phip1 / Vt),
        dpsi,
        Vt,
        coef);
}

Real sgHoleContinuityFluxFromQuasiFermiFactors(Real ni0,
                                               Real expNegPsi0,
                                               Real expPhip0,
                                               Real expPhip1,
                                               Real dpsi,
                                               Real Vt,
                                               Real coef)
{
    const Real Bu = bernoulli(dpsi / Vt);
    return coef * Bu * ni0 * expNegPsi0 * (expPhip0 - expPhip1);
}

double sgElectronFlux(double n0, double n1, double dpsi, double Vt,
                      double mu, double h)
{
    return -sgElectronContinuityFlux(n0, n1, dpsi, Vt, mu * Vt / h);
}

double sgHoleFlux(double p0, double p1, double dpsi, double Vt,
                  double mu, double h)
{
    return -sgHoleContinuityFlux(p0, p1, dpsi, Vt, mu * Vt / h);
}

} // namespace vela
