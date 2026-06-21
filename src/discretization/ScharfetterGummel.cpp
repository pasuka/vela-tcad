#include "vela/discretization/ScharfetterGummel.h"
#include "vela/discretization/Bernoulli.h"

#include <algorithm>
#include <cmath>

namespace vela {

namespace {

Real limitedExp(Real value)
{
    return std::exp(std::clamp(value, -500.0, 500.0));
}

} // namespace

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
    if (phin0 == phin1)
        return 0.0;

    const Real psi0 = psi1 - dpsi;
    return sgElectronContinuityFlux(
        ni0 * limitedExp((psi0 - phin0) / Vt),
        ni0 * limitedExp((psi1 - phin1) / Vt),
        dpsi,
        Vt,
        coef);
}

Real sgElectronContinuityFluxFromQuasiFermiStable(Real ni0,
                                                  Real psi0,
                                                  Real psi1,
                                                  Real phin0,
                                                  Real phin1,
                                                  Real Vt,
                                                  Real coef)
{
    if (phin0 == phin1)
        return 0.0;

    // Separated-factor form of coef*(B(-u)*n0 - B(u)*n1). Evaluating the
    // quasi-Fermi difference (exp(-phin0/Vt) - exp(-phin1/Vt)) directly avoids
    // the catastrophic cancellation of subtracting two large nearly-equal
    // carrier densities when (psi - phin)/Vt is large (heavy band bending).
    const Real u = (psi1 - psi0) / Vt;
    return sgElectronContinuityFluxFromQuasiFermiFactors(
        ni0,
        limitedExp(psi1 / Vt),
        limitedExp(-phin0 / Vt),
        limitedExp(-phin1 / Vt),
        u * Vt,
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

Real sgElectronContinuityFluxFromQuasiFermiVariableNi(Real ni0,
                                                      Real ni1,
                                                      Real psi0,
                                                      Real psi1,
                                                      Real phin0,
                                                      Real phin1,
                                                      Real Vt,
                                                      Real coef,
                                                      bool includeNiGradientDrift)
{
    if (phin0 == phin1)
        return 0.0;
    if (ni0 <= 0.0 || ni1 <= 0.0)
        return sgElectronContinuityFlux(
            ni0 * limitedExp((psi0 - phin0) / Vt),
            ni1 * limitedExp((psi1 - phin1) / Vt),
            psi1 - psi0,
            Vt,
            coef);

    const Real eta = (psi1 - psi0) / Vt
        + (includeNiGradientDrift ? std::log(ni1 / ni0) : 0.0);
    const Real n0 = ni0 * limitedExp((psi0 - phin0) / Vt);
    const Real n1 = ni1 * limitedExp((psi1 - phin1) / Vt);
    return coef * (bernoulli(-eta) * n0 - bernoulli(eta) * n1);
}

Real sgHoleContinuityFluxFromQuasiFermi(Real ni0,
                                        Real psi0,
                                        Real phip0,
                                        Real phip1,
                                        Real dpsi,
                                        Real Vt,
                                        Real coef)
{
    if (phip0 == phip1)
        return 0.0;

    const Real psi1 = psi0 + dpsi;
    return sgHoleContinuityFlux(
        ni0 * limitedExp((phip0 - psi0) / Vt),
        ni0 * limitedExp((phip1 - psi1) / Vt),
        dpsi,
        Vt,
        coef);
}

Real sgHoleContinuityFluxFromQuasiFermiStable(Real ni0,
                                              Real psi0,
                                              Real psi1,
                                              Real phip0,
                                              Real phip1,
                                              Real Vt,
                                              Real coef)
{
    if (phip0 == phip1)
        return 0.0;

    // Separated-factor form of coef*(B(u)*p0 - B(-u)*p1). See the electron
    // variant above for the numerical rationale.
    const Real u = (psi1 - psi0) / Vt;
    return sgHoleContinuityFluxFromQuasiFermiFactors(
        ni0,
        limitedExp(-psi0 / Vt),
        limitedExp(phip0 / Vt),
        limitedExp(phip1 / Vt),
        u * Vt,
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

Real sgHoleContinuityFluxFromQuasiFermiVariableNi(Real ni0,
                                                  Real ni1,
                                                  Real psi0,
                                                  Real psi1,
                                                  Real phip0,
                                                  Real phip1,
                                                  Real Vt,
                                                  Real coef,
                                                  bool includeNiGradientDrift)
{
    if (phip0 == phip1)
        return 0.0;
    if (ni0 <= 0.0 || ni1 <= 0.0)
        return sgHoleContinuityFlux(
            ni0 * limitedExp((phip0 - psi0) / Vt),
            ni1 * limitedExp((phip1 - psi1) / Vt),
            psi1 - psi0,
            Vt,
            coef);

    const Real eta = (psi1 - psi0) / Vt
        + (includeNiGradientDrift ? std::log(ni0 / ni1) : 0.0);
    const Real p0 = ni0 * limitedExp((phip0 - psi0) / Vt);
    const Real p1 = ni1 * limitedExp((phip1 - psi1) / Vt);
    return coef * (bernoulli(eta) * p0 - bernoulli(-eta) * p1);
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
