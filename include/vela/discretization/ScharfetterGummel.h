#pragma once

#include "vela/core/Types.h"

namespace vela {

/**
 * @brief Scharfetter-Gummel Bernoulli weights for an oriented edge 0 -> 1.
 *
 * For u = (psi_1 - psi_0) / Vt, b_plus = B(+u) and b_minus = B(-u).
 */
struct SGEdgeWeights {
    Real b_plus = 0.0;
    Real b_minus = 0.0;
};

SGEdgeWeights sgEdgeWeights(Real dpsi, Real Vt);

/**
 * @brief Electron continuity flux contribution from node 0 to node 1.
 *
 * This is the sign convention used by DDAssembler and CoupledDDAssembler
 * continuity residual rows:
 *   F_n01 = coef * (B(-u) * n0 - B(+u) * n1)
 * where coef contains mobility, thermal voltage, and edge geometry.
 */
Real sgElectronContinuityFlux(Real n0, Real n1, Real dpsi, Real Vt, Real coef);

/**
 * @brief Hole continuity flux contribution from node 0 to node 1.
 *
 * Continuity residual convention:
 *   F_p01 = coef * (B(+u) * p0 - B(-u) * p1)
 * where coef contains mobility, thermal voltage, and edge geometry.
 */
Real sgHoleContinuityFlux(Real p0, Real p1, Real dpsi, Real Vt, Real coef);

/**
 * @brief Balanced electron continuity flux for Boltzmann quasi-Fermi variables.
 *
 * Algebraically equivalent to sgElectronContinuityFlux when both edge nodes use
 * the same intrinsic density, but it evaluates equilibrium edges without the
 * subtractive cancellation in B(-u)*n0 - B(+u)*n1.
 */
Real sgElectronContinuityFluxFromQuasiFermi(Real ni0,
                                            Real psi1,
                                            Real phin0,
                                            Real phin1,
                                            Real dpsi,
                                            Real Vt,
                                            Real coef);

/**
 * @brief Cancellation-robust balanced electron continuity flux.
 *
 * Numerically equivalent to sgElectronContinuityFluxFromQuasiFermi but factors
 * out the larger carrier-density exponential, so the weighted subtraction stays
 * at O(1) magnitude. This avoids catastrophic cancellation when (psi - phin)/Vt
 * is large (heavy band bending) and the exp(psi/Vt) overflow of the separated
 * factor form when |psi| is large. Both edge potentials are passed explicitly.
 */
Real sgElectronContinuityFluxFromQuasiFermiStable(Real ni0,
                                                  Real psi0,
                                                  Real psi1,
                                                  Real phin0,
                                                  Real phin1,
                                                  Real Vt,
                                                  Real coef);

/**
 * @brief Balanced electron continuity flux with precomputed Boltzmann factors.
 *
 * Uses expPsi1 = exp(psi1/Vt) and expNegPhin0/1 = exp(-phin0/1 / Vt), allowing
 * callers to cache the per-node exponentials across all incident edges.
 */
Real sgElectronContinuityFluxFromQuasiFermiFactors(Real ni0,
                                                   Real expPsi1,
                                                   Real expNegPhin0,
                                                   Real expNegPhin1,
                                                   Real dpsi,
                                                   Real Vt,
                                                   Real coef);

/**
 * @brief Balanced electron continuity flux for edge-varying intrinsic density.
 *
 * This generalizes the quasi-Fermi form to BGN/effective-ni edges. It cancels
 * exactly for flat electron quasi-Fermi potential even when ni0 != ni1.
 *
 * When includeNiGradientDrift is false the intrinsic-density gradient term
 * log(ni1/ni0) is dropped from the Scharfetter-Gummel argument, reducing the
 * flux to the plain density-based form. This is appropriate when the per-node
 * ni variation is a material discontinuity rather than a smooth bandgap-
 * narrowing gradient.
 */
Real sgElectronContinuityFluxFromQuasiFermiVariableNi(Real ni0,
                                                      Real ni1,
                                                      Real psi0,
                                                      Real psi1,
                                                      Real phin0,
                                                      Real phin1,
                                                      Real Vt,
                                                      Real coef,
                                                      bool includeNiGradientDrift = true);

/**
 * @brief Balanced hole continuity flux for Boltzmann quasi-Fermi variables.
 *
 * Algebraically equivalent to sgHoleContinuityFlux when both edge nodes use the
 * same intrinsic density, but it evaluates equilibrium edges without the
 * subtractive cancellation in B(+u)*p0 - B(-u)*p1.
 */
Real sgHoleContinuityFluxFromQuasiFermi(Real ni0,
                                        Real psi0,
                                        Real phip0,
                                        Real phip1,
                                        Real dpsi,
                                        Real Vt,
                                        Real coef);

/**
 * @brief Cancellation-robust balanced hole continuity flux.
 *
 * See sgElectronContinuityFluxFromQuasiFermiStable for the numerical rationale.
 * Both edge potentials are passed explicitly.
 */
Real sgHoleContinuityFluxFromQuasiFermiStable(Real ni0,
                                              Real psi0,
                                              Real psi1,
                                              Real phip0,
                                              Real phip1,
                                              Real Vt,
                                              Real coef);

/**
 * @brief Balanced hole continuity flux with precomputed Boltzmann factors.
 *
 * Uses expNegPsi0 = exp(-psi0/Vt) and expPhip0/1 = exp(phip0/1 / Vt), allowing
 * callers to cache the per-node exponentials across all incident edges.
 */
Real sgHoleContinuityFluxFromQuasiFermiFactors(Real ni0,
                                               Real expNegPsi0,
                                               Real expPhip0,
                                               Real expPhip1,
                                               Real dpsi,
                                               Real Vt,
                                               Real coef);

/**
 * @brief Balanced hole continuity flux for edge-varying intrinsic density.
 *
 * This generalizes the quasi-Fermi form to BGN/effective-ni edges. It cancels
 * exactly for flat hole quasi-Fermi potential even when ni0 != ni1.
 *
 * When includeNiGradientDrift is false the intrinsic-density gradient term
 * log(ni0/ni1) is dropped from the Scharfetter-Gummel argument, reducing the
 * flux to the plain density-based form. This is appropriate when the per-node
 * ni variation is a material discontinuity rather than a smooth bandgap-
 * narrowing gradient.
 */
Real sgHoleContinuityFluxFromQuasiFermiVariableNi(Real ni0,
                                                  Real ni1,
                                                  Real psi0,
                                                  Real psi1,
                                                  Real phip0,
                                                  Real phip1,
                                                  Real Vt,
                                                  Real coef,
                                                  bool includeNiGradientDrift = true);

/**
 * @brief Scharfetter-Gummel edge fluxes for drift-diffusion.
 *
 * Computes the conventional current density on a mesh edge connecting
 * nodes i (index 0) and j (index 1).  J_ij denotes current flowing
 * from node i toward node j.
 *
 * Electron flux (conventional current from i to j):
 *   Jn_ij = mu_n * Vt / h * ( B(+u) * n_j - B(-u) * n_i )
 *
 * Hole flux (conventional current from i to j):
 *   Jp_ij = mu_p * Vt / h * ( B(-u) * p_j - B(+u) * p_i )
 *
 * where  u = (psi_j - psi_i) / Vt  and  B is the Bernoulli function.
 *
 * Note: The factor q (elementary charge) is NOT included here; callers
 * multiply by q when assembling current densities in SI units [A/m^2].
 *
 * @param n0   Carrier density at node i [m^-3]
 * @param n1   Carrier density at node j [m^-3]
 * @param dpsi Potential difference psi_j - psi_i [V]
 * @param Vt   Thermal voltage [V]
 * @param mu   Carrier mobility [m^2/V/s]
 * @param h    Edge length [m]
 * @return     Current density from i to j [m^-2 s^-1] (multiply by q for A/m^2)
 */
double sgElectronFlux(double n0, double n1, double dpsi, double Vt,
                      double mu, double h);

double sgHoleFlux    (double p0, double p1, double dpsi, double Vt,
                      double mu, double h);

} // namespace vela
