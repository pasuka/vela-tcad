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
