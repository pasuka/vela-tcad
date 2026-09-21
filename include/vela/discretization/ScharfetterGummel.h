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
 * Uses the equal-ni limit of the factorized variable-ni law, retaining
 * expm1(delta quasi-Fermi/Vt) even for sub-femtovolt conductive increments.
 * Endpoint carrier exponents use the production clamp; separated potential
 * exponentials are never materialized. Both edge potentials are explicit.
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
 * This generalizes the quasi-Fermi form to BGN/effective-ni edges. Production
 * evaluation uses a factorized expm1(delta quasi-Fermi/Vt) form, avoiding
 * cancellation of two nearly equal density flux terms. It cancels exactly for
 * flat electron quasi-Fermi potential even when ni0 != ni1.
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
 * @brief Read-only term decomposition of the production variable-ni electron SG flux.
 *
 * All values retain the caller's native solver units. The reconstructed flux
 * records the legacy subtractive form; stableFactorizedFlux mirrors the
 * production endpoint exponent clamp [-500, 500] and exact flat-quasi-Fermi
 * short circuit.
 *
 * cancellationCondition is (abs(leftTerm) + abs(rightTerm)) /
 * abs(signedDifference). Exact cancellation with nonzero terms is represented
 * by max<Real>, a finite saturation sentinel; two zero terms report 0.
 */
struct SgElectronVariableNiFluxDecomposition {
    Real ni0 = 0.0;
    Real ni1 = 0.0;
    Real n0 = 0.0;
    Real n1 = 0.0;
    Real psi0 = 0.0;
    Real psi1 = 0.0;
    Real phin0 = 0.0;
    Real phin1 = 0.0;
    Real eta = 0.0;
    Real bernoulliMinusEta = 0.0;
    Real bernoulliEta = 0.0;
    Real coef = 0.0;
    Real leftTerm = 0.0;
    Real rightTerm = 0.0;
    Real signedDifference = 0.0;
    Real reconstructedFlux = 0.0;
    Real stableFactorizedFlux = 0.0;
    Real highPrecisionReferenceFlux = 0.0;
    Real highPrecisionReferenceTermScale = 0.0;

    Real cancellationCondition = 0.0;
    bool node0ExponentClampedLow = false;
    bool node0ExponentClampedHigh = false;
    bool node1ExponentClampedLow = false;
    bool node1ExponentClampedHigh = false;
    bool includeNiGradientDrift = true;
    bool flatQuasiFermiShortCircuit = false;
};

SgElectronVariableNiFluxDecomposition
sgElectronContinuityFluxFromQuasiFermiVariableNiDecomposition(
    Real ni0,
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
 * This generalizes the quasi-Fermi form to BGN/effective-ni edges. Production
 * evaluation uses a factorized expm1(delta quasi-Fermi/Vt) form, avoiding
 * cancellation of two nearly equal density flux terms. It cancels exactly for
 * flat hole quasi-Fermi potential even when ni0 != ni1.
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
 * @brief Generalized SG flux for Fermi-Dirac electrons.
 *
 * Uses the Bessemoulin-Chatard modified thermal voltage
 *   Vt* = Vt * (eta1-eta0) / (log(n1)-log(n0)).
 * The driftPotential includes electrostatic, band-edge, and effective-ni/DOS
 * changes.  The construction reduces to classical SG in the Boltzmann limit
 * and preserves a flat quasi-Fermi level exactly.
 */
Real sgGeneralizedEinsteinFactor(
    Real density0, Real density1, Real eta0, Real eta1);

Real sgElectronFermiDiracContinuityFlux(
    Real n0, Real n1, Real eta0, Real eta1, Real driftPotential,
    Real quasiFermi0, Real quasiFermi1, Real Vt, Real coef);

/// Fermi-Dirac electron flux when the transported density has an additional
/// multiplicative quantum factor. The generalized Einstein factor is formed
/// from the underlying Fermi densities, while qmDensity enters the flux.
Real sgElectronFermiDiracQuantumContinuityFlux(
    Real qmDensity0, Real qmDensity1,
    Real fermiDensity0, Real fermiDensity1,
    Real eta0, Real eta1, Real driftPotential,
    Real quasiFermi0, Real quasiFermi1, Real Vt, Real coef);

/// Hole counterpart of sgElectronFermiDiracContinuityFlux.
Real sgHoleFermiDiracContinuityFlux(
    Real p0, Real p1, Real eta0, Real eta1, Real driftPotential,
    Real quasiFermi0, Real quasiFermi1, Real Vt, Real coef);

/// Extended-precision references for diagnostics of sub-ULP quasi-Fermi
/// drops.  The drop is supplied directly so callers can form it from
/// increment/reference components without first reconstructing an absolute
/// quasi-Fermi potential in double precision.
long double sgElectronFermiDiracQuantumContinuityFluxLongDouble(
    Real qmDensity1, Real fermiDensity0, Real fermiDensity1,
    Real eta0, Real eta1, Real driftPotential,
    long double quasiFermiDrop, Real Vt, Real coef);
long double sgHoleFermiDiracContinuityFluxLongDouble(
    Real p0, Real p1, Real eta0, Real eta1, Real driftPotential,
    long double quasiFermiDrop, Real Vt, Real coef);

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
