#pragma once

namespace vela {

/**
 * @brief Scharfetter-Gummel edge fluxes for drift-diffusion.
 *
 * Computes the conventional current density on a mesh edge connecting
 * nodes i (index 0) and j (index 1).  J_ij denotes current flowing
 * from node i toward node j.
 *
 * Electron flux (conventional current from i to j):
 *   Jn_ij = μn * Vt / h * ( B(+u) * n_j - B(-u) * n_i )
 *
 * Hole flux (conventional current from i to j):
 *   Jp_ij = μp * Vt / h * ( B(-u) * p_j - B(+u) * p_i )
 *
 * where  u = (ψ_j - ψ_i) / Vt  and  B is the Bernoulli function.
 *
 * Note: The factor q (elementary charge) is NOT included here; callers
 * multiply by q when assembling current densities in SI units [A/m²].
 *
 * @param n0   Carrier density at node i [m^-3]
 * @param n1   Carrier density at node j [m^-3]
 * @param dpsi Potential difference ψ_j - ψ_i [V]
 * @param Vt   Thermal voltage [V]
 * @param mu   Carrier mobility [m²/V/s]
 * @param h    Edge length [m]
 * @return     Current density from i to j [m^-2 s^-1] (multiply by q for A/m²)
 */
double sgElectronFlux(double n0, double n1, double dpsi, double Vt,
                      double mu, double h);

double sgHoleFlux    (double p0, double p1, double dpsi, double Vt,
                      double mu, double h);

} // namespace vela
