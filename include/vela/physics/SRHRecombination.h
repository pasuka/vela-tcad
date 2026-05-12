#pragma once

namespace vela {

/**
 * @brief Shockley-Read-Hall (SRH) net recombination rate.
 *
 * R_SRH = (n*p - ni^2) / (tau_p*(n + n1) + tau_n*(p + p1))
 *
 * Initial version: n1 = p1 = ni (trap level at mid-gap).
 *
 * A positive value indicates net recombination; negative indicates
 * net generation.
 *
 * @param n     Electron concentration [m^-3]
 * @param p     Hole concentration [m^-3]
 * @param ni    Intrinsic concentration [m^-3]
 * @param taun  Electron SRH lifetime [s]
 * @param taup  Hole SRH lifetime [s]
 * @return      Net recombination rate R [m^-3 s^-1]
 */
double srhRate(double n, double p, double ni, double taun, double taup);

} // namespace vela
