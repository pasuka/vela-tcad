#pragma once

namespace vela {

/**
 * @brief Bernoulli function B(x) = x / (exp(x) - 1).
 *
 * Piecewise implementation to handle numerical edge cases:
 *   |x| < 1e-10 : Taylor expansion  1 - x/2 + x^2/12
 *   x > 500     : x * exp(-x)   (avoids overflow of exp(x))
 *   x < -500    : -x            (exp(x) ~= 0 -> B(x) ~= -x)
 *   otherwise   : x / expm1(x)  (accurate for intermediate x)
 *
 * B(x) > 0 for all finite x.
 */
double bernoulli(double x);

} // namespace vela
