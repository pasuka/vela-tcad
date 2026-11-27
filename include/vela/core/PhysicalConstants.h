#pragma once

#include <numbers>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

namespace vela::constants {

// A binary has one immutable constants contract. Never combine it with material
// scalars generated under another profile. Fresh CMake builds default to the
// Sentaurus profile; an explicitly cached CODATA build remains reproducible.
#if defined(VELA_CONSTANTS_SENT_U2022_12) && VELA_CONSTANTS_SENT_U2022_12
inline constexpr std::string_view profileName = "sentaurus_u2022_12";
inline constexpr double q = 1.602192e-19;       // C
inline constexpr double kb = 1.380662e-23;      // J/K
inline constexpr double eps0 = 8.8542e-12;      // F/m (manual: 8.8542e-14 F/cm)
inline constexpr double h = 6.626176e-34;       // J*s
inline constexpr double m0 = 9.109534e-31;      // kg
// The documented DOS evaluation convention, not a geometry/trigonometry pi.
inline constexpr double dosPi = 3.141592654;
#else
inline constexpr std::string_view profileName = "codata2018";
inline constexpr double q = 1.602176634e-19;
inline constexpr double kb = 1.380649e-23;
inline constexpr double eps0 = 8.8541878128e-12;
inline constexpr double h = 6.62607015e-34;
inline constexpr double m0 = 9.1093837015e-31;
inline constexpr double dosPi = std::numbers::pi;
#endif
inline constexpr double T0 = 300.0;
inline constexpr double Vt_300 = kb * T0 / q;

constexpr long double extendedVacuumPermittivity()
{
#if defined(VELA_CONSTANTS_SENT_U2022_12) && VELA_CONSTANTS_SENT_U2022_12
    return 8.8542e-12L;
#else
    return 8.8541878128e-12L;
#endif
}

template <typename T>
struct ArithmeticConstants {
    T q;
    T kb;
    T h;
    T m0;
    T dosPi;
};

template <typename T>
constexpr ArithmeticConstants<T> arithmeticConstants()
{
    if constexpr (std::is_same_v<T, long double>) {
#if defined(VELA_CONSTANTS_SENT_U2022_12) && VELA_CONSTANTS_SENT_U2022_12
        return {1.602192e-19L, 1.380662e-23L, 6.626176e-34L, 9.109534e-31L, 3.141592654L};
#else
        return {1.602176634e-19L, 1.380649e-23L, 6.62607015e-34L, 9.1093837015e-31L,
                std::numbers::pi_v<long double>};
#endif
    } else {
        return {q, kb, h, m0, dosPi};
    }
}

template <typename T>
constexpr T thermalVoltage(T temperature)
{
    constexpr auto constants = arithmeticConstants<T>();
    return constants.kb * temperature / constants.q;
}

inline void requireProfile(std::string_view requested)
{
    if (requested != profileName)
        throw std::invalid_argument("Physical constants profile mismatch: input '" +
            std::string(requested) + "', binary '" + std::string(profileName) + "'.");
}

} // namespace vela::constants
