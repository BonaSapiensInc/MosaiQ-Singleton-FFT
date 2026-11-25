#pragma once

#include <cstddef>

namespace mosaiq::fft::detail {

/// Supported Singleton radices through the 61-smooth prime radix domain.
[[nodiscard]] constexpr bool is_supported_radix(std::size_t radix) noexcept {
    switch (radix) {
        case 2:
        case 3:
        case 4:
        case 5:
        case 7:
        case 11:
        case 13:
        case 17:
        case 19:
        case 23:
        case 29:
        case 31:
        case 37:
        case 41:
        case 43:
        case 47:
        case 53:
        case 59:
        case 61:
            return true;
        default:
            return false;
    }
}

inline constexpr std::size_t k_max_prime_radix = 61U;

}  // namespace mosaiq::fft::detail
