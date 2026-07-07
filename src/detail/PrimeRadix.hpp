#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <numbers>
#include <span>

namespace mosaiq::fft::detail {

template <std::floating_point T>
[[nodiscard]] constexpr T two_pi() noexcept {
    return static_cast<T>(2) * std::numbers::pi_v<T>;
}

template <std::floating_point T>
[[nodiscard]] std::complex<T> twiddle(std::size_t j, std::size_t k, std::size_t radix, int sign) noexcept {
    const T angle = static_cast<T>(sign) * two_pi<T>() * static_cast<T>(j * k) / static_cast<T>(radix);
    return {std::cos(angle), std::sin(angle)};
}

// ---------------------------------------------------------------------------
// Prime Radix Butterfly Expansion (Lemma: explicit unroll through radix 64)
// Supported radices: 2, 3, 4, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61
// ---------------------------------------------------------------------------

template <std::floating_point T, std::size_t P>
struct PrimeRadix {
    static_assert(P >= 2U, "Prime radix must be at least 2");

    static void apply(std::span<std::complex<T>, P> x, int sign) noexcept {
        std::array<std::complex<T>, P> out{};
        for (std::size_t k = 0U; k < P; ++k) {
            for (std::size_t j = 0U; j < P; ++j) {
                out[k] += x[j] * twiddle<T>(j, k, P, sign);
            }
        }
        for (std::size_t k = 0U; k < P; ++k) {
            x[k] = out[k];
        }
    }
};

template <std::floating_point T>
struct PrimeRadix<T, 2> {
    static void apply(std::span<std::complex<T>, 2U> x, int /*sign*/) noexcept {
        const std::complex<T> t0 = x[0];
        const std::complex<T> t1 = x[1];
        x[0] = t0 + t1;
        x[1] = t0 - t1;
    }
};

template <std::floating_point T>
void dispatch_prime_radix(std::span<std::complex<T>> buffer, std::size_t radix, int sign) noexcept {
    switch (radix) {
        case 2:  PrimeRadix<T, 2>::apply(std::span<std::complex<T>, 2U>{buffer.data(), 2U}, sign); break;
        case 3:  PrimeRadix<T, 3>::apply(std::span<std::complex<T>, 3U>{buffer.data(), 3U}, sign); break;
        case 4:  PrimeRadix<T, 4>::apply(std::span<std::complex<T>, 4U>{buffer.data(), 4U}, sign); break;
        case 5:  PrimeRadix<T, 5>::apply(std::span<std::complex<T>, 5U>{buffer.data(), 5U}, sign); break;
        case 7:  PrimeRadix<T, 7>::apply(std::span<std::complex<T>, 7U>{buffer.data(), 7U}, sign); break;
        case 11: PrimeRadix<T, 11>::apply(std::span<std::complex<T>, 11U>{buffer.data(), 11U}, sign); break;
        case 13: PrimeRadix<T, 13>::apply(std::span<std::complex<T>, 13U>{buffer.data(), 13U}, sign); break;
        case 17: PrimeRadix<T, 17>::apply(std::span<std::complex<T>, 17U>{buffer.data(), 17U}, sign); break;
        case 19: PrimeRadix<T, 19>::apply(std::span<std::complex<T>, 19U>{buffer.data(), 19U}, sign); break;
        case 23: PrimeRadix<T, 23>::apply(std::span<std::complex<T>, 23U>{buffer.data(), 23U}, sign); break;
        case 29: PrimeRadix<T, 29>::apply(std::span<std::complex<T>, 29U>{buffer.data(), 29U}, sign); break;
        case 31: PrimeRadix<T, 31>::apply(std::span<std::complex<T>, 31U>{buffer.data(), 31U}, sign); break;
        case 37: PrimeRadix<T, 37>::apply(std::span<std::complex<T>, 37U>{buffer.data(), 37U}, sign); break;
        case 41: PrimeRadix<T, 41>::apply(std::span<std::complex<T>, 41U>{buffer.data(), 41U}, sign); break;
        case 43: PrimeRadix<T, 43>::apply(std::span<std::complex<T>, 43U>{buffer.data(), 43U}, sign); break;
        case 47: PrimeRadix<T, 47>::apply(std::span<std::complex<T>, 47U>{buffer.data(), 47U}, sign); break;
        case 53: PrimeRadix<T, 53>::apply(std::span<std::complex<T>, 53U>{buffer.data(), 53U}, sign); break;
        case 59: PrimeRadix<T, 59>::apply(std::span<std::complex<T>, 59U>{buffer.data(), 59U}, sign); break;
        case 61: PrimeRadix<T, 61>::apply(std::span<std::complex<T>, 61U>{buffer.data(), 61U}, sign); break;
        default: break;
    }
}

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
