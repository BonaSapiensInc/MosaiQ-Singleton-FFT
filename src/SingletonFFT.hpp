#pragma once

#include <array>
#include <complex>
#include <concepts>
#include <cstddef>
#include <span>
#include <vector>

namespace mosaiq::fft {

/// Bona Sapiens Determinism Contract: type-safe floating-point FFT scalar.
template <typename T>
concept FloatingPoint = std::floating_point<T>;

/// Prime factorization result for mixed-radix decomposition.
struct Factorization {
    std::vector<std::size_t> factors;
    std::size_t length{0U};
};

/// Compile-time factorization buffer (no heap allocation).
template <std::size_t MaxFactors = 32U>
struct ConstFactorization {
    std::array<std::size_t, MaxFactors> factors{};
    std::size_t count{0U};
    std::size_t length{0U};

    [[nodiscard]] constexpr bool operator==(const Factorization& other) const noexcept {
        if (length != other.length || count != other.factors.size()) {
            return false;
        }
        for (std::size_t index = 0U; index < count; ++index) {
            if (factors[index] != other.factors[index]) {
                return false;
            }
        }
        return true;
    }
};

/// Runtime factorizer: decomposes \f$N\f$ into Singleton-style radices.
[[nodiscard]] Factorization factorize(std::size_t n);

namespace detail {

template <std::size_t N>
[[nodiscard]] consteval ConstFactorization<> factorize_length() {
    if (N < 2U) {
        return ConstFactorization<>{};
    }

    ConstFactorization<> result{};
    result.length = N;
    std::size_t remainder = N;

    while ((remainder % 4U) == 0U) {
        if ((remainder / 4U) == 2U) {
            break;
        }
        result.factors[result.count++] = 4U;
        remainder /= 4U;
    }

    while ((remainder % 2U) == 0U) {
        result.factors[result.count++] = 2U;
        remainder /= 2U;
    }

    for (std::size_t prime = 3U; prime * prime <= remainder; prime += 2U) {
        while ((remainder % prime) == 0U) {
            result.factors[result.count++] = prime;
            remainder /= prime;
        }
    }

    if (remainder > 1U) {
        result.factors[result.count++] = remainder;
    }

    return result;
}

}  // namespace detail

/// Compile-time factorizer for fixed transform lengths known at compile time.
template <std::size_t N>
[[nodiscard]] consteval ConstFactorization<> factorize_ct() {
    return detail::factorize_length<N>();
}

/// In-place forward mixed-radix FFT (Singleton / Cooley-Tukey).
/// Convention: \f$X_k = \sum_{n=0}^{N-1} x_n \exp(-2\pi i kn / N)\f$.
template <FloatingPoint T>
void fft(std::span<std::complex<T>> data);

/// In-place inverse mixed-radix FFT with \f$1/N\f$ normalization.
template <FloatingPoint T>
void ifft(std::span<std::complex<T>> data);

}  // namespace mosaiq::fft
