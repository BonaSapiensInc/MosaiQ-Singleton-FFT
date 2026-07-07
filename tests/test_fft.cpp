#include <cstddef>
#include <cstdint>
#include <complex>
#include <random>
#include <type_traits>
#include <span>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "SingletonFFT.hpp"

namespace {

constexpr double k_double_tolerance = 1e-12;

template <typename T>
constexpr double tolerance_for() {
    if constexpr (std::is_same_v<T, float>) {
        return 1e-5;
    }
    return k_double_tolerance;
}

template <typename T>
std::vector<std::complex<T>> random_signal(std::size_t length, std::uint64_t seed) {
    std::vector<std::complex<T>> signal(length);
    std::mt19937_64 engine{seed};
    std::uniform_real_distribution<T> distribution(static_cast<T>(-1), static_cast<T>(1));

    for (auto& sample : signal) {
        sample = {distribution(engine), distribution(engine)};
    }

    return signal;
}

template <typename T>
void require_round_trip(std::size_t length, std::uint64_t seed) {
    auto original = random_signal<T>(length, seed);
    auto buffer = original;

    mosaiq::fft::fft(std::span<std::complex<T>>{buffer});
    mosaiq::fft::ifft(std::span<std::complex<T>>{buffer});

    const double tolerance = tolerance_for<T>();
    for (std::size_t index = 0U; index < length; ++index) {
        REQUIRE(buffer[index].real() == Catch::Approx(original[index].real()).margin(tolerance));
        REQUIRE(buffer[index].imag() == Catch::Approx(original[index].imag()).margin(tolerance));
    }
}

template <typename T>
void require_factorization_product(std::size_t length) {
    const mosaiq::fft::Factorization factors = mosaiq::fft::factorize(length);
    std::size_t product = 1U;
    for (const std::size_t radix : factors.factors) {
        product *= radix;
    }
    REQUIRE(product == length);
}

}  // namespace

TEST_CASE("factorize decomposes supported composite lengths", "[factorization]") {
    require_factorization_product<double>(64U);
    require_factorization_product<double>(81U);
    require_factorization_product<double>(442U);
}

TEST_CASE("constexpr factorize matches runtime factorize for 442", "[factorization]") {
    constexpr auto compile_time = mosaiq::fft::factorize_ct<442U>();
    const auto runtime = mosaiq::fft::factorize(442U);

    REQUIRE(compile_time == runtime);
}

TEST_CASE("forward-inverse FFT round trip (double)", "[fft][double]") {
    SECTION("pure power of two: 64") {
        require_round_trip<double>(64U, 0xC0FFEEULL);
    }

    SECTION("pure power of two: 128") {
        require_round_trip<double>(128U, 0xBADA55ULL);
    }

    SECTION("pure power of three: 81") {
        require_round_trip<double>(81U, 0xDEADBEEFULL);
    }

    SECTION("mixed prime length: 2 * 13 * 17 = 442") {
        require_round_trip<double>(442U, 0xFACEFEEDULL);
    }
}

TEST_CASE("forward-inverse FFT round trip (float)", "[fft][float]") {
    require_round_trip<float>(442U, 0x12345678ULL);
}

TEST_CASE("impulse response preserves unit energy", "[fft][impulse]") {
    std::vector<std::complex<double>> impulse(61U, {0.0, 0.0});
    impulse.front() = {1.0, 0.0};

    mosaiq::fft::fft(std::span{impulse});
    mosaiq::fft::ifft(std::span{impulse});

    REQUIRE(impulse.front().real() == Catch::Approx(1.0).margin(k_double_tolerance));
    REQUIRE(impulse.front().imag() == Catch::Approx(0.0).margin(k_double_tolerance));

    for (std::size_t index = 1U; index < impulse.size(); ++index) {
        REQUIRE(impulse[index].real() == Catch::Approx(0.0).margin(k_double_tolerance));
        REQUIRE(impulse[index].imag() == Catch::Approx(0.0).margin(k_double_tolerance));
    }
}
