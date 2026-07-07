#include "SingletonFFT.hpp"

#include "detail/PrimeRadix.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace mosaiq::fft {
namespace {

constexpr int k_forward_sign = -1;
constexpr int k_inverse_sign = 1;

template <typename T>
void validate_length(std::span<std::complex<T>> data) {
    if (data.size() < 2U) {
        throw std::invalid_argument("mosaiq::fft: transform length must be at least 2");
    }
}

template <typename T>
[[nodiscard]] std::vector<std::complex<T>> make_twiddles(std::size_t length, int sign) {
    std::vector<std::complex<T>> twiddles(length);
    for (std::size_t index = 0U; index < length; ++index) {
        const T angle = static_cast<T>(sign) * detail::two_pi<T>() * static_cast<T>(index) /
                        static_cast<T>(length);
        twiddles[index] = {std::cos(angle), std::sin(angle)};
    }
    return twiddles;
}

/// KissFFT-style factor buffer: [radix, remainder, radix, remainder, ..., 1].
[[nodiscard]] std::vector<std::size_t> build_stage_factors(std::size_t length) {
    std::vector<std::size_t> buffer;
    std::size_t remainder = length;
    std::size_t radix = 4U;
    const std::size_t floor_sqrt =
        static_cast<std::size_t>(std::floor(std::sqrt(static_cast<double>(remainder))));

    do {
        while ((remainder % radix) != 0U) {
            switch (radix) {
                case 4U:
                    radix = 2U;
                    break;
                case 2U:
                    radix = 3U;
                    break;
                default:
                    radix += 2U;
                    break;
            }
            if (radix > floor_sqrt) {
                radix = remainder;
            }
        }

        remainder /= radix;
        buffer.push_back(radix);
        buffer.push_back(remainder);
    } while (remainder > 1U);

    return buffer;
}

template <typename T>
void butterfly_radix2(std::complex<T>* output,
                      std::size_t stride,
                      std::span<const std::complex<T>> twiddles,
                      int group_count) {
    std::complex<T>* leg_one = output + group_count;
    std::size_t twiddle_index = 0U;

    for (int group = 0; group < group_count; ++group) {
        std::array<std::complex<T>, 2U> scratch{};
        scratch[0] = output[group];
        scratch[1] = leg_one[group] * twiddles[twiddle_index];
        detail::PrimeRadix<T, 2>::apply(std::span{scratch}, 0);
        output[group] = scratch[0];
        leg_one[group] = scratch[1];
        twiddle_index += stride;
    }
}

template <typename T>
void butterfly_prime_radix(std::complex<T>* output,
                           std::size_t stride,
                           std::span<const std::complex<T>> twiddles,
                           std::size_t length,
                           int group_count,
                           int radix,
                           int sign) {
    std::array<std::complex<T>, detail::k_max_prime_radix> scratch{};

    for (int group = 0; group < group_count; ++group) {
        for (int branch = 0; branch < radix; ++branch) {
            scratch[static_cast<std::size_t>(branch)] = output[group + branch * group_count];
        }

        for (int branch = 1; branch < radix; ++branch) {
            int twiddle_index = static_cast<int>(stride) * group * branch;
            twiddle_index %= static_cast<int>(length);
            scratch[static_cast<std::size_t>(branch)] *=
                twiddles[static_cast<std::size_t>(twiddle_index)];
        }

        detail::dispatch_prime_radix(std::span<std::complex<T>>{scratch}.first(static_cast<std::size_t>(radix)),
                                     static_cast<std::size_t>(radix),
                                     sign);

        for (int output_bin = 0; output_bin < radix; ++output_bin) {
            output[group + output_bin * group_count] =
                scratch[static_cast<std::size_t>(output_bin)];
        }
    }
}

template <typename T>
void butterfly_generic(std::complex<T>* output,
                       std::size_t stride,
                       std::span<const std::complex<T>> twiddles,
                       std::size_t length,
                       int group_count,
                       int radix) {
    std::vector<std::complex<T>> scratch(static_cast<std::size_t>(radix));

    for (int group = 0; group < group_count; ++group) {
        int index = group;
        for (int branch = 0; branch < radix; ++branch) {
            scratch[static_cast<std::size_t>(branch)] = output[index];
            index += group_count;
        }

        index = group;
        for (int output_bin = 0; output_bin < radix; ++output_bin) {
            int twiddle_index = 0;
            output[index] = scratch[0U];
            for (int branch = 1; branch < radix; ++branch) {
                twiddle_index += static_cast<int>(stride) * index;
                if (twiddle_index >= static_cast<int>(length)) {
                    twiddle_index -= static_cast<int>(length);
                }
                output[index] += scratch[static_cast<std::size_t>(branch)] *
                                 twiddles[static_cast<std::size_t>(twiddle_index)];
            }
            index += group_count;
        }
    }
}

template <typename T>
void butterfly_dispatch(std::complex<T>* output,
                        std::size_t stride,
                        std::span<const std::complex<T>> twiddles,
                        std::size_t length,
                        int group_count,
                        int radix,
                        int sign) {
    if (radix == 2) {
        butterfly_radix2(output, stride, twiddles, group_count);
        return;
    }

    if (detail::is_supported_radix(static_cast<std::size_t>(radix))) {
        butterfly_prime_radix(output, stride, twiddles, length, group_count, radix, sign);
        return;
    }

    butterfly_generic(output, stride, twiddles, length, group_count, radix);
}

template <typename T>
void transform_pass(std::complex<T>* output,
                    const std::complex<T>* input,
                    std::size_t stride,
                    std::size_t input_stride,
                    std::span<const std::size_t> stage_factors,
                    std::span<const std::complex<T>> twiddles,
                    std::size_t length,
                    int sign) {
    const int radix = static_cast<int>(stage_factors.front());
    const int group_count = static_cast<int>(stage_factors[1U]);
    const std::size_t next_stride = stride * static_cast<std::size_t>(radix);
    std::complex<T>* const pass_begin = output;
    std::complex<T>* const pass_end = output + static_cast<std::size_t>(radix * group_count);

    if (group_count == 1) {
        std::complex<T>* cursor = output;
        const std::complex<T>* sample = input;
        do {
            *cursor = *sample;
            sample += stride * input_stride;
            ++cursor;
        } while (cursor != pass_end);
    } else {
        std::complex<T>* cursor = output;
        const std::complex<T>* sample = input;
        do {
            transform_pass(cursor,
                           sample,
                           next_stride,
                           input_stride,
                           stage_factors.subspan(2U),
                           twiddles,
                           length,
                           sign);
            sample += stride * input_stride;
            cursor += static_cast<std::size_t>(group_count);
        } while (cursor != pass_end);
    }

    butterfly_dispatch(pass_begin, stride, twiddles, length, group_count, radix, sign);
}

template <typename T>
void transform(std::span<std::complex<T>> data, int sign, bool normalize) {
    validate_length(data);

    const std::vector<std::size_t> stage_factors = build_stage_factors(data.size());
    if (stage_factors.empty()) {
        throw std::invalid_argument("mosaiq::fft: unable to factorize transform length");
    }

    const std::vector<std::complex<T>> twiddles = make_twiddles<T>(data.size(), sign);
    std::vector<std::complex<T>> scratch(data.size());

    transform_pass(scratch.data(),
                   data.data(),
                   1U,
                   1U,
                   std::span<const std::size_t>{stage_factors.data(), stage_factors.size()},
                   std::span<const std::complex<T>>{twiddles.data(), twiddles.size()},
                   data.size(),
                   sign);

    for (std::size_t index = 0U; index < data.size(); ++index) {
        data[index] = scratch[index];
    }

    if (normalize) {
        const T scale = static_cast<T>(1) / static_cast<T>(data.size());
        for (auto& sample : data) {
            sample *= scale;
        }
    }
}

}  // namespace

Factorization factorize(std::size_t n) {
    if (n < 2U) {
        return Factorization{.factors = {}, .length = n};
    }

    Factorization result{};
    result.length = n;
    std::size_t remainder = n;

    while ((remainder % 4U) == 0U) {
        if ((remainder / 4U) == 2U) {
            break;
        }
        result.factors.push_back(4U);
        remainder /= 4U;
    }

    while ((remainder % 2U) == 0U) {
        result.factors.push_back(2U);
        remainder /= 2U;
    }

    for (std::size_t prime = 3U; prime * prime <= remainder; prime += 2U) {
        while ((remainder % prime) == 0U) {
            result.factors.push_back(prime);
            remainder /= prime;
        }
    }

    if (remainder > 1U) {
        if (!detail::is_supported_radix(remainder)) {
            throw std::invalid_argument("mosaiq::fft: unsupported prime radix in factorization");
        }
        result.factors.push_back(remainder);
    }

    for (const std::size_t radix : result.factors) {
        if (!detail::is_supported_radix(radix)) {
            throw std::invalid_argument("mosaiq::fft: unsupported radix in factorization");
        }
    }

    std::size_t product = 1U;
    for (const std::size_t radix : result.factors) {
        product *= radix;
    }

    if (product != n) {
        throw std::invalid_argument("mosaiq::fft: internal factorization mismatch");
    }

    return result;
}

template <FloatingPoint T>
void fft(std::span<std::complex<T>> data) {
    transform(data, k_forward_sign, false);
}

template <FloatingPoint T>
void ifft(std::span<std::complex<T>> data) {
    transform(data, k_inverse_sign, true);
}

template void fft<float>(std::span<std::complex<float>> data);
template void fft<double>(std::span<std::complex<double>> data);
template void fft<long double>(std::span<std::complex<long double>> data);

template void ifft<float>(std::span<std::complex<float>> data);
template void ifft<double>(std::span<std::complex<double>> data);
template void ifft<long double>(std::span<std::complex<long double>> data);

}  // namespace mosaiq::fft
