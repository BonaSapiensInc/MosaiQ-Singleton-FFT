#include "SingletonFFT.hpp"

#include "detail/PrimeRadix.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

namespace mosaiq::fft {
namespace {

constexpr int k_forward_sign = -1;
constexpr int k_inverse_sign = 1;

struct PrimeMetrics {
    std::ptrdiff_t maxn{0};  // total count of prime factors
    std::ptrdiff_t maxf{0};  // largest prime factor (<= 61)
    std::ptrdiff_t maxp{0};  // product of square-free prime factors
};

/// Computes dynamic workspace bounds (maxn, maxf, maxp) per I. G. Kim (2012)
/// while strictly enforcing the 61-smooth prime radix semigroup.
[[nodiscard]] PrimeMetrics analyze_61_smooth(std::size_t n) {
    if (n < 2U) {
        throw std::invalid_argument("mosaiq::fft: transform length must be at least 2");
    }

    static constexpr std::size_t k_primes[] = {
        2U, 3U, 5U, 7U, 11U, 13U, 17U, 19U, 23U,
        29U, 31U, 37U, 41U, 43U, 47U, 53U, 59U, 61U
    };

    std::size_t rem = n;
    std::size_t total_factors = 0U;
    std::size_t largest_prime = 2U;
    std::size_t square_free_prod = 1U;

    for (const std::size_t p : k_primes) {
        std::size_t count = 0U;
        while ((rem % p) == 0U) {
            rem /= p;
            ++count;
        }
        if (count > 0U) {
            total_factors += count;
            largest_prime = p;
            if ((count % 2U) != 0U) {
                square_free_prod *= p;
            }
        }
    }

    if (rem != 1U) {
        throw std::invalid_argument(
            "mosaiq::fft: transform length exceeds 61-smooth prime radix domain");
    }

    // Radix-4 pairs consume workspace up to factor 4 when n is a pure power of 2
    const std::size_t effective_maxf = std::max<std::size_t>(largest_prime, 4U);

    std::ptrdiff_t maxn = static_cast<std::ptrdiff_t>(total_factors);
    std::ptrdiff_t maxf = static_cast<std::ptrdiff_t>(effective_maxf);
    std::ptrdiff_t maxp = static_cast<std::ptrdiff_t>(square_free_prod);

    if (maxp <= maxn) {
        maxp = maxn + 1;
    }
    if ((maxn % 2) == 0) {
        maxn += 1;  // Singleton invariant: maxn assumed odd
    }

    return PrimeMetrics{.maxn = maxn, .maxf = maxf, .maxp = maxp};
}

/// Faithful C++20 translation of R. C. Singleton (1968) / I. G. Kim (2006, 2012)
/// non-recursive, true in-place multivariate mixed-radix FFT kernel (`cfft.F`).
/// Supports general strided dimensions via (ntot, n, nspan).
template <typename T>
void singleton_cfft_strided(std::span<std::complex<T>> data,
                            std::size_t ntot_u,
                            std::size_t n_u,
                            std::size_t nspan_u,
                            int isgn) {
    if (n_u < 2U) {
        return;
    }

    const PrimeMetrics bounds = analyze_61_smooth(n_u);
    std::ptrdiff_t maxn = bounds.maxn;
    std::ptrdiff_t maxf = bounds.maxf;
    std::ptrdiff_t maxp = bounds.maxp;

    // Dynamic O(maxf + maxp) auxiliary workspace — zero O(N) data scratch buffer!
    std::vector<std::ptrdiff_t> nfac_vec(static_cast<std::size_t>(maxn + 2), 0);
    std::vector<T> at_vec(static_cast<std::size_t>(maxf + 1), T{0});
    std::vector<T> ck_vec(static_cast<std::size_t>(maxf + 1), T{0});
    std::vector<T> bt_vec(static_cast<std::size_t>(maxf + 1), T{0});
    std::vector<T> sk_vec(static_cast<std::size_t>(maxf + 1), T{0});
    std::vector<std::ptrdiff_t> np_vec(static_cast<std::size_t>(maxp + 2), 0);

    // 1-based index accessors matching Fortran `cfft.F` semantics cleanly
    auto re = [&](std::ptrdiff_t idx) noexcept -> T& {
        return reinterpret_cast<T(&)[2]>(data[static_cast<std::size_t>(idx - 1)])[0];
    };
    auto im = [&](std::ptrdiff_t idx) noexcept -> T& {
        return reinterpret_cast<T(&)[2]>(data[static_cast<std::size_t>(idx - 1)])[1];
    };
    auto nfac = [&](std::ptrdiff_t idx) noexcept -> std::ptrdiff_t& {
        return nfac_vec[static_cast<std::size_t>(idx)];
    };
    auto at = [&](std::ptrdiff_t idx) noexcept -> T& {
        return at_vec[static_cast<std::size_t>(idx)];
    };
    auto bt = [&](std::ptrdiff_t idx) noexcept -> T& {
        return bt_vec[static_cast<std::size_t>(idx)];
    };
    auto ck = [&](std::ptrdiff_t idx) noexcept -> T& {
        return ck_vec[static_cast<std::size_t>(idx)];
    };
    auto sk = [&](std::ptrdiff_t idx) noexcept -> T& {
        return sk_vec[static_cast<std::size_t>(idx)];
    };
    auto np = [&](std::ptrdiff_t idx) noexcept -> std::ptrdiff_t& {
        return np_vec[static_cast<std::size_t>(idx)];
    };

    const std::ptrdiff_t ntot = static_cast<std::ptrdiff_t>(ntot_u);
    const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(n_u);
    const std::ptrdiff_t nspan = static_cast<std::ptrdiff_t>(nspan_u);

    std::ptrdiff_t inc = isgn;
    const T pi = std::numbers::pi_v<T>;
    T rad = static_cast<T>(2) * pi;
    T s72 = std::sin(static_cast<T>(0.4) * pi);
    const T c72 = std::cos(static_cast<T>(0.4) * pi);
    T s120 = std::sqrt(static_cast<T>(0.75));

    if (isgn < 0) {
        s72 = -s72;
        s120 = -s120;
        rad = -rad;
        inc = -inc;
    }

    std::ptrdiff_t nt = inc * ntot;
    const std::ptrdiff_t ks = inc * nspan;
    std::ptrdiff_t kspan = ks;
    std::ptrdiff_t nn = nt - inc;
    std::ptrdiff_t jc = ks / n;
    const T radf = rad * static_cast<T>(jc) * static_cast<T>(0.5);
    std::ptrdiff_t i = 0;
    std::ptrdiff_t jf = 0;

    // -----------------------------------------------------------------------
    // Symmetric factor decomposition of n: [square factors | square-free | reversed square]
    // -----------------------------------------------------------------------
    std::ptrdiff_t m = 0;
    std::ptrdiff_t k = n;

    while ((k - (k / 16) * 16) == 0) {
        ++m;
        nfac(m) = 4;
        k /= 16;
    }

    std::ptrdiff_t j = 3;
    std::ptrdiff_t jj = 9;
    for (;;) {
        while ((k % jj) == 0) {
            ++m;
            nfac(m) = j;
            k /= jj;
        }
        j += 2;
        jj = j * j;
        if (jj > k) {
            break;
        }
    }

    std::ptrdiff_t kt = 0;
    if (k > 4) {
        if ((k - (k / 4) * 4) == 0) {
            ++m;
            nfac(m) = 2;
            k /= 4;
        }
        kt = m;
        j = 2;
        for (;;) {
            if ((k % j) == 0) {
                ++m;
                nfac(m) = j;
                k /= j;
            }
            j = ((j + 1) / 2) * 2 + 1;
            if (j > k) {
                break;
            }
        }
    } else {
        kt = m;
        nfac(m + 1) = k;
        if (k != 1) {
            ++m;
        }
    }

    if (kt != 0) {
        j = kt;
        for (std::ptrdiff_t jx = j; jx >= 1; --jx) {
            ++m;
            nfac(m) = nfac(jx);
        }
    }

    // -----------------------------------------------------------------------
    // Main non-recursive mixed-radix butterfly passes
    // -----------------------------------------------------------------------
    std::ptrdiff_t kspnn = 0;
    std::ptrdiff_t k1 = 0, k2 = 0, k3 = 0, k4 = 0, kk = 0;
    T sd{}, cd{}, c1{}, c2{}, c3{}, s1{}, s2{}, s3{};
    T aa{}, bb{}, ak{}, bk{}, aj{}, bj{}, akp{}, akm{}, bkp{}, bkm{}, ajp{}, ajm{}, bjp{}, bjm{};

    for (;;) {
        sd = radf / static_cast<T>(kspan);
        const T sin_sd = std::sin(sd);
        cd = static_cast<T>(2) * sin_sd * sin_sd;
        sd = std::sin(sd + sd);
        kk = 1;
        ++i;

        if (nfac(i) != 2) {
            if (nfac(i) != 4) {
                // ---> Transform for odd prime factors (3, 5, 7..61)
                k = nfac(i);
                kspnn = kspan;
                kspan /= k;

                if (k == 3) {
                    for (;;) {
                        k1 = kk + kspan;
                        k2 = k1 + kspan;
                        ak = re(kk);
                        bk = im(kk);
                        aj = re(k1) + re(k2);
                        bj = im(k1) + im(k2);
                        re(kk) = ak + aj;
                        im(kk) = bk + bj;
                        ak = static_cast<T>(-0.5) * aj + ak;
                        bk = static_cast<T>(-0.5) * bj + bk;
                        aj = (re(k1) - re(k2)) * s120;
                        bj = (im(k1) - im(k2)) * s120;
                        re(k1) = ak - bj;
                        im(k1) = bk + aj;
                        re(k2) = ak + bj;
                        im(k2) = bk - aj;
                        kk = k2 + kspan;
                        if (kk < nn) {
                            continue;
                        }
                        kk -= nn;
                        if (kk > kspan) {
                            break;
                        }
                    }
                } else if (k == 5) {
                    c2 = c72 * c72 - s72 * s72;
                    s2 = static_cast<T>(2) * c72 * s72;
                    for (;;) {
                        k1 = kk + kspan;
                        k2 = k1 + kspan;
                        k3 = k2 + kspan;
                        k4 = k3 + kspan;
                        akp = re(k1) + re(k4);
                        akm = re(k1) - re(k4);
                        bkp = im(k1) + im(k4);
                        bkm = im(k1) - im(k4);
                        ajp = re(k2) + re(k3);
                        ajm = re(k2) - re(k3);
                        bjp = im(k2) + im(k3);
                        bjm = im(k2) - im(k3);
                        aa = re(kk);
                        bb = im(kk);
                        re(kk) = aa + akp + ajp;
                        im(kk) = bb + bkp + bjp;
                        ak = akp * c72 + ajp * c2 + aa;
                        bk = bkp * c72 + bjp * c2 + bb;
                        aj = akm * s72 + ajm * s2;
                        bj = bkm * s72 + bjm * s2;
                        re(k1) = ak - bj;
                        re(k4) = ak + bj;
                        im(k1) = bk + aj;
                        im(k4) = bk - aj;
                        ak = akp * c2 + ajp * c72 + aa;
                        bk = bkp * c2 + bjp * c72 + bb;
                        aj = akm * s2 - ajm * s72;
                        bj = bkm * s2 - bjm * s72;
                        re(k2) = ak - bj;
                        re(k3) = ak + bj;
                        im(k2) = bk + aj;
                        im(k3) = bk - aj;
                        kk = k4 + kspan;
                        if (kk < nn) {
                            continue;
                        }
                        kk -= nn;
                        if (kk > kspan) {
                            break;
                        }
                    }
                } else {
                    if (k != jf) {
                        jf = k;
                        s1 = rad / static_cast<T>(k);
                        c1 = std::cos(s1);
                        s1 = std::sin(s1);
                        ck(jf) = static_cast<T>(1);
                        sk(jf) = static_cast<T>(0);
                        j = 1;
                        for (;;) {
                            ck(j) = ck(k) * c1 + sk(k) * s1;
                            sk(j) = ck(k) * s1 - sk(k) * c1;
                            ck(k - 1) = ck(j);
                            sk(k - 1) = -sk(j);
                            --k;
                            ++j;
                            if (j >= k) {
                                break;
                            }
                        }
                    }
                    for (;;) {
                        k1 = kk;
                        k2 = kk + kspnn;
                        aa = re(kk);
                        bb = im(kk);
                        ak = aa;
                        bk = bb;
                        j = 1;
                        k1 += kspan;
                        for (;;) {
                            k2 -= kspan;
                            ++j;
                            at(j) = re(k1) + re(k2);
                            ak += at(j);
                            bt(j) = im(k1) + im(k2);
                            bk += bt(j);
                            ++j;
                            at(j) = re(k1) - re(k2);
                            bt(j) = im(k1) - im(k2);
                            k1 += kspan;
                            if (k1 >= k2) {
                                break;
                            }
                        }
                        re(kk) = ak;
                        im(kk) = bk;
                        k1 = kk;
                        k2 = kk + kspnn;
                        j = 1;
                        for (;;) {
                            k1 += kspan;
                            k2 -= kspan;
                            jj = j;
                            ak = aa;
                            bk = bb;
                            aj = static_cast<T>(0);
                            bj = static_cast<T>(0);
                            k = 1;
                            for (;;) {
                                ++k;
                                ak += at(k) * ck(jj);
                                bk += bt(k) * ck(jj);
                                ++k;
                                aj += at(k) * sk(jj);
                                bj += bt(k) * sk(jj);
                                jj += j;
                                if (jj > jf) {
                                    jj -= jf;
                                }
                                if (k >= jf) {
                                    break;
                                }
                            }
                            k = jf - j;
                            re(k1) = ak - bj;
                            im(k1) = bk + aj;
                            re(k2) = ak + bj;
                            im(k2) = bk - aj;
                            ++j;
                            if (j >= k) {
                                break;
                            }
                        }
                        kk += kspnn;
                        if (kk <= nn) {
                            continue;
                        }
                        kk -= nn;
                        if (kk > kspan) {
                            break;
                        }
                    }
                }

                // Multiply by rotation factor (for odd factors)
                if (i == m) {
                    break;
                }
                kk = jc + 1;
                for (;;) {
                    c2 = static_cast<T>(1) - cd;
                    s1 = sd;
                    for (;;) {
                        c1 = c2;
                        s2 = s1;
                        kk += kspan;
                        for (;;) {
                            ak = re(kk);
                            re(kk) = c2 * ak - s2 * im(kk);
                            im(kk) = s2 * ak + c2 * im(kk);
                            kk += kspnn;
                            if (kk <= nt) {
                                continue;
                            }
                            ak = s1 * s2;
                            s2 = s1 * c2 + c1 * s2;
                            c2 = c1 * c2 - ak;
                            kk = kk - nt + kspan;
                            if (kk > kspnn) {
                                break;
                            }
                        }
                        c2 = c1 - (cd * c1 + sd * s1);
                        s1 = s1 + (sd * c1 - cd * s1);
                        // Truncation error renormalization
                        c1 = static_cast<T>(0.5) / (c2 * c2 + s1 * s1) + static_cast<T>(0.5);
                        s1 *= c1;
                        c2 *= c1;
                        kk = kk - kspnn + jc;
                        if (kk > kspan) {
                            break;
                        }
                    }
                    kk = kk - kspan + jc + inc;
                    if (kk > (jc + jc)) {
                        break;
                    }
                }
                continue;
            }

            // ---> Transform for factor of 4
            kspnn = kspan;
            kspan /= 4;
            for (;;) {
                c1 = static_cast<T>(1);
                s1 = static_cast<T>(0);
                for (;;) {
                    k1 = kk + kspan;
                    k2 = k1 + kspan;
                    k3 = k2 + kspan;
                    akp = re(kk) + re(k2);
                    akm = re(kk) - re(k2);
                    ajp = re(k1) + re(k3);
                    ajm = re(k1) - re(k3);
                    re(kk) = akp + ajp;
                    ajp = akp - ajp;
                    bkp = im(kk) + im(k2);
                    bkm = im(kk) - im(k2);
                    bjp = im(k1) + im(k3);
                    bjm = im(k1) - im(k3);
                    im(kk) = bkp + bjp;
                    bjp = bkp - bjp;
                    if (isgn < 0) {
                        akp = akm + bjm;
                        akm = akm - bjm;
                        bkp = bkm - ajm;
                        bkm = bkm + ajm;
                    } else {
                        akp = akm - bjm;
                        akm = akm + bjm;
                        bkp = bkm + ajm;
                        bkm = bkm - ajm;
                    }
                    if (s1 == static_cast<T>(0)) {
                        re(k1) = akp;
                        im(k1) = bkp;
                        re(k2) = ajp;
                        im(k2) = bjp;
                        re(k3) = akm;
                        im(k3) = bkm;
                    } else {
                        re(k1) = akp * c1 - bkp * s1;
                        im(k1) = akp * s1 + bkp * c1;
                        re(k2) = ajp * c2 - bjp * s2;
                        im(k2) = ajp * s2 + bjp * c2;
                        re(k3) = akm * c3 - bkm * s3;
                        im(k3) = akm * s3 + bkm * c3;
                    }
                    kk = k3 + kspan;
                    if (kk <= nt) {
                        continue;
                    }
                    c2 = c1 - (cd * c1 + sd * s1);
                    s1 = (sd * c1 - cd * s1) + s1;
                    c1 = static_cast<T>(0.5) / (c2 * c2 + s1 * s1) + static_cast<T>(0.5);
                    s1 *= c1;
                    c1 *= c2;
                    c2 = c1 * c1 - s1 * s1;
                    s2 = static_cast<T>(2) * c1 * s1;
                    c3 = c2 * c1 - s2 * s1;
                    s3 = c2 * s1 + s2 * c1;
                    kk = kk - nt + jc;
                    if (kk > kspan) {
                        break;
                    }
                }
                kk = kk - kspan + inc;
                if (kk > jc) {
                    break;
                }
            }
            if (kspan == jc) {
                break;
            }
            continue;
        }

        // ---> Transform for factor of 2
        kspan /= 2;
        k1 = kspan + 2;
        for (;;) {
            k2 = kk + kspan;
            ak = re(k2);
            bk = im(k2);
            re(k2) = re(kk) - ak;
            im(k2) = im(kk) - bk;
            re(kk) += ak;
            im(kk) += bk;
            kk = k2 + kspan;
            if (kk <= nn) {
                continue;
            }
            kk -= nn;
            if (kk > jc) {
                break;
            }
        }
        if (kk > kspan) {
            break;
        }
        for (;;) {
            c1 = static_cast<T>(1) - cd;
            s1 = sd;
            for (;;) {
                k2 = kk + kspan;
                ak = re(kk) - re(k2);
                bk = im(kk) - im(k2);
                re(kk) += re(k2);
                im(kk) += im(k2);
                re(k2) = c1 * ak - s1 * bk;
                im(k2) = s1 * ak + c1 * bk;
                kk = k2 + kspan;
                if (kk < nt) {
                    continue;
                }
                k2 = kk - nt;
                c1 = -c1;
                kk = k1 - k2;
                if (kk > k2) {
                    continue;
                }
                ak = c1 - (cd * c1 + sd * s1);
                s1 = (sd * c1 - cd * s1) + s1;
                c1 = static_cast<T>(0.5) / (ak * ak + s1 * s1) + static_cast<T>(0.5);
                s1 *= c1;
                c1 *= ak;
                kk += jc;
                if (kk >= k2) {
                    break;
                }
            }
            k1 = k1 + inc + inc;
            kk = (k1 - kspan) / 2 + jc;
            if (kk > (jc + jc)) {
                break;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Two-stage in-place permutation to normal order
    // Stage 1: Permutation for square factors of n
    // -----------------------------------------------------------------------
    np(1) = ks;
    if (kt != 0) {
        k = kt + kt + 1;
        if (m < k) {
            --k;
        }
        j = 1;
        np(k + 1) = jc;
        for (;;) {
            np(j + 1) = np(j) / nfac(j);
            np(k) = np(k + 1) * nfac(j);
            ++j;
            --k;
            if (j >= k) {
                break;
            }
        }
        k3 = np(k + 1);
        kspan = np(2);
        kk = jc + 1;
        k2 = kspan + 1;
        j = 1;

        if (n == ntot) {
            for (;;) {
                for (;;) {
                    std::swap(re(kk), re(k2));
                    std::swap(im(kk), im(k2));
                    kk += inc;
                    k2 += kspan;
                    if (k2 >= ks) {
                        break;
                    }
                }
                bool restart_outer = false;
                for (;;) {
                    k2 -= np(j);
                    ++j;
                    k2 += np(j + 1);
                    if (k2 > np(j)) {
                        continue;
                    }
                    j = 1;
                    for (;;) {
                        if (kk < k2) {
                            restart_outer = true;
                            break;
                        }
                        kk += inc;
                        k2 += kspan;
                        if (k2 >= ks) {
                            break;
                        }
                    }
                    if (restart_outer || kk >= ks) {
                        break;
                    }
                }
                if (!restart_outer) {
                    break;
                }
            }
            jc = k3;
        } else {
            for (;;) {
                k = kk + jc;
                std::ptrdiff_t kkx = kk;
                for (; kkx <= k - 1; kkx += inc) {
                    std::swap(re(kkx), re(k2));
                    std::swap(im(kkx), im(k2));
                    k2 += inc;
                }
                kk = kkx + ks - jc;
                k2 = k2 + ks - jc;
                if (kk < nt) {
                    continue;
                }
                k2 = k2 - nt + kspan;
                kk = kk - nt + jc;
                if (k2 < ks) {
                    continue;
                }
                bool restart_outer = false;
                for (;;) {
                    k2 -= np(j);
                    ++j;
                    k2 += np(j + 1);
                    if (k2 > np(j)) {
                        continue;
                    }
                    j = 1;
                    for (;;) {
                        if (kk < k2) {
                            restart_outer = true;
                            break;
                        }
                        kk += jc;
                        k2 += kspan;
                        if (k2 >= ks) {
                            break;
                        }
                    }
                    if (restart_outer || kk >= ks) {
                        break;
                    }
                }
                if (!restart_outer) {
                    break;
                }
            }
            jc = k3;
        }
    }

    if ((2 * kt + 1) >= m) {
        return;
    }

    // -----------------------------------------------------------------------
    // Stage 2: Cycle-following permutation for square-free factors of n
    // -----------------------------------------------------------------------
    kspnn = np(kt + 1);
    j = m - kt;
    nfac(j + 1) = 1;
    for (;;) {
        nfac(j) *= nfac(j + 1);
        --j;
        if (j == kt) {
            break;
        }
    }

    ++kt;
    nn = nfac(kt) - 1;
    if (nn > maxp) {
        throw std::logic_error("mosaiq::fft: square-free cycle workspace bound exceeded");
    }

    jj = 0;
    j = 0;
    for (;;) {
        k2 = nfac(kt);
        k = kt + 1;
        kk = nfac(k);
        ++j;
        if (j > nn) {
            break;
        }
        for (;;) {
            jj += kk;
            if (jj < k2) {
                np(j) = jj;
                break;
            }
            jj -= k2;
            k2 = kk;
            ++k;
            kk = nfac(k);
        }
    }

    j = 0;
    for (;;) {
        ++j;
        kk = np(j);
        if (kk < 0) {
            continue;
        }
        if (kk != j) {
            for (;;) {
                k = kk;
                kk = np(k);
                np(k) = -kk;
                if (kk == j) {
                    break;
                }
            }
            k3 = kk;
            continue;
        }
        np(j) = -j;
        if (j == nn) {
            break;
        }
    }

    maxf *= inc;
    for (;;) {
        j = k3 + 1;
        nt -= kspnn;
        const std::ptrdiff_t ii = nt - inc + 1;
        if (nt < 0) {
            break;
        }
        for (;;) {
            --j;
            if (np(j) < 0) {
                continue;
            }
            jj = jc;
            for (;;) {
                kspan = jj;
                if (jj > maxf) {
                    kspan = maxf;
                }
                jj -= kspan;
                k = np(j);
                kk = jc * k + ii + jj;
                k1 = kk + kspan;
                k2 = 0;
                for (;;) {
                    ++k2;
                    at(k2) = re(k1);
                    bt(k2) = im(k1);
                    k1 -= inc;
                    if (k1 == kk) {
                        break;
                    }
                }
                for (;;) {
                    k1 = kk + kspan;
                    k2 = k1 - jc * (k + np(k));
                    k = -np(k);
                    for (;;) {
                        re(k1) = re(k2);
                        im(k1) = im(k2);
                        k1 -= inc;
                        k2 -= inc;
                        if (k1 == kk) {
                            break;
                        }
                    }
                    kk = k2;
                    if (k == j) {
                        break;
                    }
                }
                k1 = kk + kspan;
                k2 = 0;
                for (;;) {
                    ++k2;
                    re(k1) = at(k2);
                    im(k1) = bt(k2);
                    k1 -= inc;
                    if (k1 == kk) {
                        break;
                    }
                }
                if (jj == 0) {
                    break;
                }
            }
            if (j == 1) {
                break;
            }
        }
    }
}

template <typename T>
void transform(std::span<std::complex<T>> data, int sign, bool normalize) {
    if (data.size() < 2U) {
        throw std::invalid_argument("mosaiq::fft: transform length must be at least 2");
    }

    singleton_cfft_strided(data, data.size(), data.size(), data.size(), sign);

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
