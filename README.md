# MosaiQ-Singleton-FFT

![C++20](https://img.shields.io/badge/C++-20-blue.svg)
![CMake](https://img.shields.io/badge/CMake-3.20+-064F8C.svg)
![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-blue.svg)
![Commercial: \$199](https://img.shields.io/badge/Commercial-%24199-brightgreen.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21236967.svg)](https://doi.org/10.5281/zenodo.21236967)

A strictly deterministic C++20 realization of R. C. Singleton's **true in-place**, non-recursive mixed-radix FFT, modernized by BonaSapiens Inc. over the **61-smooth** prime radix domain. The public API is span-typed; the transform kernel carries no \(O(N)\) data scratch buffer and no recursive Cooley–Tukey / KissFFT fallback.

---

## Overview

`MosaiQ-Singleton-FFT` factors a transform length \(N\) into Singleton's symmetric square / square-free schedule and executes non-recursive butterfly passes with a two-stage cycle-following permutation. Auxiliary workspace is bounded by \(O(\mathrm{maxf} + \mathrm{maxp})\), not \(O(N)\).

Supported radices (61-smooth domain):

\[
\{2,\,3,\,4,\,5,\,7,\,11,\,13,\,17,\,19,\,23,\,29,\,31,\,37,\,41,\,43,\,47,\,53,\,59,\,61\}
\]

Forward transforms use the convention

\[
X_k = \sum_{n=0}^{N-1} x_n \exp\!\left(-\frac{2\pi i\, kn}{N}\right).
\]

Inverse transforms apply the conjugate exponential and scale by \(1/N\).

| Component | Location |
|-----------|----------|
| Public API | `src/SingletonFFT.hpp` |
| Singleton transform engine | `src/SingletonFFT.cpp` |
| 61-smooth radix predicate | `src/detail/PrimeRadix.hpp` |
| Validation suite | `tests/test_fft.cpp` |
| Legacy provenance archive | `docs/legacy/` |

---

## Lineage & True In-Place Architecture

Unlike generic wrappers that silently fall back to recursive Cooley–Tukey / KissFFT passes with \(O(N)\) heap scratch allocations, `MosaiQ-Singleton-FFT` is a pure C++20 realization of the non-recursive, **true in-place** multivariate mixed-radix algorithm:

* **R. C. Singleton** (SRI, 1968/1969) — Symmetric square/square-free factorization and two-stage cycle-following in-place permutation.
* **I. G. Kim** (*J. Korean Magn. Soc.* **22**(2), 33, 2012) — Structured `goto`-free reformulation with dynamic \(O(\mathrm{maxf} + \mathrm{maxp})\) auxiliary workspace bounds (`maxn`, `maxf`, `maxp`), eliminating \(O(N)\) copy-in/copy-out overhead across multidimensional grids. DOI: [10.4283/JKMS.2012.22.2.033](https://doi.org/10.4283/JKMS.2012.22.2.033).
* **BonaSapiens C++20 Engine (2026)** — Strictly typed `std::span` interface over the **61-smooth** prime radix domain with recurrence renormalization of rotation factors.

---

## Mathematical Note: 61-Smooth Grids & The Riemann Hypothesis

*(Commemorative Refactoring Note — November 2026)*

This unification into a pure non-recursive Singleton pipeline commemorates recent mathematical discourse surrounding the **Riemann Hypothesis (RH)**, highlighting two asymptotic intersections between discrete Fourier grid geometry and analytic number theory:

1. **Deterministic Immunity to GRH:**
   Supporting arbitrary large prime lengths \(p\) in \(O(p \log p)\) time via Rader's algorithm requires finding a primitive root modulo \(p\). Guaranteeing a deterministic polynomial-time bound \(O(\log^6 p)\) for the least primitive root relies on the **Generalized Riemann Hypothesis (GRH)** (Ankeny–Burgess–Wang theorem). By bounding the prime radix semigroup strictly to the 18 primes \(\le 61\) (`61-smooth` numbers) and purging ad-hoc radix-2 recursive branches, `MosaiQ-Singleton-FFT` guarantees deterministic complexity and invariant in-place permutation cycles unconditionally—independent of GRH.
2. **Smooth-Number Grid Gaps & Zeta Zeros:**
   For any target physical grid dimension \(N\) (e.g., Brillouin-zone \(k\)-meshes or non-equilibrium time contours), the minimal zero-padding overhead to the nearest 61-smooth integer \(N' = \min\{m \ge N : P^+(m) \le 61\}\) is governed asymptotically by the Dickman–de Bruijn function \(\rho(u)\), while local fluctuations in smooth-number gaps are controlled by the non-trivial zeros of the Riemann zeta function \(\zeta(s)\). Within practical supercomputing memory bounds, 61-smooth integers form an ultra-dense lattice requiring virtually zero padding waste.

---

## The Bona Sapiens Determinism Contract

This codebase adheres to a fixed set of engineering invariants. Deviations are treated as defects.

1. **`std::span` at every public boundary.** All transform entry points accept `std::span<std::complex<T>>`. Buffer extent is validated before any arithmetic is performed.
2. **True in-place kernel.** The Singleton engine uses only \(O(\mathrm{maxf} + \mathrm{maxp})\) auxiliary workspace; no \(O(N)\) data scratch buffer and no recursive stage stack.
3. **Absolute type safety.** Scalar types are constrained by the `FloatingPoint` concept (`std::floating_point<T>`). Instantiations outside this set do not compile.
4. **61-smooth domain firewall.** Lengths with a prime factor \(> 61\) are rejected at the API boundary.
5. **Numerical contract.** Double-precision round-trip accuracy is verified at \(\lVert x - \mathrm{IFFT}(\mathrm{FFT}(x)) \rVert_\infty < 10^{-12}\) across power-of-two, power-of-three, mixed-prime, and square/square-free composite lengths (including \(N \in \{61, 210, 442, 1980\}\)).

---

## Usage

```cpp
#include <complex>
#include <span>
#include <vector>

#include "SingletonFFT.hpp"

int main() {
    std::vector<std::complex<double>> signal(442);
    // ... populate signal ...

    mosaiq::fft::fft(std::span{signal});   // forward
    mosaiq::fft::ifft(std::span{signal});  // inverse, normalized by 1/N

    return 0;
}
```

Compile-time factorization is available for fixed \(N\):

```cpp
constexpr auto factors = mosaiq::fft::factorize_ct<442>();
```

---

## CMake Integration

Consume this library from another CMake project via `FetchContent`:

```cmake
cmake_minimum_required(VERSION 3.20)
project(MySpectralApp LANGUAGES CXX)

include(FetchContent)

FetchContent_Declare(
    mosaiq_singleton_fft
    GIT_REPOSITORY https://github.com/BonaSapiensInc/MosaiQ-Singleton-FFT.git
    GIT_TAG        main
)

FetchContent_MakeAvailable(mosaiq_singleton_fft)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE mosaiq::fft)
```

The exported alias `mosaiq::fft` maps to the static target `mosaiq_singleton_fft`. Include paths are propagated automatically; add `#include "SingletonFFT.hpp"` in consumer translation units.

> **Note:** AGPLv3 applies to FetchContent consumers that redistribute networked or SaaS binaries unless a commercial license is obtained.

---

## Build & Test

Requirements: CMake 3.20 or later, a C++20-compliant toolchain (GCC 11+, Clang 14+, or MSVC 19.29+).

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

The test binary `mosaiq_singleton_fft_tests` (Catch2 v3) exercises:

- Runtime and `consteval` factorization
- Forward–inverse round-trip at \(10^{-12}\) tolerance (double)
- Mixed-prime lengths (`442`), 61-smooth boundaries (`61`, `210`, `1980`)
- Rejection of unsupported primes (`67`)
- Impulse-response preservation

Strict warnings are enforced (`-Wall -Wextra -Wpedantic -Werror` on GCC/Clang; `/W4 /WX` on MSVC).

---

## Ecosystem

`MosaiQ-Singleton-FFT` serves as the foundational zero-dependency spectral primitive for the **BonaSapiens MosaiQ** computational physics suite (`MosaiQ-Wigner`, `MosaiQ-Hodge`, `MosaiQ-Faddeev`, and `MosaiQ-Keldysh`).

---

## Citation

If you use `MosaiQ-Singleton-FFT` in scientific research, cite the archived release by its Zenodo DOI. This assigns a permanent, versioned identifier to the software artifact and satisfies reproducibility requirements for peer-reviewed publication.

```bibtex
@software{mosaiq_singleton_fft,
  author       = {Kim, In-Gee},
  title        = {MosaiQ-Singleton-FFT: Modern C++20 Mixed-Radix FFT Engine},
  month        = jul,
  year         = {2026},
  publisher    = {Zenodo},
  version      = {v0.1.0},
  doi          = {10.5281/zenodo.21236967},
  url          = {https://doi.org/10.5281/zenodo.21236967}
}
```

Related historical paper:

* I. G. Kim, *J. Korean Magn. Soc.* **22**(2), 33 (2012). DOI: [10.4283/JKMS.2012.22.2.033](https://doi.org/10.4283/JKMS.2012.22.2.033)

---

## License

Dual-licensed under **GNU Affero General Public License v3.0** ([`LICENSE`](LICENSE)) and a **\$199 Commercial License** ([`LICENSE-COMMERCIAL.md`](LICENSE-COMMERCIAL.md)).

* Open-source / academic / networked SaaS: AGPLv3 copyleft applies.
* Closed-source product integration without AGPL disclosure: purchase the commercial license.
