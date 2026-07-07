# MosaiQ-Singleton-FFT

![C++20](https://img.shields.io/badge/C++-20-blue.svg)
![CMake](https://img.shields.io/badge/CMake-3.20+-064F8C.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

A modernized, strictly deterministic C++20 implementation of the mixed-radix Cooley–Tukey fast Fourier transform. This library is the Bona Sapiens reference engine for in-place spectral analysis: factorization, butterfly stages, and prime-radix kernels are implemented with explicit compile-time structure and zero tolerance for undefined behavior at API boundaries.

---

## Overview

`MosaiQ-Singleton-FFT` decomposes a transform length \(N\) into small radices and executes a recursive decimation-in-frequency schedule. Unlike legacy Singleton implementations limited to hand-coded radices 2–5, this engine **explicitly unrolls prime radix butterflies through 61** via template metaprogramming (`PrimeRadix<T, P>`), yielding fully specialized \(O(P^2)\) micro-kernels at each stage without runtime dispatch overhead in the hot path.

Supported radices:

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
| Transform engine | `src/SingletonFFT.cpp` |
| Prime radix expansion | `src/detail/PrimeRadix.hpp` |
| Validation suite | `tests/test_fft.cpp` |

---

## The Bona Sapiens Determinism Contract

This codebase adheres to a fixed set of engineering invariants. Deviations are treated as defects.

1. **`std::span` at every boundary.** All transform entry points accept `std::span<std::complex<T>>`. Buffer extent is validated before any arithmetic is performed.
2. **Zero raw pointers in the public and butterfly layers.** Radix kernels operate on `std::span<std::complex<T>, P>`. No pointer arithmetic is exposed to callers.
3. **Absolute type safety.** Scalar types are constrained by the `FloatingPoint` concept (`std::floating_point<T>`). Instantiations outside this set do not compile.
4. **Numerical contract.** Double-precision round-trip accuracy is verified at \(\lVert x - \mathrm{IFFT}(\mathrm{FFT}(x)) \rVert_\infty < 10^{-12}\) across power-of-two, power-of-three, and mixed-prime lengths (including \(N = 442 = 2 \times 13 \times 17\)).

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
    GIT_REPOSITORY https://github.com/Bona-Sapiens/MosaiQ-Singleton-FFT.git
    GIT_TAG        main
)

FetchContent_MakeAvailable(mosaiq_singleton_fft)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE mosaiq::fft)
```

The exported alias `mosaiq::fft` maps to the static target `mosaiq_singleton_fft`. Include paths are propagated automatically; add `#include "SingletonFFT.hpp"` in consumer translation units.

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
- Mixed-prime lengths that engage radices 13 and 17
- Impulse-response preservation

Strict warnings are enforced (`-Wall -Wextra -Wpedantic -Werror` on GCC/Clang; `/W4 /WX` on MSVC).

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
