# Legacy Reference Archive

This directory holds historical provenance material for the MosaiQ-Singleton-FFT lineage.

## Archived Material

| File | Description |
|------|-------------|
| `JKoreanMagSoc_22_33_(2012).pdf` | In-Gee Kim, *J. Korean Magn. Soc.* 22(2), 33 (2012) — Singleton mixed-radix FFT with dynamic memory allocation (migrated from `../legacy_singleton/`) |
| `cfft.F` | Legacy complex FFT driver (Fortran) |
| `primfact.F` | Legacy prime factorization routine (Fortran) |

## Migration Instructions

Legacy material from `../legacy_singleton/` has been migrated. Additional documents may be added as follows:

| Source | Suggested Destination |
|--------|----------------------|
| Singleton (1979) IEEE DSP Program paper | `SINGLETON_1979_MIXED_RADIX_FFT.pdf` |
| Original Fortran Y / F90 translations | `legacy/fortran/` |
| In-Gee Kim Singleton-FFT repository notes | `legacy/ingeekim-singleton-fft/` |

## Policy

- Verify copyright and redistribution rights before committing PDFs.
- Prefer immutable filenames with year prefix for traceability.
- Do not store credentials, build artifacts, or binary object files here.

See `docs/02_SINGLETON_CFFT_STRUCTURE.md` §8 for the architectural context.
