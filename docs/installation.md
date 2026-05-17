---
layout: default
title: Installation
nav_order: 2
description: "Building the GRAB binary that implements SPAsqr."
has_children: false
---

# **Installation**

SPA<sub>SQR</sub> ships as the `SPAsqr` method of the
[**GRAB**](https://github.com/qhengncsu/GRAB-feat-cpp) command-line
binary. GRAB is a pure C++17 application; **all third-party
dependencies are vendored** in `third_party/` (Eigen, Boost
headers-only subset, zlib, zstd, libdeflate, plink2/pgenlib, bgen,
htslib). The only requirements on your system are a C++17 compiler
(g++ / clang++ / MinGW g++) and `make`. The build pulls nothing from
the system (no `apt-get`, `brew`, `vcpkg`, Conan, …) and produces a
single statically linked executable that runs unmodified on **Linux,
macOS, and Windows** (MSYS2 / MinGW).

## Quick install

```bash
git clone https://github.com/qhengncsu/GRAB-feat-cpp.git
cd GRAB-feat-cpp
make -j$(nproc)             # or sysctl -n hw.ncpu on macOS
./build/grab --help SPAsqr  # smoke-test
```

The build produces `build/grab` (or `build/grab.exe` on Windows). Copy
it to anywhere on your `$PATH` and you're done — there is no install
step, no shared library, and no headers exposed to users.

## Build matrix

| Platform | x86_64                  | arm64 (Apple Silicon, ARM Linux) |
| -------- | ----------------------- | -------------------------------- |
| Linux    | AVX-512 / AVX2 / scalar | NEON via Eigen + scalar          |
| macOS    | AVX-512 / AVX2 / scalar | NEON via Eigen + scalar          |
| Windows  | AVX-512 / AVX2 / scalar | (untested)                       |

The `Makefile` auto-detects platform (`uname -s`), architecture
(`uname -m`), and AVX2 availability. `GRAB_MARCH=-march=native` is the
default for best SIMD; override with `GRAB_MARCH=-march=x86-64-v2` for
a portable distribution binary.

## Verifying the build

```bash
./build/grab --help
./build/grab --help SPAsqr
./build/grab --version   # build banner + linked-library versions
```

If `--help SPAsqr` prints the SPA<sub>SQR</sub> help screen, you're ready to move
on to [Workflow 1: LOCO PGS + SPA<sub>SQR</sub>]({{ site.baseurl }}/docs/workflow-1.html).

## External tools used by the pipeline

GRAB itself has no runtime dependencies, but the surrounding pipeline
calls out to two external programs that you install separately:

- [**LDAK-KVIK**](https://dougspeed.com/ldak-kvik/) (or
  [REGENIE](https://rgcgithub.github.io/regenie/)) — for LOCO
  polygenic scores in Workflow 1.
- [**PLINK 2**](https://www.cog-genomics.org/plink/2.0/) — for sparse
  GRM construction (used in Workflow 2). Required only if you want a
  GRM-aware variance correction.

Both are downloaded as single binaries; no source build needed.

> **Note**
> If you already have a working SPA<sub>SQR</sub> binary (e.g., installed by a
> collaborator), check the version with `grab --version` to confirm it
> includes the `SPAsqr` method and `--spasqr-mode wald` support — both
> are required for the full workflow documented on this site.
