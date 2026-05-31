---
layout: default
title: Installation
nav_order: 2
description: "Installing the GRAB binary (SPAsqr) and the external tools used by the workflows."
has_children: false
---

# **Installation**

SPA<sub>SQR</sub> ships as the `SPAsqr` method of the [**GRAB**](https://github.com/GeneticAnalysisinBiobanks/GRAB) C++ program. To install GRAB, there are mainly two options.

## Option 1 — Build from source (recommended)

GRAB is self-contained — all third-party libraries are bundled in the source tree, so the only requirements are a (reasonably new) C++17 compiler and GNU `make`. Building from source produces a binary that uses the best available SIMD instruction set on your CPU (AVX2 / AVX-512 on x86, NEON on ARM).

```bash
git clone --depth=1 https://github.com/GeneticAnalysisinBiobanks/GRAB.git
cd GRAB
make -j
build/grab2 --version
```

## Option 2 — Prebuilt binary

Download the archive for your platform from the [GRAB Releases page](https://github.com/GeneticAnalysisinBiobanks/GRAB/releases/latest). Builds are provided for `linux-x86_64`, `linux-aarch64`, `macos-arm64`, and `windows-x86_64`. 

## External tools

Our workflows call three external programs. Easiest install path for each:

- [**LDAK-KVIK**](https://dougspeed.com/ldak-kvik/) — download the prebuilt binary from the LDAK site.
- [**REGENIE**](https://rgcgithub.github.io/regenie/) — `conda install -c conda-forge -c bioconda regenie`.
- [**PLINK 2**](https://www.cog-genomics.org/plink/2.0/) — download the prebuilt binary from the PLINK 2 site.
