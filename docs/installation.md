---
layout: default
title: Installation
nav_order: 2
description: "Installing the GRAB binary (SPAsqr) and the external tools used by the workflows."
has_children: false
---

# **Installation**

SPA<sub>SQR</sub> is available as the `SPAsqr` method in [**GRAB**](https://github.com/GeneticAnalysisinBiobanks/GRAB). You can build GRAB from source or download the executable for your platform.

## Option 1 — Build from source (recommended)

To build GRAB, you need a compiler that supports C++17 and GNU `make`. The source includes all third-party libraries. Building on your own machine lets GRAB use the best SIMD instruction set your CPU supports (AVX2 or AVX-512 on x86, NEON on ARM).

```bash
git clone --depth=1 https://github.com/GeneticAnalysisinBiobanks/GRAB.git
cd GRAB
make -j
build/grab2 --version
```

## Option 2 — Prebuilt binary

Download the archive for your platform from the [GRAB releases page](https://github.com/GeneticAnalysisinBiobanks/GRAB/releases/latest). Prebuilt binaries are available for `linux-x86_64`, `linux-aarch64`, `macos-arm64`, and `windows-x86_64`.

## External tools

The tools you need depend on your workflow. Use LDAK-KVIK or REGENIE to train the LOCO PGS. Workflow 2 also uses PLINK 2 to build a sparse GRM.

- [**LDAK-KVIK**](https://dougspeed.com/ldak-kvik/) — download the prebuilt binary from the LDAK site.
- [**REGENIE**](https://rgcgithub.github.io/regenie/) — `conda install -c conda-forge -c bioconda regenie`.
- [**PLINK 2**](https://www.cog-genomics.org/plink/2.0/) — download the prebuilt binary from the PLINK 2 site.
