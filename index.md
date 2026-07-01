---
layout: default
title: Overview
nav_order: 1
description: "SPAsqr: saddlepoint-approximated smoothed quantile regression for biobank-scale GWAS."
permalink: /
---

# **SPA<sub>SQR</sub>**

SPA<sub>SQR</sub> is a fast, robust framework for **quantile regression
GWAS** on quantitative traits at biobank scale. It is implemented as the
`SPAsqr` method of the
[**GRAB**](https://github.com/GeneticAnalysisinBiobanks/GRAB) command-line
binary — a single C++ executable that runs on Linux, macOS, and Windows.

## What is quantile regression GWAS?

- **Linear GWAS** models the **mean** of a trait — one average effect per
  variant.
- **Quantile regression (QR)** models the trait's **quantiles** — so a
  variant's effect can change across the phenotype distribution.
- This reveals **quantile-dependent effects** (subgroup-specific, heterogeneous effects) that mean-based GWAS
  dilutes or misses entirely.

## Why SPA<sub>SQR</sub>?

- **Detects effects linear GWAS misses** — tests many quantiles at once
  and combines them into one $p$-value via the Cauchy combination test.
- **Calibrated for rare variants** — saddlepoint approximation keeps
  $p$-values accurate for rare variants.
- **Handles relatedness** — uses LOCO polygenic scores as an offset, with
  an optional sparse genetic relationship matrix for further variance calibration.
- **More power from smoothing** — convolution-smoothed QR lowers
  variance, offering greater association power than non-smooth QR.
- **Biobank-scale speed** — multithreaded C++ implementation that rivals REGENIE and
  LDAK-KVIK in efficiency.

![SPASQR association results for glucose in UK Biobank](docs/assets/images/glu_combined_leadb.png)

## Pipeline

1. **Workflow 1 — LOCO PGS + SPA<sub>SQR</sub>.** Train 
   LOCO polygenic scores with
   [LDAK-KVIK](https://dougspeed.com/ldak-kvik/) or
   [REGENIE](https://rgcgithub.github.io/regenie/), then run
   `grab --method SPAsqr` with `--pred-list`. The recommended path for
   cohorts with low relatedness.

2. **Workflow 2 — LOCO PGS + GRM + SPA<sub>SQR</sub>.** In addition to the
   LOCO PGS, build a sparse genetic relationship matrix using
   [PLINK 2](https://www.cog-genomics.org/plink/2.0/), and pass it via `--sp-grm-plink2` so that the
   score-statistic variance takes relatedness into account. Recommended for
   cohorts with strong relatedness.

3. **(Optional) Effect-size estimation.** Re-run with
   `--spasqr-mode wald` to obtain effect size estimates and their standard errors.

## Where to go next

- [Installation]({{ site.baseurl }}/docs/installation.html) 
- [Workflow 1: LOCO PGS + SPA<sub>SQR</sub>]({{ site.baseurl }}/docs/workflow-1.html)
- [Workflow 2: LOCO PGS + GRM + SPA<sub>SQR</sub>]({{ site.baseurl }}/docs/workflow-2.html)
- [LDAK-KVIK runs faster one-trait-at-a-time]({{ site.baseurl }}/docs/all-vs-per-trait.html) 
- [Effect-size estimation]({{ site.baseurl }}/docs/effect-size-estimation.html) 
- [FAQ]({{ site.baseurl }}/docs/faq.html)
