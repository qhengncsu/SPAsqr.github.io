---
layout: default
title: Overview
nav_order: 1
description: "SPAsqr: saddlepoint-approximated smoothed quantile regression for biobank-scale GWAS."
permalink: /
---

# **Documentation for SPA<sub>SQR</sub>**

SPA<sub>SQR</sub> is a saddlepoint-approximated, **smoothed quantile
regression** (SQR) framework for genome-wide association studies on
quantitative traits. It performs association testing across multiple
quantile levels $\tau$ simultaneously, combines them via the Cauchy
combination test (CCT), accommodates **leave-one-chromosome-out (LOCO)
polygenic scores** as offsets and an optional **sparse GRM** for
variance calibration, and applies **saddlepoint approximation** in the
tails so that rare variants are well-calibrated even under heavy-tailed
or skewed phenotypes.

SPA<sub>SQR</sub> is implemented as the `SPAsqr` method of the
[**GRAB**](https://github.com/qhengncsu/GRAB-feat-cpp) command-line
binary — a single statically linked C++17 executable that runs
unmodified on Linux, macOS, and Windows.

## Why smoothed quantile regression?

Conventional linear GWAS targets the **conditional mean** of $Y$, which
loses power when the genetic effect concentrates in the tails of the
phenotype distribution, when $Y$ is non-Gaussian, or when the effect is
quantile-dependent (e.g. heteroskedastic / dispersion effects).
SPA<sub>SQR</sub> targets the **conditional quantiles** $Q_\tau(Y \mid
G, X)$ at a user-specified grid of $\tau$ levels (default
$\{0.1, 0.3, 0.5, 0.7, 0.9\}$) and combines them into a single
$p$-value via CCT. Smoothing the check loss with a Gaussian kernel of
bandwidth $h$ lowers the variance of the rank score, which translates
to a smaller denominator in the score statistic and hence higher power
than non-smooth quantile regression.

## What SPA<sub>SQR</sub> does, in one paragraph

For each chromosome $c$, SPA<sub>SQR</sub> fits a **null SQR model**
$Q_\tau(Y \mid X, \hat Y_{-c}) = X^\top \beta_\tau + \hat Y_{-c}$ on
the chromosome-specific LOCO PGS as an offset, then computes a
**variance-stabilized score statistic** $S_j = G_j^\top R$ for every
variant $j$ on that chromosome. The variance of $S_j$ is
$\widehat\sigma_g^{\,2}(G_j)\, R^\top \Phi\, R$ where $\Phi$ is the
sparse GRM (or $I_n$ for unrelated samples). A saddlepoint
approximation is applied whenever $|S_j|$ exceeds a configurable
$z$-threshold so that tail $p$-values stay calibrated for rare or
unbalanced traits. Per-$\tau$ $p$-values are combined into a single
**P_CCT** via the Cauchy combination test.

## Pipeline

1. **Workflow 1 — LOCO PGS + SPA<sub>SQR</sub>.** Train chromosome-specific
   LOCO polygenic scores with
   [LDAK-KVIK](https://dougspeed.com/ldak-kvik/) or
   [REGENIE](https://rgcgithub.github.io/regenie/), then run
   `grab --method SPAsqr` with `--pred-list`. The recommended path for
   essentially unrelated cohorts.

2. **Workflow 2 — LOCO PGS + GRM + SPA<sub>SQR</sub>.** In addition to the
   LOCO PGS, build a sparse genetic relationship matrix using
   [PLINK 2](https://www.cog-genomics.org/plink/2.0/) (preferred since
   late 2025) or GCTA, and pass it via `--sp-grm-plink2` so that the
   score-statistic variance is GRM-aware. Recommended whenever the
   cohort retains first- or second-degree relatives.

3. **(Optional) Effect-size estimation.** Re-run with
   `--spasqr-mode wald` to obtain per-marker per-$\tau$
   $\hat\beta_G$ and SE via M-estimation sandwich variance.

## Where to go next

- [Installation]({{ site.baseurl }}/docs/installation.html) — building the GRAB binary.
- [Workflow 1: LOCO PGS + SPA<sub>SQR</sub>]({{ site.baseurl }}/docs/workflow-1.html)
- [Workflow 2: LOCO PGS + GRM + SPA<sub>SQR</sub>]({{ site.baseurl }}/docs/workflow-2.html)
- [All traits at once vs one trait at a time]({{ site.baseurl }}/docs/all-vs-per-trait.html) — scaling SPA<sub>SQR</sub> across many phenotypes.
- [Effect-size estimation]({{ site.baseurl }}/docs/effect-size-estimation.html) — Wald mode.
- [Strategies for improving statistical power]({{ site.baseurl }}/docs/strategies.html)
- [FAQ]({{ site.baseurl }}/docs/faq.html)

## Citation

If you use SPA<sub>SQR</sub> in your work, please cite the
SPA<sub>SQR</sub> manuscript (link to be added on publication).
