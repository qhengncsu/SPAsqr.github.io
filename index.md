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
[**GRAB**](https://github.com/qiang-heng/GRAB-feat-cpp) command-line
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

1. **Step 0a — LOCO polygenic scores.** Train chromosome-specific PGS
   on the phenotype using
   [LDAK-KVIK](https://dougspeed.com/kvik/) (recommended) or
   [REGENIE](https://rgcgithub.github.io/regenie/). Optional — but
   strongly recommended for power.

2. **Step 0b — Sparse GRM.** Construct a sparse genetic relationship
   matrix using
   [PLINK 2](https://www.cog-genomics.org/plink/2.0/) (preferred since
   late 2025) or GCTA. Optional — supply when the cohort has
   non-trivial relatedness.

3. **Steps 1–2 — Run SPA<sub>SQR</sub>.** A single
   `grab --method SPAsqr` call reads the phenotype + covariates,
   fits the null SQR model per chromosome with bandwidth
   $h = \mathrm{IQR}(\tilde Y - \hat Y_{-c}) / k$, applies SPA in the
   tails, and writes one tab-delimited result file per phenotype.

4. **(Optional) Effect-size estimation.** Re-run with
   `--spasqr-mode wald` to obtain per-marker per-$\tau$
   $\hat\beta_G$ and SE via M-estimation sandwich variance.

## Where to go next

- [Installation]({{ site.baseurl }}/docs/installation.html) — building the GRAB binary.
- [Step 0a — LOCO polygenic scores]({{ site.baseurl }}/docs/step-0a-loco-pgs.html)
- [Step 0b — Sparse GRM]({{ site.baseurl }}/docs/step-0b-sparse-grm.html)
- [Running SPA<sub>SQR</sub>]({{ site.baseurl }}/docs/running-spasqr.html) — the main usage page.
- [Effect-size estimation]({{ site.baseurl }}/docs/effect-size-estimation.html) — Wald mode.
- [Strategies for improving statistical power]({{ site.baseurl }}/docs/strategies.html)
- [FAQ]({{ site.baseurl }}/docs/faq.html)

## Citation

If you use SPA<sub>SQR</sub> in your work, please cite the
SPA<sub>SQR</sub> manuscript (link to be added on publication).
