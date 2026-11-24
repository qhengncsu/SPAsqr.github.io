---
layout: default
title: Overview
nav_order: 1
description: "SPAsqr: saddlepoint-approximated smoothed quantile regression for biobank-scale GWAS."
permalink: /
---

# **SPA<sub>SQR</sub>**

SPA<sub>SQR</sub> uses smoothed quantile regression to test genetic variants for association with quantitative traits in biobank-scale studies. It is fast and robust, and runs as the `SPAsqr` method in [**GRAB**](https://github.com/GeneticAnalysisinBiobanks/GRAB), a single C++ executable for Linux, macOS, and Windows.

## What is quantile regression GWAS?

A linear GWAS estimates the average effect of each variant on a trait. Quantile regression (QR) estimates its effect at different points in the trait distribution. For example, a variant may have a stronger effect at the upper quantiles than near the median. Looking across quantiles can reveal heterogeneous or subgroup-specific effects that are hard to detect from the mean alone. QR GWAS results are also generally less sensitive than linear GWAS results to monotonic transformations of the trait, such as rank-based inverse normal transformation.

## Why SPA<sub>SQR</sub>?

- **Association across quantiles.** SPA<sub>SQR</sub> tests for association at several quantiles and combines their $p$-values into one $p$-value per variant using the Cauchy combination test.
- **Rare variants.** A saddlepoint approximation keeps association $p$-values calibrated for rare variants.
- **Relatedness.** Leave-one-chromosome-out (LOCO) polygenic scores (PGS) serve as an offset. An optional sparse genetic relationship matrix (GRM) provides further variance calibration.
- **Power.** Convolution smoothing reduces variance and improves power compared with unsmoothed quantile regression.
- **Speed.** The multithreaded C++ implementation is comparable in efficiency to REGENIE and LDAK-KVIK.

![SPASQR association results for glucose in UK Biobank](docs/assets/images/glu_combined_leadb.png)

## Pipeline

1. **Workflow 1 — LOCO PGS + SPA<sub>SQR</sub>.** For cohorts with little relatedness, we recommend training the LOCO PGS with [LDAK-KVIK](https://dougspeed.com/ldak-kvik/) or [REGENIE](https://rgcgithub.github.io/regenie/), then passing the prediction list to `grab --method SPAsqr` with `--pred-list`.

2. **Workflow 2 — LOCO PGS + GRM + SPA<sub>SQR</sub>.** For cohorts with strong relatedness, we recommend adding a sparse GRM to account for relatedness when estimating the variance of the score statistic. Build the matrix with [PLINK 2](https://www.cog-genomics.org/plink/2.0/) and pass it to GRAB with `--sp-grm-plink2`. For admixed cohorts, we recommend using the LOCO PGS alone because the GRM is confounded with population structure.

3. **Effect-size estimation (optional).** For selected variants, run GRAB with `--spasqr-mode wald` to estimate the effect and its standard error at each quantile and to test whether the effect varies across quantiles.

## Where to go next

- [Installation]({{ site.baseurl }}/docs/installation.html) 
- [Workflow 1: LOCO PGS + SPA<sub>SQR</sub>]({{ site.baseurl }}/docs/workflow-1.html)
- [Workflow 2: LOCO PGS + GRM + SPA<sub>SQR</sub>]({{ site.baseurl }}/docs/workflow-2.html)
- [Effect-size estimation]({{ site.baseurl }}/docs/effect-size-estimation.html) 
- [FAQ]({{ site.baseurl }}/docs/faq.html)

## Getting help
{: #ai-agents-and-llms-can-help-you }

If you run into trouble installing GRAB or running an example, try sharing the tutorial URL and the error message with Claude Code or Codex. When we tested them, each tool was able to install GRAB and run every example using just the URL.

We welcome bug reports and suggestions at heng.qiang6@gmail.com.
