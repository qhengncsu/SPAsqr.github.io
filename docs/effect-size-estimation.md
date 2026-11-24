---
layout: default
title: Effect-size estimation
nav_order: 5
description: "Per-marker, per-tau γ̂(τ), SE, and a heterogeneity test via SPAsqr Wald mode."
has_children: false
---

# **Effect-size estimation: `--spasqr-mode wald`**

Workflows 1 and 2 use **score mode** to screen the whole genome. For each chromosome and each quantile $\tau$, GRAB fits a null smoothed quantile regression model without the variant. It uses that model to test $H_0: \gamma(\tau) = 0$ for every variant on the chromosome. Reusing the null model makes it practical to test millions of variants, and the saddlepoint approximation keeps rare-variant $p$-values calibrated. The score test gives a $p$-value but does not estimate the effect size.

Use **Wald mode** (`--spasqr-mode wald`) to estimate effects at a short list of variants, such as genome-wide significant hits from score mode or variants of particular interest. For each variant and quantile, GRAB fits the full model with the variant as a covariate. It reports the estimated effect $\hat\gamma(\tau)$, its standard error, and a Wald $p$-value. It also tests whether the effect is the same across quantiles. These additional model fits make a genome-wide run computationally infeasible.

## Complete pipeline

This example uses the phenotype file after RINT (`simu_geno_int.txt`) and the prediction list (`simu_geno_ldak_pred.list`) from [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html).

```bash
# 1. List the variants to estimate, one ID per line
cat > simu_geno_wald_extract <<EOF
SNP_1031
SNP_1040
SNP_1428
SNP_187
SNP_2170
SNP_2287
SNP_3240
SNP_4380
EOF

# 2. Run SPAsqr in Wald mode
./grab2 --method SPAsqr --spasqr-mode wald --spasqr-h-scale 5 \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --extract simu_geno_wald_extract \
    --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
    --threads 8 \
    --out spasqr_effect
```

## Step by step

### 1. Variant list

Create a text file with one variant ID per line. The IDs must match the `ID` column in the `.bim` file. The example list is available in [`data/simu_geno_wald_extract`](https://github.com/qhengncsu/SPAsqr.github.io/tree/main/data).

### 2. Wald-mode run

Start with the command from Workflow 1. Use `--spasqr-mode wald` to estimate effect sizes and `--extract` to select the variants in your list. We use a narrower bandwidth (`--spasqr-h-scale 5`) to reduce smoothing bias. The other options, including `--spasqr-taus`, `--pred-list`, and `--sp-grm-*`, work the same way as in score mode.

## Output format

GRAB writes one file per trait. Here the files are `spasqr_effect.Quantitative1.SPAsqr` and `spasqr_effect.Quantitative2.SPAsqr`. With nine quantiles, the default output has 41 columns. The first 10 are fixed. The next three sets of nine give the effects, standard errors, and Wald $p$-values at each quantile. The final four columns report the results of the heterogeneity test.

```
CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  LOG10P_HWE  LOG10P_CCT
LOG10P_tau0.1  ...  LOG10P_tau0.9
BETA_tau0.1    ...  BETA_tau0.9
SE_tau0.1      ...  SE_tau0.9
BETA_HOM  SE_HOM  HET_Q  LOG10P_HET
```

As in score mode, `--output-z` adds the optional `Z_tau` block, for 50 columns in total.

| Columns | Meaning |
| --- | --- |
| `BETA_tau<τ>` | Effect size at quantile $\tau$, on the scale of the trait after GRAB's transformation (the RINT scale in this example). |
| `SE_tau<τ>` | Standard error of `BETA_tau<τ>`. |
| `LOG10P_tau<τ>` | Two-sided Wald $p$-value at quantile $\tau$, reported as $-\log_{10}p$. `LOG10P_CCT` combines them across quantiles. |
| `BETA_HOM`, `SE_HOM` | Common effect estimate and its standard error under the hypothesis that the effect is the same at all quantiles. |
| `HET_Q`, `LOG10P_HET` | Wald statistic for that hypothesis and its $p$-value, reported as $-\log_{10}p$. With $K$ quantiles, `HET_Q` follows a $\chi^2$ distribution with $K-1$ degrees of freedom when the effect is the same at all quantiles. |
| `Z_tau<τ>` | Written with `--output-z`. Wald $Z$-score. |

The heterogeneity test assumes that the effect estimates are approximately normally distributed. For rare variants this approximation can give poorly calibrated $p$-values. GRAB therefore runs the test only for variants with MAF at or above `--spasqr-het-min-maf` (default 0.01). Below that threshold, GRAB reports `NA` for `BETA_HOM`, `SE_HOM`, `HET_Q`, and `LOG10P_HET`. It still reports the per-quantile results as usual.

To see how the estimated effect varies across quantiles, plot `BETA_tau` against the quantile level and draw error bars at ±1.96 × `SE_tau` around each estimate.

## Other options

| Flag | Default | What it does |
| --- | --- | --- |
| `--spasqr-taus` | `0.1,0.3,0.5,0.7,0.9` | Quantile levels, comma-separated (max 20). |
| `--spasqr-h-scale` | `5` | Bandwidth divisor, `h = IQR / k`. Wald mode defaults to a narrower bandwidth than score mode (`5` vs `3`) to reduce smoothing bias. |
| `--spasqr-het-min-maf` | `0.01` | Minimum MAF for the heterogeneity test. `0` runs the test for every variant. |
| `--output-z` | off | Also write the `Z_tau` columns. |
