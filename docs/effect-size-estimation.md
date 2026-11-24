---
layout: default
title: Effect-size estimation
nav_order: 5
description: "Per-marker, per-tau γ̂(τ), SE, and a heterogeneity test via SPAsqr Wald mode."
has_children: false
---

# **Effect-size estimation: `--spasqr-mode wald`**

SPA<sub>SQR</sub>'s default mode, used in Workflows 1 and 2, is **genome-wide screening**: for each chromosome and each quantile level $\tau$ it fits one null smoothed quantile regression model without the variant, and then tests $H_0: \gamma(\tau) = 0$ for every variant on that chromosome with a score test. Because the null model does not depend on the variant, this is fast enough for millions of variants, and the saddlepoint approximation keeps the $p$-values calibrated for rare variants. However, a score test does not produce an estimate of the effect size $\gamma(\tau)$.

To obtain effect sizes, SPA<sub>SQR</sub> provides a **Wald mode** (`--spasqr-mode wald`). For each variant and each $\tau$, it fits the full smoothed quantile regression model including the variant, and reports the estimated effect $\hat\gamma(\tau)$ together with its standard error and a Wald $p$-value. It also tests whether the effect is the same at all quantiles (a heterogeneity test). Since this requires one model fit per variant, Wald mode is intended for a short list of variants rather than the whole genome, typically the genome-wide-significant hits from score mode or a curated set of variants of interest.

## Complete pipeline

Here we assume that the RINT-transformed phenotype file `simu_geno_int.txt` and the prediction list `simu_geno_ldak_pred.list` from [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html) are already in place.

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

We first write a plain-text file with one variant ID per line, matching the `ID` column of the `.bim` file. The example file is available at [`data/simu_geno_wald_extract`](https://github.com/qhengncsu/SPAsqr.github.io/tree/main/data).

### 2. Wald-mode run

We then rerun GRAB with three changes from the Workflow 1 command: `--spasqr-mode wald` switches to effect-size estimation, `--spasqr-h-scale 5` sets a narrower bandwidth to reduce smoothing bias, and `--extract` restricts the run to the listed variants. All other flags, including `--spasqr-taus`, `--pred-list`, and `--sp-grm-*`, work exactly as in score mode.

## Output format

GRAB writes one file per trait, `spasqr_effect.Quantitative1.SPAsqr` and `spasqr_effect.Quantitative2.SPAsqr`. With nine quantiles the default output has 41 columns: 10 fixed columns, three blocks of nine per-quantile columns, and four columns for the heterogeneity test.

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

The heterogeneity test relies on the asymptotic normality of the effect estimates and may not be well calibrated for rare variants, so it is run only for variants with MAF at least `--spasqr-het-min-maf` (default 0.01). Below this threshold, `BETA_HOM`, `SE_HOM`, `HET_Q`, and `LOG10P_HET` are `NA`; the per-quantile results are unaffected.

Plotting `BETA_tau` against the quantile level with error bars of ±1.96 × `SE_tau` gives the quantile-effect profile of a variant.

## Other options

| Flag | Default | What it does |
| --- | --- | --- |
| `--spasqr-taus` | `0.1,0.3,0.5,0.7,0.9` | Quantile levels, comma-separated (max 20). |
| `--spasqr-h-scale` | `5` | Bandwidth divisor, `h = IQR / k`. Wald mode defaults to a narrower bandwidth than score mode (`5` vs `3`) to reduce smoothing bias. |
| `--spasqr-het-min-maf` | `0.01` | Minimum MAF for the heterogeneity test. `0` runs the test for every variant. |
| `--output-z` | off | Also write the `Z_tau` columns. |
