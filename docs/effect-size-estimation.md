---
layout: default
title: Effect-size estimation
nav_order: 5
description: "Per-marker, per-tau γ̂(τ) and SE via SPAsqr Wald mode."
has_children: false
---

# **Effect-size estimation: `--spasqr-mode wald`**

Score mode (Workflows 1–2) gives calibrated $p$-values and $Z$-scores but no effect sizes. **Wald mode** fits the smoothed QR model per variant and reports $\hat\gamma(\tau)$ with its standard error. Use it on a short list of variants, typically the genome-wide-significant hits from score mode.

## Complete pipeline

This pipeline assumes the RINT-transformed phenotypes and the prediction list from [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html) are already in place.

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
./grab2 --method SPAsqr --spasqr-mode wald \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --extract simu_geno_wald_extract \
    --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
    --pheno-transform int \
    --threads 8 \
    --out spasqr_effect
```

## Step by step

### 1. Variant list

We first write a plain-text file with one variant ID per line, matching the `ID` column of the `.bim` file. The example file is available at [`data/simu_geno_wald_extract`](https://github.com/qhengncsu/SPAsqr.github.io/tree/main/data).

### 2. Wald-mode run

We then rerun GRAB with two changes from the Workflow 1 command: `--spasqr-mode wald` switches to effect-size estimation, and `--extract` restricts the run to the listed variants. All other flags, including `--spasqr-taus`, `--pheno-transform`, `--pred-list`, and `--sp-grm-*`, work exactly as in score mode.

## Output format

GRAB writes one file per trait, `spasqr_effect.Quantitative1.SPAsqr` and `spasqr_effect.Quantitative2.SPAsqr`. With nine quantiles each file has 46 columns: 10 fixed columns followed by four blocks of nine per-quantile columns.

```
CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  P_CCT
P_tau0.1     ...  P_tau0.9
Z_tau0.1     ...  Z_tau0.9
BETA_tau0.1  ...  BETA_tau0.9
SE_tau0.1    ...  SE_tau0.9
```

| Columns | Meaning |
| --- | --- |
| `BETA_tau<τ>` | Effect size at quantile $\tau$, on the transformed scale (RINT scale under `--pheno-transform int`). |
| `SE_tau<τ>` | Standard error of `BETA_tau<τ>`. |
| `Z_tau<τ>`, `P_tau<τ>` | Wald $Z$-score and two-sided normal-approximation $p$-value. |

## Other options

| Flag | Default | What it does |
| --- | --- | --- |
| `--spasqr-taus` | `0.1,0.3,0.5,0.7,0.9` | Quantile levels, comma-separated (max 20). |
| `--spasqr-h-scale` | `5` | Bandwidth divisor, `h = IQR / k`. Wald mode defaults to a narrower bandwidth than score mode (`5` vs `3`) to reduce smoothing bias. |
