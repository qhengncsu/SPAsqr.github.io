---
layout: default
title: Effect-size estimation
nav_order: 5
description: "Per-marker, per-tau γ̂(τ) and SE via SPAsqr Wald mode."
has_children: false
---

# **Effect-size estimation — `--spasqr-mode wald`**

SPA<sub>SQR</sub>'s default mode is **genome-wide screening**: it tests $H_0: \gamma(\tau) = 0$ at every $\tau$ with a score test, fitting one null model per chromosome and quantile. It gives calibrated $p$-values and signed $Z$-scores per marker, but **not** effect-size estimates. To estimate effect sizes for a list of SNPs, use **Wald mode**.

## Example

Prepare a plain-text file (`simu_geno_wald_extract`) with one variant ID per line — typically the genome-wide-significant hits from score mode, or a curated list. An 8-variant example is at [`data/simu_geno_wald_extract`](https://github.com/qhengncsu/SPAsqr.github.io/tree/main/data):

```
$ cat simu_geno_wald_extract
SNP_1031
SNP_1040
SNP_1428
SNP_187
SNP_2170
SNP_2287
SNP_3240
SNP_4380
```

We may estimate the effect sizes for those variants on `Quantitative1,Quantitative2` via:

```bash
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

The outputs `spasqr_effect.Quantitative1.SPAsqr` and `spasqr_effect.Quantitative2.SPAsqr`. With nine quantiles each file has 46 columns:

```
CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  P_CCT
P_tau0.1    ...  P_tau0.9
Z_tau0.1    ...  Z_tau0.9
BETA_tau0.1 ...  BETA_tau0.9
SE_tau0.1   ...  SE_tau0.9
```

For each requested $\tau$:

- `BETA_tau<val>` — the estimated effect size at quantile $\tau$, on the transformed scale (e.g. INT scale under `--pheno-transform int`).
- `SE_tau<val>` — standard error of `BETA_tau<val>`.
- `Z_tau<val>` and `P_tau<val>` — Wald Z-score and two-sided $p$-values obtained via normal approximation.


Wald mode defaults to a narrower bandwidth than score mode (`--spasqr-h-scale 5` vs `3`) to keep smoothing bias low at the cost of slightly slower convergence; override via `--spasqr-h-scale <k>` when desired.
