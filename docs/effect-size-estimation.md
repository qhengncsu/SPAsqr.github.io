---
layout: default
title: Effect-size estimation
nav_order: 6
description: "Per-marker, per-tau γ̂_τ and SE via SPAsqr Wald mode."
has_children: false
---

# **Effect-size estimation — `--spasqr-mode wald`**

The default `--spasqr-mode score` is optimized for **genome-wide screening**: it tests $H_0: \gamma_\tau = 0$ at every $\tau$ via a rank-score statistic computed from a single null-model fit per chromosome, and returns calibrated $p$-values and signed $Z$-scores per marker — but **not** the per-allele effect estimate $\hat\gamma_\tau$ itself. For follow-up analyses on a small candidate set (top GWAS hits, prior-literature variants), switch to **Wald mode**. Wald mode re-fits the *full* model

$$
Q_\tau(Y \mid X, G) \;=\; X^\top \beta_\tau \;+\; G\, \gamma_{\tau}
$$

once **per (marker, $\tau$)** and reports $\hat\gamma_{\tau}$, its sandwich-variance standard error, the Wald statistic, and the two-sided $p$-value. Wald mode is meant for a short, curated list of SNPs passed via `--extract`, not for a full genome-wide scan.

## Example

The user prepares a plain-text file (here `simu_geno_wald_extract`) listing one variant ID per line — typically the genome-wide-significant hits from the score-mode output of [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html), or a curated set from prior literature. This documentation repository provides an 8-variant example list at [`data/simu_geno_wald_extract`](https://github.com/qhengncsu/SPAsqr.github.io/tree/main/data):

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

We then estimate the per-allele effect sizes for those variants on both traits via:

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

The outputs `spasqr_effect.Quantitative1.SPAsqr` and `spasqr_effect.Quantitative2.SPAsqr` share the score-mode output format, with two additional columns per $\tau$:

```
CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  P_CCT
P_tau0.1    ...  P_tau0.9
Z_tau0.1    ...  Z_tau0.9
BETA_tau0.1 ...  BETA_tau0.9
SE_tau0.1   ...  SE_tau0.9
```

For each requested $\tau$:

- `BETA_tau<val>` — the smoothed-QR estimate $\hat\gamma_\tau$ of the per-allele effect on the pheno-transform scale (e.g. INT scale under `--pheno-transform int`).
- `SE_tau<val>` — the sandwich-variance standard error.
- `Z_tau<val>` and `P_tau<val>` — the Wald statistic $\hat\gamma_\tau / \widehat{\mathrm{SE}}$ and the two-sided normal $p$-value.


Wald mode defaults to a narrower bandwidth than score mode (`--spasqr-h-scale 10` vs `3`) to keep smoothing bias on $\hat\gamma_\tau$ low at the cost of slightly slower convergence; override via `--spasqr-h-scale <k>` when desired.
