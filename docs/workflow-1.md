---
layout: default
title: "Workflow 1: LOCO PGS + SPA<sub>SQR</sub>"
nav_order: 3
description: "End-to-end recipe for running SPAsqr with a LOCO PGS (without a sparse GRM)."
has_children: false
---

# **Workflow 1: LOCO PGS + SPA<sub>SQR</sub>**

SPA<sub>SQR</sub> uses **leave-one-chromosome-out (LOCO) polygenic scores (PGS)** to account for relatedness and substantially improve power. Before fitting the null model for a chromosome, it subtracts the LOCO PGS trained on all other chromosomes from the trait. LOCO PGS from linear models work well as offsets. Here we use [**LDAK-KVIK**](https://dougspeed.com/ldak-kvik/) or [**REGENIE**](https://rgcgithub.github.io/regenie/) to train them.

## Inputs

The [`data/`](https://github.com/qhengncsu/SPAsqr.github.io/tree/main/data) folder contains all the example inputs.

```
simu_geno.{bed,bim,fam}    PLINK 1 genotypes (5000 subjects, 5000 variants, 22 autosomes, 1250 families of 4)
simu_geno.pheno            phenotypes + covariates
```

```
$ head -3 simu_geno.pheno
FID     IID     MALE  PC1         PC2         PC3          PC4         Quantitative1  Quantitative2
S00001  S00001  0     0.0065558   -0.0190989  0.00331922   0.00574267  0.447611215    -1.38874574
S00002  S00002  1     0.00947819  -0.0120386  -0.0226929   0.0132888   -1.28469274    0.626376502
```

We test two traits, `Quantitative1` and `Quantitative2`, adjusting for `MALE` and `PC1`–`PC4`. GRAB adds the intercept automatically. LDAK-KVIK and REGENIE expect the phenotype file to start with `FID IID`. If you use GRAB alone, a single `IID` column is also sufficient.

## Complete pipeline

```bash
# 1. Rank-based inverse-normal-transform the traits
./grab2 --int-pheno --pheno simu_geno.pheno --pheno-name Quantitative1,Quantitative2 --out simu_geno_int

# 2. Train the LOCO PGS with LDAK-KVIK
./ldak6.2.linux \
    --kvik-step1 ldak_step1 \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --mpheno ALL \
    --covar simu_geno.pheno   --covar-names MALE,PC1,PC2,PC3,PC4 \
    --max-threads 8

# 3. Build the prediction list
cat > simu_geno_ldak_pred.list <<EOF
Quantitative1   $(pwd)/ldak_step1.step1.pheno1.loco.prs
Quantitative2   $(pwd)/ldak_step1.step1.pheno2.loco.prs
EOF

# 4. Run SPAsqr
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
    --threads 8 \
    --out spasqr_results
```

## Step by step

### 1. Rank-based inverse normal transformation

```bash
./grab2 --int-pheno --pheno simu_geno.pheno --pheno-name Quantitative1,Quantitative2 --out simu_geno_int
```

We first apply a rank-based inverse normal transformation (RINT) to each trait. In our UK Biobank analyses, we generally found more associations after RINT, possibly because the LOCO PGS explained more of the variance in the transformed trait. The command writes `simu_geno_int.txt` with the original `FID IID` columns and the transformed trait values.

```
$ head -3 simu_geno_int.txt
FID     IID     Quantitative1   Quantitative2
S00001  S00001  0.40780638     -1.39757115
S00002  S00002  -1.28887908    0.636848104
```

### 2. Train the LOCO PGS

```bash
./ldak6.2.linux \
    --kvik-step1 ldak_step1 \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --mpheno ALL \
    --covar simu_geno.pheno   --covar-names MALE,PC1,PC2,PC3,PC4 \
    --max-threads 8
```

LDAK-KVIK trains the LOCO PGS on the transformed traits and writes one file per trait. The file name follows the trait order in the phenotype file. Here, `ldak_step1.step1.pheno1.loco.prs` is for Quantitative1 and `ldak_step1.step1.pheno2.loco.prs` is for Quantitative2. Each row represents a participant. Each chromosome column gives that participant's PGS trained on all other chromosomes.

```
$ head -3 ldak_step1.step1.pheno1.loco.prs
FID     IID     Chr1     Chr2     Chr3     ...  Chr22
S00001  S00001  0.2954   0.3367   0.2901        0.3047
S00002  S00002  -0.3886  -0.4240  -0.3571       -0.4109
```

We use the same genotype file for both steps in this example. In practice, the LOCO PGS is usually trained on a few hundred thousand genotyped SNPs, and association testing uses the full set of imputed variants.

### 3. Build the prediction list

```bash
cat > simu_geno_ldak_pred.list <<EOF
Quantitative1   $(pwd)/ldak_step1.step1.pheno1.loco.prs
Quantitative2   $(pwd)/ldak_step1.step1.pheno2.loco.prs
EOF
```

The **prediction list** tells GRAB where to find each trait's LOCO PGS file. It uses the same two-column format as REGENIE's `pred.list`: the trait name followed by the absolute path to the file. Here we create the list by hand and use `$(pwd)` to insert the full path to the current directory.

### 4. Run SPA<sub>SQR</sub>

```bash
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
    --threads 8 \
    --out spasqr_results
```

Use the same phenotype file you used to train the LOCO PGS (see [LOCO PGS scale](#loco-pgs-scale)), and specify the quantiles to test with `--spasqr-taus`. A single GRAB command fits the null models, tests the variants, and writes a result file for each trait. In this example, the files are `spasqr_results.Quantitative1.SPAsqr` and `spasqr_results.Quantitative2.SPAsqr`.

## Output format

Each row gives the results for one variant. With nine quantiles, the default output has 28 columns. The first 10 contain variant information and the combined $p$-value. Two sets of nine columns follow, one for the quantile-specific $p$-values (`LOG10P_tau`) and one for the score statistics (`Z_Norm_tau`). Use `--output-z` to add the nine `Z_tau` columns, bringing the total to 37.

GRAB reports all $p$-values as $-\log_{10}p$, so a larger value indicates stronger evidence of association. The genome-wide threshold $p < 5\times10^{-8}$ corresponds to `LOG10P > 7.30`. To recover the $p$-value, use `10^(-LOG10P)`. The example below shows the rows for `SNP_1428` and `SNP_2170` in `spasqr_results.Quantitative1.SPAsqr`.

```
CHROM  POS     ID        REF  ALT  MISS_RATE  ALT_FREQ  MAC   LOG10P_HWE  LOG10P_CCT
7      62000   SNP_1428  A    G    0          0.4624    4624  0.904       7.178
10     122000  SNP_2170  A    G    0          0.267     2670  0.714       5.639
```

```
LOG10P_tau0.1  LOG10P_tau0.2  LOG10P_tau0.3  LOG10P_tau0.4  LOG10P_tau0.5  LOG10P_tau0.6  LOG10P_tau0.7  LOG10P_tau0.8  LOG10P_tau0.9
5.682          7.823          7.706          7.060          6.497          6.105          5.847          5.741          5.827
6.374          6.140          5.180          4.184          3.491          3.200          3.258          3.480          3.434
```

```
Z_Norm_tau0.1  Z_Norm_tau0.2  Z_Norm_tau0.3  Z_Norm_tau0.4  Z_Norm_tau0.5  Z_Norm_tau0.6  Z_Norm_tau0.7  Z_Norm_tau0.8  Z_Norm_tau0.9
+4.74          +5.66          +5.61          +5.35          +5.11          +4.94          +4.82          +4.77          +4.81
-5.08          -4.96          -4.51          -3.99          -3.60          -3.42          -3.45          -3.59          -3.56
```

| Columns | Meaning |
| --- | --- |
| `LOG10P_HWE` | $-\log_{10}$ of the Hardy–Weinberg equilibrium $p$-value. |
| `LOG10P_CCT` | $-\log_{10}$ of the Cauchy-combined $p$-value across all quantiles. Computed in the log domain, so it stays finite even when $p$ would underflow below $10^{-308}$. |
| `LOG10P_tau<τ>` | $-\log_{10}$ of the saddlepoint $p$-value at each quantile. |
| `Z_Norm_tau<τ>` | Raw score statistic $S/\sqrt{\operatorname{Var}(S)}$, before saddlepoint correction. |
| `Z_tau<τ>` | Written with `--output-z`. Signed $Z$-score consistent with `LOG10P_tau`: `sign(Z_Norm_tau) × Φ⁻¹(1 − p/2)`. Nearly identical to `Z_Norm_tau` for common variants; differs for rare ones. |

## Other options

### REGENIE instead of LDAK-KVIK

REGENIE writes both the LOCO PGS and the prediction list, so one command replaces steps 2 and 3.

```bash
# 2. Train the LOCO PGS with REGENIE (writes simu_geno_regenie_pred.list)
./regenie \
    --step 1 \
    --bed simu_geno \
    --phenoFile simu_geno_int.txt --phenoColList Quantitative1,Quantitative2 \
    --covarFile simu_geno.pheno   --covarColList MALE,PC1,PC2,PC3,PC4 \
    --bsize 1000 --threads 8 \
    --out simu_geno_regenie

# 3. Run SPAsqr with REGENIE's prediction list
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_regenie_pred.list \
    --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
    --threads 8 \
    --out spasqr_results
```

In REGENIE's `.loco` files, each row is a chromosome and each column is a participant, the reverse of the layout in LDAK-KVIK. REGENIE joins `FID` and `IID` with an underscore. GRAB recognizes the format from the header and automatically handles the difference in [LOCO PGS scale](#loco-pgs-scale).

### LOCO PGS scale

The phenotype and the LOCO PGS must be in the same scale. Give GRAB the same phenotype file and covariates you used for training the LOCO PGS, especially if you applied RINT before training.

- **LDAK-KVIK** reports the LOCO PGS on the training phenotype scale. GRAB uses the phenotype as supplied.
- **REGENIE** reports the LOCO PGS in units of the covariate-adjusted residual standard deviation. GRAB automatically rescales the phenotype to match.

Without `--pred-list`, GRAB applies RINT by default. To use the phenotype as supplied, set `--pheno-transform raw`.

### Raw phenotypes (without RINT)

To analyze the traits on their original scale, skip step 1. Use `simu_geno.pheno` to train the LOCO PGS with LDAK-KVIK or REGENIE, then pass the same file to GRAB.

```bash
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno.pheno --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
    --threads 8 \
    --out spasqr_results
```

Make sure the prediction list points to the LOCO PGS trained on the raw traits.

The table below summarizes how to choose the phenotype file and transformation option for analyses with or without RINT.

| Apply RINT? | LOCO PGS supplied via `--pred-list` | File to pass to `--pheno` | `--pheno-transform` |
|---|---|---|---|
| Yes | LDAK-KVIK, trained on `simu_geno_int.txt` | `simu_geno_int.txt` | Omit (GRAB defaults to `raw`) |
| Yes | REGENIE, trained on `simu_geno_int.txt` | `simu_geno_int.txt` | Omit (GRAB defaults to `standardize`) |
| Yes | None | `simu_geno.pheno` | Omit (GRAB defaults to `int`) |
| No | LDAK-KVIK, trained on `simu_geno.pheno` | `simu_geno.pheno` | Omit (GRAB defaults to `raw`) |
| No | REGENIE, trained on `simu_geno.pheno` | `simu_geno.pheno` | Omit (GRAB defaults to `standardize`) |
| No | None | `simu_geno.pheno` | Specify `--pheno-transform raw` |

### `grab2 --method SPAsqr` flags

Required:

| Flag | What it does |
| --- | --- |
| `--method SPAsqr` | Selects SPA<sub>SQR</sub>. |
| `--bfile` | PLINK 1 prefix. `--pfile` (PLINK 2), `--vcf`, and `--bgen` are also accepted; exactly one is needed. |
| `--pheno` | Phenotype file starting with `FID IID` or `IID` (`#FID`/`#IID` headers also work). |
| `--out` | Output prefix. |

Optional:

| Flag | Default | What it does |
| --- | --- | --- |
| `--pred-list` | — | LOCO PGS prediction list. Omit to run without an offset. |
| `--pheno-transform` | depends on the LOCO PGS | Transformation GRAB applies to the trait: `int`, `standardize`, or `raw` (none). Defaults to `raw` with an LDAK-KVIK LOCO PGS and to `int` without `--pred-list`; always `standardize` with a REGENIE LOCO PGS. See [LOCO PGS scale](#loco-pgs-scale). |
| `--pheno-name` | all trait columns | Traits to test, comma-separated. |
| `--covar` | — | Covariate file; may be the same file as `--pheno`. |
| `--covar-name` | — | Covariate columns, comma-separated. |
| `--spasqr-taus` | `0.1,0.3,0.5,0.7,0.9` | Quantile levels to test (max 20). |
| `--spasqr-h-scale` | `3` | Bandwidth divisor in `h = IQR / k`. A larger divisor gives less smoothing. |
| `--output-z` | off | Also write the `Z_tau` columns. |
| `--threads` | `1` | Number of threads. |

SNP filters:

| Flag | Default | What it does |
| --- | --- | --- |
| `--maf` | `1e-5` | Minimum minor allele frequency. |
| `--mac` | `10` | Minimum minor allele count. |
| `--geno` | `0.1` | Maximum per-variant missingness. |
| `--extract` | — | File of variant IDs to test, one per line. |
| `--chr` | all autosomes | Chromosomes to test, comma-separated. |
