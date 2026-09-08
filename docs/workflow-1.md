---
layout: default
title: "Workflow 1: LOCO PGS + SPA<sub>SQR</sub>"
nav_order: 3
description: "End-to-end recipe for running SPAsqr with a LOCO polygenic score (without a sparse GRM)."
has_children: false
---

# **Workflow 1: LOCO PGS + SPA<sub>SQR</sub>**

SPA<sub>SQR</sub> uses **leave-one-chromosome-out (LOCO) polygenic scores (PGS)** as an offset: for each chromosome, the PGS built from all *other* chromosomes is subtracted from the trait before fitting the null model. This controls for relatedness and substantially improves power. LOCO PGS from linear-model software work well for this purpose, so we build them with [**LDAK-KVIK**](https://dougspeed.com/ldak-kvik/) or [**REGENIE**](https://rgcgithub.github.io/regenie/).

## Inputs

All files are in the [`data/`](https://github.com/qhengncsu/SPAsqr.github.io/tree/main/data) folder of this site.

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

We test `Quantitative1` and `Quantitative2` adjusting for `MALE` and `PC1`–`PC4`. GRAB adds the intercept automatically. The phenotype file must start with `FID IID` for LDAK-KVIK and REGENIE (GRAB alone also accepts a single `IID` column).

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
    --pheno-transform int \
    --threads 8 \
    --out spasqr_results
```

## Step by step

### 1. Rank-based inverse normal transformation

```bash
./grab2 --int-pheno --pheno simu_geno.pheno --pheno-name Quantitative1,Quantitative2 --out simu_geno_int
```

We first apply a rank-based inverse normal transformation (RINT) to each trait. In UK Biobank, RINT generally yields more associations than raw traits, possibly because the LOCO PGS explains more variance of the RINT-transformed trait. The command writes `simu_geno_int.txt`, which keeps the `FID IID` columns and replaces each trait column with its RINT version:

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

We then train the LOCO PGS on the RINT-transformed traits with LDAK-KVIK. It writes one LOCO PGS file per trait, named by the trait's **position** in the phenotype file: `ldak_step1.step1.pheno1.loco.prs` for Quantitative1 and `ldak_step1.step1.pheno2.loco.prs` for Quantitative2. Each file has one row per subject and one column per chromosome, holding that subject's PGS built from all other chromosomes:

```
$ head -3 ldak_step1.step1.pheno1.loco.prs
FID     IID     Chr1     Chr2     Chr3     ...  Chr22
S00001  S00001  0.2954   0.3367   0.2901        0.3047
S00002  S00002  -0.3886  -0.4240  -0.3571       -0.4109
```

For simplicity we use the same genotype file for LOCO PGS training and association testing. In practice the two stages usually differ: LOCO PGS training typically uses a few hundred thousand genotyped SNPs, while association testing uses the full imputed set.

### 3. Build the prediction list

```bash
cat > simu_geno_ldak_pred.list <<EOF
Quantitative1   $(pwd)/ldak_step1.step1.pheno1.loco.prs
Quantitative2   $(pwd)/ldak_step1.step1.pheno2.loco.prs
EOF
```

GRAB locates the LOCO PGS files through a **prediction list**, a two-column text file that pairs each trait name with the absolute path to its LOCO PGS file. We write it by hand here, using `$(pwd)` to expand the current directory into an absolute path. The format follows REGENIE's `pred.list`.

### 4. Run SPA<sub>SQR</sub>

```bash
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
    --pheno-transform int \
    --threads 8 \
    --out spasqr_results
```

Finally we run SPA<sub>SQR</sub>. Null-model fitting and association testing happen in a single call. `--spasqr-taus` sets the quantile levels to test, and `--pheno-transform int` makes GRAB apply RINT to the trait it reads. Here it is redundant, since `simu_geno_int.txt` is already RINT-transformed (see [Raw phenotypes](#raw-phenotypes-without-rint) if you skip RINT). GRAB writes one result file per trait: `spasqr_results.Quantitative1.SPAsqr` and `spasqr_results.Quantitative2.SPAsqr`.

## Output format

Each result file has one row per variant. With nine quantiles it has 37 columns: 10 fixed columns with variant information and the combined $p$-value, followed by three blocks of nine per-quantile columns (`P_tau`, `Z_tau`, `Z_Norm_tau`).

```
CHROM  POS     ID        REF  ALT  MISS_RATE  ALT_FREQ  MAC   HWE_P    P_CCT
7      62000   SNP_1428  A    G    0          0.4624    4624  0.1246   3.96e-08
10     122000  SNP_2170  A    G    0          0.267     2670  0.1933   1.60e-06
```

```
P_tau0.1   P_tau0.2   P_tau0.3   P_tau0.4   P_tau0.5   P_tau0.6   P_tau0.7   P_tau0.8   P_tau0.9
1.82e-06   8.23e-09   8.31e-09   3.32e-08   1.21e-07   3.21e-07   6.56e-07   9.79e-07   9.72e-07
3.62e-08   8.94e-07   1.83e-06   3.21e-06   2.77e-06   2.20e-06   2.06e-06   2.05e-06   1.27e-06
```

```
Z_tau0.1  Z_tau0.2  Z_tau0.3  Z_tau0.4  Z_tau0.5  Z_tau0.6  Z_tau0.7  Z_tau0.8  Z_tau0.9
+4.77     +5.76     +5.76     +5.52     +5.29     +5.11     +4.97     +4.89     +4.89
-5.51     -4.92     -4.78     -4.66     -4.70     -4.74     -4.75     -4.75     -4.84
```

| Columns | Meaning |
| --- | --- |
| `P_CCT` | Cauchy-combined $p$-value across all quantiles. |
| `P_tau<τ>` | Saddlepoint $p$-value at each quantile. |
| `Z_tau<τ>` | Signed $Z$-score consistent with `P_tau`: `sign(Z_Norm_tau) × Φ⁻¹(1 − P_tau/2)`. |
| `Z_Norm_tau<τ>` | Raw score statistic $S/\sqrt{\operatorname{Var}(S)}$, before saddlepoint correction. Nearly identical to `Z_tau` for common variants; differs for rare ones. |

## Other options

### REGENIE instead of LDAK-KVIK

REGENIE can replace steps 2 and 3 with a single call, because it writes the prediction list itself.

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
    --pheno-transform int \
    --threads 8 \
    --out spasqr_results
```

REGENIE's `.loco` files are the transpose of LDAK-KVIK's: one row per chromosome and one column per subject, with `FID` and `IID` joined by an underscore. GRAB detects which format it is reading from the header.

### Raw phenotypes (without RINT)

To avoid RINT, skip step 1 and feed the raw `simu_geno.pheno` to LDAK-KVIK or REGENIE directly. Both tools regress out the covariates and standardize the residuals internally, so the LOCO PGS live on a standardized scale. Pass `--pheno-transform standardize` so that GRAB puts the trait on the same scale:

```bash
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno.pheno --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --pheno-transform standardize \
    --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
    --threads 8 \
    --out spasqr_results
```

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
| `--pred-list` | — | LOCO PGS prediction list. Omit to run without an offset (valid but much less powerful). |
| `--pheno-transform` | `int` | Transformation GRAB applies to the trait: `int` or `standardize`. **Must match the transformation applied to the trait during LOCO PGS training.** Redundant when the input trait is already RINT-transformed, as in this workflow. |
| `--pheno-name` | all trait columns | Traits to test, comma-separated. |
| `--covar` | — | Covariate file; may be the same file as `--pheno`. |
| `--covar-name` | — | Covariate columns, comma-separated. |
| `--spasqr-taus` | `0.1,0.3,0.5,0.7,0.9` | Quantile levels to test (max 20). |
| `--spasqr-h-scale` | `3` | Bandwidth divisor, `h = IQR / k`. Larger = less smoothing. |
| `--threads` | `1` | Number of threads. |

SNP filters:

| Flag | Default | What it does |
| --- | --- | --- |
| `--maf` | `1e-5` | Minimum minor allele frequency. |
| `--mac` | `10` | Minimum minor allele count. |
| `--geno` | `0.1` | Maximum per-variant missingness. |
| `--extract` | — | File of variant IDs to test, one per line. |
| `--chr` | all autosomes | Chromosomes to test, comma-separated. |
