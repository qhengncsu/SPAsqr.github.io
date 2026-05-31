---
layout: default
title: "Workflow 1: LOCO PGS + SPA<sub>SQR</sub>"
nav_order: 3
description: "End-to-end recipe for running SPAsqr with a LOCO polygenic score (without a sparse GRM)."
has_children: false
---

# **Workflow 1: LOCO PGS + SPA<sub>SQR</sub>**

In this section we give a detailed tutorial on how to run SPA<sub>SQR</sub> association testing using **leave-one-chromosome-out (LOCO) polygenic scores (PGS)** as an offset. LOCO PGS are per-subject predictions of the trait built from variants on every chromosome *except* the one currently being tested. SPA<sub>SQR</sub> subtracts the chromosome-specific LOCO PGS from the trait before fitting the null smoothed quantile regression model: for each chromosome $c$ and quantile level $\tau$,

$$
\widetilde{\beta}_{h,\tau} \;=\; \arg\min_{\beta} \; \frac{1}{n} \sum_{i=1}^n \ell_{h,\tau}\!\left( \tilde Y_i - \hat Y_{-c,i} - X_i^{\!\top} \beta \right),
$$

where $\tilde Y$ is the trait after rank-based inverse normal transformation, $\hat Y_{-c}$ is the LOCO PGS for chromosome $c$, $X$ is the covariate matrix, and $\ell_{h,\tau}$ is the Gaussian-kernel convolution smoothed quantile check loss with bandwidth $h = \mathrm{IQR}(\tilde Y - \hat Y_{-c})/k$. Absorbing the LOCO PGS as an offset removes the polygenic background driven by other chromosomes, which controls for relatedness and substantially improves the statistical power of the per-variant score tests.

Although SPA<sub>SQR</sub> is a quantile GWAS method, LOCO PGS computed using linear GWAS software work very well for our purpose. We therefore outsource the task of PGS construction to either [**LDAK-KVIK**](https://dougspeed.com/ldak-kvik/) or [**REGENIE**](https://rgcgithub.github.io/regenie/).


## Data that you will need

Throughout this page we assume a working directory containing the genotype file and the phenotype/covariate file:

```
simu_geno.{bed,bim,fam}    PLINK 1 genotype fileset
simu_geno.pheno            phenotype + covariate columns
grab2                      GRAB binary
ldak6.2.linux              LDAK binary
regenie                    REGENIE binary
plink2                     PLINK 2 binary (used in Workflow 2 to build the sparse GRM)
```

The phenotype file `simu_geno.pheno` is in PLINK format, with `FID`/`IID` in the first two columns and phenotype/covariate data in the remaining columns:

```
$ head simu_geno.pheno
FID     IID     MALE  PC1         PC2         PC3          PC4         Quantitative1  Quantitative2
S00001  S00001  0     0.0065558   -0.0190989  0.00331922   0.00574267  0.447611215    -1.38874574
S00002  S00002  1     0.00947819  -0.0120386  -0.0226929   0.0132888   -1.28469274    0.626376502
```

A single file may carry both phenotype and covariate columns (as here), or the two may live in separate files supplied through `--pheno` and `--covar`. GRAB also accepts a single `IID` key column in place of the `FID IID` pair; LDAK-KVIK and REGENIE require the `FID IID` pair. The bundled `simu_geno.*` PLINK1 file set in this tutorial is a simulated 5000-subject × 5000-variant × 22-autosome data set (1250 families × 4 individual per family). `Quantitative1` and `Quantitative2` in the phenotype file carry genuine polygenic signals (heritability ≈ 0.30, 500 causal SNPs per trait). All input files are available in the [`data/`](https://github.com/qhengncsu/SPAsqr.github.io/tree/main/data) folder of this documentation repository; download them to replicate this tutorial.

In what follows, we wish to run QR GWAS for the two quantitative traits `Quantitative1` and `Quantitative2`, adjusting for the covariates `MALE`, `PC1`, `PC2`, `PC3`, `PC4`. GRAB automatically adds an intercept to the covariate matrix, so there is no need to include one in the phenotype file. For simplicity, the same `simu_geno.{bed,bim,fam}` fileset is used for both LOCO PGS construction and SPA<sub>SQR</sub> association testing. In practice, the two stages typically use different variant subsets: LOCO PGS is usually trained on directly genotyped (unimputed) variants — a few hundred thousand high-quality SNPs — since imputed variants add measurement noise and slow down training without materially improving the accuracy of the LOCO PGS. Association testing, by contrast, is run on the full imputed dataset (millions of variants, including rare ones) to maximize discovery.

For GRAB to utilize the LOCO PGS, we also need a small text file called a **prediction list** — a two-column table pairing each phenotype name (column 1) with the absolute path to its LOCO PGS file (column 2). The prediction list is what GRAB reads (via the `--pred-list` argument) so that it knows which PGS file is paired with which trait. This file format mimics the workflow of REGENIE.


## Inverse normal transformation

Before computing the PGS we recommend applying a **rank-based inverse normal transformation (INT)** to each trait's non-missing values, because in UK Biobank real data analysis we find that doing so generally yields more associations than skipping it; this is likely because LOCO PGS are more accurate when computed on an INT-transformed trait. GRAB offers a utility for applying INT to selected columns of the phenotype file:

```bash
./grab2 --int-pheno --pheno simu_geno.pheno --pheno-name Quantitative1,Quantitative2 --out simu_geno_int
```

The output `simu_geno_int.txt` retains the `FID IID` key columns and replaces each requested trait column with its INT-transformed version.

```
$ head simu_geno_int.txt
FID     IID     Quantitative1   Quantitative2
S00001  S00001  0.40780638     -1.39757115
S00002  S00002  -1.28887908    0.636848104
```

## Computing the LOCO PGS with LDAK-KVIK

We may compute LDAK LOCO PGS for both `Quantitative1` and `Quantitative2` using `MALE`, `PC1`, `PC2`, `PC3`, `PC4` as covariates via:

```bash
./ldak6.2.linux \
    --kvik-step1 ldak_step1 \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --mpheno ALL \
    --covar simu_geno.pheno   --covar-names MALE,PC1,PC2,PC3,PC4 \
    --max-threads 8
```

This writes one LOCO PGS file per phenotype, named by the trait's **position** in the phenotype file rather than its column name:

```
ldak_step1.step1.pheno1.loco.prs        (LOCO PGS for Quantitative1)
ldak_step1.step1.pheno2.loco.prs        (LOCO PGS for Quantitative2)
```

Each row of a LOCO PGS file is one subject's LOCO PGS, with one column per chromosome:

```
$ head -3 ldak_step1.step1.pheno1.loco.prs
FID     IID     Chr1     Chr2     Chr3     Chr4     ...  Chr22
S00001  S00001  0.2954   0.3367   0.2901   0.2952        0.3047
S00002  S00002  -0.3886  -0.4240  -0.3571  -0.4036       -0.4109
```


As mentioned,for GRAB to be able to properly utilize the LOCO PGS computed by LDAK-KVIK, we must manually create a prediction list:

```bash
cat > simu_geno_ldak_pred.list <<EOF
Quantitative1   $(pwd)/ldak_step1.step1.pheno1.loco.prs
Quantitative2   $(pwd)/ldak_step1.step1.pheno2.loco.prs
EOF
```

`$(pwd)` expands to the current working directory so here each entry ends up as the absolute path for the correspondling LOCO PGS. 


## Computing the LOCO PGS with REGENIE

Alternatively, we may compute LOCO PGS via REGENIE:

```bash
./regenie \
    --step 1 \
    --bed simu_geno \
    --phenoFile simu_geno_int.txt --phenoColList Quantitative1,Quantitative2 \
    --covarFile simu_geno.pheno   --covarColList MALE,PC1,PC2,PC3,PC4 \
    --bsize 1000 --threads 8 \
    --out simu_geno_regenie
```

This generates two LOCO PGS files plus REGENIE's own prediction list:

```
simu_geno_regenie_1.loco          (LOCO PGS for Quantitative1)
simu_geno_regenie_2.loco          (LOCO PGS for Quantitative2)
simu_geno_regenie_pred.list       (pairs Quantitative1 / Quantitative2 with their .loco files)
```

REGENIE's `.loco` format is the **transpose** of LDAK-KVIK's: each row is one chromosome and each column is one subject (with FID and IID joined into a single `FID_IID` token). REGENIE writes 22 (or 23 if `ChrX` is present in the genotype) chromosome rows; each row holds the leave-one-chromosome-out PGS for the corresponding chromosome:

```
$ head -5 simu_geno_regenie_1.loco
FID_IID  S00001_S00001  S00002_S00002  S00003_S00003  S00004_S00004  ...
1        -0.0089        -0.1859        -0.0925        -0.1799        ...
2        -0.0091        -0.1506        -0.1031        -0.1426        ...
3        -0.0602        -0.1983        -0.1319        -0.2331        ...
4        -0.0156        -0.2151        -0.0762        -0.2677        ...
```

GRAB auto-detects the format from each LOCO file's header and distinguishes whether the LOCO OGS file is produced by LDAK-KVIK or REGENIE.
Unlike LDAK-KVIK, REGENIE produces `simu_geno_regenie_pred.list` automatically, in the exact format `grab2 --pred-list` expects:

```
$ cat simu_geno_regenie_pred.list
Quantitative1   /abs/path/to/simu_geno_regenie_1.loco
Quantitative2   /abs/path/to/simu_geno_regenie_2.loco
```

## Running association testing with GRAB

Once the prediction list is ready, null model fitting (SPA<sub>SQR</sub> step 1) and association testing (SPA<sub>SQR</sub> step 2) are completed in a single `grab2` call:

```bash
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
    --threads 8 \
    --pheno-transform int \
    --out spasqr_results
```

Pass either the manually-created `simu_geno_ldak_pred.list` (LDAK-KVIK) or the REGENIE-auto-emitted `simu_geno_regenie_pred.list` (REGENIE) to `--pred-list`. Below are the required and optional flags.

**Required:**

| Flag | What it does |
| --- | --- |
| `--method` | Selects the GRAB method to run; use `SPAsqr` to trigger SPA<sub>SQR</sub> and all of the `--spasqr-*` options below. |
| `--bfile` | PLINK 1 genotype fileset prefix (e.g. `simu_geno` for `simu_geno.{bed,bim,fam}`). PLINK 2 (`--pfile`), VCF (`--vcf`), and BGEN (`--bgen`) are also accepted — exactly one of the four is needed. |
| `--pheno` | Phenotype file (e.g. `simu_geno_int.txt`). Starts with `FID IID` (or `IID`). `#FID #IID` and `#IID` also works. |
| `--out` | Output prefix (e.g. `spasqr_results`). Each trait gets its own tab-delimited result file. |

**Optional:**

| Flag | Default | What it does |
| --- | --- | --- |
| `--pred-list` | — | Prediction list (e.g. `simu_geno_ldak_pred.list` for LDAK-KVIK or `simu_geno_regenie_pred.list` for REGENIE). Omit to run with no LOCO offset — valid but much less powerful. |
| `--pheno-transform` | `int` | One of `int` / `standardize`. **Must match the transform used during PGS construction.** If `simu_geno_int.txt` is used to compute LOGO PGS, leave it at the default. With raw `simu_geno.pheno` fed to LDAK-KVIK or REGENIE, set this to `standardize`. |
| `--pheno-name` | all trait columns | Comma-separated list of trait columns to analyze (e.g. `Quantitative1,Quantitative2`). |
| `--covar` + `--covar-name` | — | Covariate file + comma-separated list of covariate columns (e.g. `MALE,PC1,PC2,PC3,PC4`). May point to the same file as `--pheno`. |
| `--spasqr-taus` | `0.1,0.3,0.5,0.7,0.9` | Quantile levels at which to test, comma-separated (max 20 levels). |
| `--spasqr-h-scale` | `3` (score mode) | Bandwidth divisor: $h = \mathrm{IQR}(\tilde Y - \hat Y_{-c}) / \text{scale}$. Larger value → less smoothing. |
| `--threads` | `1` | Number of threads used for parallel computing (e.g. `8`). |

**SNP filters:**

| Flag | Default | What it does |
| --- | --- | --- |
| `--maf` | `1e-5` | Minimum minor allele frequency (e.g. `0.01`). |
| `--mac` | `10` | Minimum minor allele count. |
| `--geno` | `0.1` | Maximum per-variant missingness fraction. |
| `--extract` | — | Restrict testing to the variant IDs listed in a file, one ID per line (e.g. `--extract snps.txt`). |
| `--chr` | all autosomes | Comma-separated chromosomes to test (e.g. `1,2,5`). |

Variants that fail the above QC constraints are omitted from GWAS results.

GRAB writes one output file per phenotype:

```
spasqr_results.Quantitative1.SPAsqr
spasqr_results.Quantitative2.SPAsqr
```

Each file has one row per variant and 28 columns. Below are the top two hits for `Quantitative1`, split into three blocks so the rows fit on the page — in the file these all sit on a single tab-separated line.

Columns 1–10 — variant info, QC fields, and the combined $p$-value `P_CCT`:

```
CHROM  POS     ID        REF  ALT  MISS_RATE  ALT_FREQ  MAC   HWE_P    P_CCT
7      62000   SNP_1428  A    G    0          0.4624    4624  0.1246   3.96e-08
10     122000  SNP_2170  A    G    0          0.267     2670  0.1933   1.60e-06
```

Columns 11–19 — per-quantile $p$-values, $\tau = 0.1, 0.2, \ldots, 0.9$:

```
P_tau0.1   P_tau0.2   P_tau0.3   P_tau0.4   P_tau0.5   P_tau0.6   P_tau0.7   P_tau0.8   P_tau0.9
1.82e-06   8.23e-09   8.31e-09   3.32e-08   1.21e-07   3.21e-07   6.56e-07   9.79e-07   9.72e-07
3.62e-08   8.94e-07   1.83e-06   3.21e-06   2.77e-06   2.20e-06   2.06e-06   2.05e-06   1.27e-06
```

Columns 20–28 — per-quantile $Z$-scores (signed), same $\tau$ order:

```
Z_tau0.1  Z_tau0.2  Z_tau0.3  Z_tau0.4  Z_tau0.5  Z_tau0.6  Z_tau0.7  Z_tau0.8  Z_tau0.9
+4.77     +5.76     +5.76     +5.52     +5.29     +5.11     +4.97     +4.89     +4.89
-5.51     -4.92     -4.78     -4.66     -4.70     -4.74     -4.75     -4.75     -4.84
```

`P_CCT` is the main result: the nine per-$\tau$ $p$-values combined into one via the Cauchy combination test. The `P_tauX` and `Z_tauX` columns let you see which quantiles the signal is coming from — e.g. SNP_1428 has positive $Z$ at every $\tau$, SNP_2170 has negative $Z$ at every $\tau$.


### End-to-end recipes (with INT)

LDAK-KVIK LOCO PGS with INT:

```bash
# 1. INT-transform the selected traits
./grab2 --int-pheno --pheno simu_geno.pheno --pheno-name Quantitative1,Quantitative2 --out simu_geno_int

# 2. Train the LOCO PGS on the INT-transformed Y
./ldak6.2.linux \
    --kvik-step1 ldak_step1 \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --mpheno ALL \
    --covar simu_geno.pheno   --covar-names MALE,PC1,PC2,PC3,PC4 \
    --max-threads 8

# 3. Build the LDAK pred-list
cat > simu_geno_ldak_pred.list <<EOF
Quantitative1   $(pwd)/ldak_step1.step1.pheno1.loco.prs
Quantitative2   $(pwd)/ldak_step1.step1.pheno2.loco.prs
EOF

# 4. Run SPAsqr; --pheno-transform int is the default and matches the INT-trained PGS
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --out spasqr_results
```

REGENIE LOCO PGS with INT:

```bash
# 1. INT-transform the selected traits
./grab2 --int-pheno --pheno simu_geno.pheno --pheno-name Quantitative1,Quantitative2 --out simu_geno_int

# 2. Train the LOCO PGS on the INT-transformed Y
./regenie \
    --step 1 \
    --bed simu_geno \
    --phenoFile simu_geno_int.txt --phenoColList Quantitative1,Quantitative2 \
    --covarFile simu_geno.pheno   --covarColList MALE,PC1,PC2,PC3,PC4 \
    --bsize 1000 --threads 8 \
    --out simu_geno_regenie

# 3. Run SPAsqr with REGENIE's native pred-list
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_regenie_pred.list \
    --out spasqr_results
```


## Skipping the INT

We may also skip the INT and feed the raw `simu_geno.pheno` directly to LDAK-KVIK or REGENIE. Before fitting the LOCO PGS, both LDAK-KVIK and REGENIE internally regress the covariates out of the trait and then standardize the residuals to mean zero and unit variance, so the LOCO PGS still live on a **standardized scale**, instead of the scale of the raw `Y` column. Thus, we should pass `--pheno-transform standardize` to GRAB at association testing so that the trait and the LOCO PGS live on the same scale.

LDAK-KVIK LOCO PGS without INT:

```bash
# 1. Train the LOCO PGS on raw Y (--mpheno selects specific phenotype columns by position)
./ldak6.2.linux \
    --kvik-step1 ldak_step1 \
    --bfile simu_geno \
    --pheno simu_geno.pheno --mpheno ALL \
    --covar simu_geno.pheno --covar-names MALE,PC1,PC2,PC3,PC4 \
    --max-threads 8

# 2. Build the LDAK pred-list
cat > simu_geno_ldak_pred.list <<EOF
Quantitative1   $(pwd)/ldak_step1.step1.pheno1.loco.prs
Quantitative2   $(pwd)/ldak_step1.step1.pheno2.loco.prs
EOF

# 3. Run SPAsqr with --pheno-transform standardize
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno.pheno --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --pheno-transform standardize \
    --out spasqr_results
```

REGENIE LOCO PGS without INT:

```bash
# 1. Train the LOCO PGS on raw Y
./regenie \
    --step 1 \
    --bed simu_geno \
    --phenoFile simu_geno.pheno --phenoColList Quantitative1,Quantitative2 \
    --covarFile simu_geno.pheno --covarColList MALE,PC1,PC2,PC3,PC4 \
    --bsize 1000 --threads 8 \
    --out simu_geno_regenie

# 2. Run SPAsqr with REGENIE's native pred-list and --pheno-transform standardize
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno.pheno --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_regenie_pred.list \
    --pheno-transform standardize \
    --out spasqr_results
```
