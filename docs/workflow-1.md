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
Although SPA<sub>SQR</sub> is a quantile GWAS method, LOCO PGS computed using *linear* GWAS software work very well for our purpose. We therefore outsource the task of PGS construction to either [**LDAK-KVIK**](https://dougspeed.com/ldak-kvik/) or [**REGENIE**](https://rgcgithub.github.io/regenie/).


## Data that you will need

Throughout this page we use the same set of files:

```
geno.{bed,bim,fam}    PLINK fileset of variants to test
pheno.txt             phenotype file, two traits bmi and ldl
covar.txt             covariate file with columns age and sex
```

For concreteness, `pheno.txt` has the following first few rows

```
$ head pheno.txt
FID    IID         bmi      ldl
fam1   sample1    24.3    138.4
fam1   sample2    31.7    102.5
fam2   sample3    27.5    156.2
fam3   sample4    22.1    119.0
```

and `covar.txt` similarly:

```
$ head covar.txt
FID    IID         age    sex
fam1   sample1      52      1
fam1   sample2      45      0
fam2   sample3      61      1
fam3   sample4      38      0
```

The FID and IID columns of `pheno.txt` and `covar.txt` should match those in `geno.fam`. These simulated data are available in our Github Repository for users to replicate this tutorial. Note that GRAB automatically adds an intercept to the covariates, so there is no need to include an intercept column in `covar.txt`.

For simplicity, this tutorial uses the same `geno.{bed,bim,fam}` fileset for both LOCO PGS construction and SPA<sub>SQR</sub> association testing. In practice, the two stages typically use different variant subsets: LOCO PGS is usually trained on directly genotyped (unimputed) variants — a few hundred thousand high-quality SNPs — since imputed variants add measurement noise and slow down training without materially improving the accuracy of the LOCO PGS. Association testing, by contrast, is run on the full imputed dataset (tens of millions of variants, including rare ones) to maximize discovery. 

For GRAB to utilize the LOCO PGS, we also need a small text file called a **prediction list** — a two-column table pairing each phenotype name (column 1) with the absolute path to its LOCO PGS file (column 2). The prediction list is what GRAB reads (via the `--pred-list` argument) so that it knows which PGS file is paired with which trait. This file format mimics the workflow of REGENIE.

## Inverse normal transformation

Before computing the PGS we recommend applying a **rank-based inverse normal transformation (INT)** to each trait's non-missing values, because in real data analysis we find that doing so generally yields more associations than skipping it; this is likely because the LOCO PGS are more accurate when computed on an INT-transformed trait. GRAB offers a utility for applying INT on the phenotype file:

```bash
./grab --int-pheno --pheno pheno.txt --out pheno_int
```

The output `pheno_int.txt` contains the normalized phenotypes, with missing values unchanged:

```
$ head pheno_int.txt
FID    IID         bmi          ldl
fam1   sample1    -0.299         0.522
fam1   sample2     1.049        -1.064
fam2   sample3     0.299         1.281
fam3   sample4    -0.842        -0.674
```

## Computing the LOCO PGS with LDAK-KVIK

We can compute LOCO PGS via LDAK-KVIK step 1 with:

```bash
./ldak6.2.linux \
    --kvik-step1 ldak_step1 \
    --bfile geno \
    --pheno pheno_int.txt \
    --mpheno ALL \
    --covar covar.txt \
    --max-threads 8
```

With `--mpheno ALL`, LDAK-KVIK produces one LOCO PGS file per phenotype column in `pheno_int.txt`, named by the trait's **position** in the phenotype file rather than its column name:

```
ldak_step1.step1.pheno1.loco.prs        (LOCO PGS for bmi)
ldak_step1.step1.pheno2.loco.prs        (LOCO PGS for ldl)
```

Within each LOCO PGS file, every row is one subject's LOCO PGS across the 22 autosomes:

```
$ head -3 ldak_step1.step1.pheno1.loco.prs
FID    IID         Chr1     Chr2     Chr3   ...   Chr22
fam1   sample1     0.124   -0.083    0.211         0.196
fam1   sample2    -0.241    0.018   -0.135        -0.057
```

Unlike REGENIE, LDAK-KVIK does not emit a prediction list pairing each phenotype with its LOCO PGS file. One needs to create `ldak_pred_list.txt` manually with the following shell snippet:

```bash
cat > ldak_pred_list.txt <<EOF
bmi	$(pwd)/ldak_step1.step1.pheno1.loco.prs
ldl	$(pwd)/ldak_step1.step1.pheno2.loco.prs
EOF
```

`$(pwd)` expands to the current working directory so each entry ends up as an absolute path.

## Computing the LOCO PGS with REGENIE

We may also compute LOCO PGS via REGENIE Step 1:

```bash
./regenie \
    --step 1 \
    --bed geno \
    --phenoFile pheno_int.txt \
    --covarFile covar.txt \
    --bsize 1000 \
    --out regenie_step1
```

This produces one LOCO PGS file per trait plus REGENIE's own prediction list:

```
regenie_step1_1.loco          (LOCO PGS for bmi)
regenie_step1_2.loco          (LOCO PGS for ldl)
regenie_step1_pred.list       (pairs bmi and ldl with their .loco files)
```

REGENIE's `.loco` format is the **transpose** of LDAK-KVIK's: each row is one chromosome and each column is one subject (with FID and IID joined into a single `FID_IID` token):

```
$ head -3 regenie_step1_1.loco
FID_IID       fam1_sample1   fam1_sample2   fam2_sample3   ...
1             0.0497         -0.0624         -0.0156
2             0.0502         -0.0558         -0.0152
```
GRAB auto-detects the format from the file header and can distinguish whether the LOCO PGS files are produced by LDAK-KVIK or REGENIE.

Unlike with LDAK-KVIK, REGENIE produces `regenie_step1_pred.list` which can be passed directly to GRAB's `--pred-list` option:

```
$ cat regenie_step1_pred.list
bmi	/abs/path/to/regenie_step1_1.loco
ldl	/abs/path/to/regenie_step1_2.loco
```

## Running association testing with GRAB

Once the prediction list is ready, null model fitting (SPA<sub>SQR</sub> step 1) and association testing (SPA<sub>SQR</sub> 2) is completed in a single `grab` call:

```bash
./grab --method SPAsqr \
    --bfile geno \
    --pheno pheno_int.txt \
    --covar covar.txt \
    --pred-list ldak_pred_list.txt \
    --out spasqr_results
```

Below are the required and optional flags.

**Required:**

| Flag | What it does |
| --- | --- |
| `--method` | Selects the GRAB method to run; use `SPAsqr` to trigger SPA<sub>SQR</sub> and all of the `--spasqr-*` options below. |
| `--bfile` | PLINK 1 genotype fileset prefix (e.g. `geno` for `geno.{bed,bim,fam}`). PLINK 2 (`--pfile`), VCF (`--vcf`), and BGEN (`--bgen`) are also accepted — exactly one of the four is needed. |
| `--pheno` | Phenotype file (e.g. `pheno_int.txt`). Starts with `FID` and `IID`. |
| `--covar` | Covariate file (e.g. `covar.txt`). Starts with `FID` and `IID`. GRAB adds an intercept automatically — do not include one in the file. |
| `--out` | Output prefix (e.g. `spasqr_results`). Each trait gets its own tab-delimited result file. |

**Optional:**

| Flag | Default | What it does |
| --- | --- | --- |
| `--pred-list` | — | Prediction list (e.g. `ldak_pred_list.txt` for LDAK-KVIK or `regenie_step1_pred.list` for REGENIE). Omit to run with no LOCO offset — valid but much less powerful. |
| `--pheno-transform` | `int` | One of `int` / `standardize`. **Must match the transform used during PGS construction.** With `pheno_int.txt` and the INT workflow, leave it at the default. With raw `pheno.txt` fed to LDAK-KVIK or REGENIE, set this to `standardize`. |
| `--pheno-name` | all trait columns | Comma-separated list of trait columns to analyze (e.g. `bmi,ldl`). Omit to test every trait column found in `--pheno`. |
| `--covar-name` | all covar columns | Comma-separated list of covariate columns to use (e.g. `age,sex`). Omit to use all columns. |
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

Variants that fail the above QC constraints are omitted from GWAS results with $p$-values and $Z$-scores filled by NA.

GRAB writes one output file per phenotype:

```
spasqr_results.bmi.SPAsqr
spasqr_results.ldl.SPAsqr
```

Each file lists per-variant statistics including per-quantile $p$-values, the Cauchy-combined $P_\mathrm{CCT}$, and per-$\tau$ $Z$-scores:

```
$ head -3 spasqr_results.bmi.SPAsqr
CHROM  POS      ID          REF  ALT  MISS_RATE   ALT_FREQ  MAC    HWE_P     P_CCT      P_tau0.1  P_tau0.3   P_tau0.5   P_tau0.7    P_tau0.9    Z_tau0.1    Z_tau0.3   Z_tau0.5  Z_tau0.7  Z_tau0.9
1      1171417  rs6603782   C    T    0.0137558   0.326478  25748  0.536365  0.0071169  0.923124  0.397414   0.0379922  0.00415006  0.00223711  -0.0964998  -0.846248  -2.07494  -2.86651  -3.0567
1      2236359  rs60363208  G    A    0.00435185  0.173843  13841  0.412157  0.227344   0.273908  0.0829882  0.152614   0.509465    0.70161     1.09411     1.73361    1.43036   0.65967   -0.383148
```

The first nine columns are the variant identifier and standard variant-level QC fields. `P_CCT` is the Cauchy-combined $p$-value across all $\tau$ levels and is the main genome-wide significance result. The `P_tauX` and `Z_tauX` columns give the per-quantile $p$-values and $Z$-scores; comparing them across $\tau$ yields insight into the heterogeneity of effect sizes across the phenotype distribution.

### End-to-end recipes (with INT)

LDAK-KVIK LOCO PGS with INT:

```bash
# 1. INT-transform the phenotype
./grab --int-pheno --pheno pheno.txt --out pheno_int

# 2. Train the LOCO PGS on the INT-transformed Y
./ldak6.2.linux \
    --kvik-step1 ldak_step1 \
    --bfile geno \
    --pheno pheno_int.txt --mpheno ALL \
    --covar covar.txt \
    --max-threads 8

# 3. Build the pred-list
cat > ldak_pred_list.txt <<EOF
bmi	$(pwd)/ldak_step1.step1.pheno1.loco.prs
ldl	$(pwd)/ldak_step1.step1.pheno2.loco.prs
EOF

# 4. Run SPAsqr; --pheno-transform int is the default and matches the INT-trained PGS
./grab --method SPAsqr \
    --bfile geno \
    --pheno pheno_int.txt \
    --covar covar.txt \
    --pred-list ldak_pred_list.txt \
    --out spasqr_results
```

REGENIE LOCO PGS with INT:

```bash
# 1. INT-transform the phenotype
./grab --int-pheno --pheno pheno.txt --out pheno_int

# 2. Train the LOCO PGS on the INT-transformed Y
./regenie \
    --step 1 \
    --bed geno \
    --phenoFile pheno_int.txt \
    --covarFile covar.txt \
    --bsize 1000 \
    --out regenie_step1

# 3. Run SPAsqr with REGENIE's native pred-list
./grab --method SPAsqr \
    --bfile geno \
    --pheno pheno_int.txt \
    --covar covar.txt \
    --pred-list regenie_step1_pred.list \
    --out spasqr_results
```

## Skipping the INT

We may also skip the INT and feed the raw `pheno.txt` directly to LDAK-KVIK or REGENIE. Before fitting the LOCO PGS, both LDAK-KVIK and REGENIE internally regress the covariates out of the trait and then standardize the residuals to mean zero and unit variance, so the LOCO PGS still live on a **standardized scale**, instead of the scale of the raw `Y` column. Thus, we should pass `--pheno-transform standardize` to GRAB at association testing so that the trait and the LOCO PGS live on the same scale.

LDAK-KVIK LOCO PGS without INT:

```bash
# 1. Train the LOCO PGS on raw Y
./ldak6.2.linux \
    --kvik-step1 ldak_step1 \
    --bfile geno \
    --pheno pheno.txt --mpheno ALL \
    --covar covar.txt \
    --max-threads 8

# 2. Build the pred-list
cat > ldak_pred_list.txt <<EOF
bmi	$(pwd)/ldak_step1.step1.pheno1.loco.prs
ldl	$(pwd)/ldak_step1.step1.pheno2.loco.prs
EOF

# 3. Run SPAsqr with --pheno-transform standardize
./grab --method SPAsqr \
    --bfile geno \
    --pheno pheno.txt \
    --covar covar.txt \
    --pred-list ldak_pred_list.txt \
    --pheno-transform standardize \
    --out spasqr_results
```

REGENIE LOCO PGS without INT:

```bash
# 1. Train the LOCO PGS on raw Y
./regenie \
    --step 1 \
    --bed geno \
    --phenoFile pheno.txt \
    --covarFile covar.txt \
    --bsize 1000 \
    --out regenie_step1

# 2. Run SPAsqr with REGENIE's native pred-list and --pheno-transform standardize
./grab --method SPAsqr \
    --bfile geno \
    --pheno pheno.txt \
    --covar covar.txt \
    --pred-list regenie_step1_pred.list \
    --pheno-transform standardize \
    --out spasqr_results
```
