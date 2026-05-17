---
layout: default
title: "Workflow 1: LOCO PGS + SPA<sub>SQR</sub>"
nav_order: 3
description: "End-to-end recipe for running SPAsqr with a LOCO polygenic score (without a sparse GRM)."
has_children: false
---

# **Workflow 1: LOCO PGS + SPA<sub>SQR</sub>**

In this section we give a detailed tutorial on how to run SPA<sub>SQR</sub> association testing using a **leave-one-chromosome-out (LOCO) polygenic score (PGS)** as an offset. A LOCO PGS is a per-subject prediction of the trait built from variants on every chromosome *except* the one currently being tested. SPA<sub>SQR</sub> subtracts the chromosome-specific LOCO PGS from the trait before fitting the null smoothed quantile regression model on each chromosome, which helps control for relatedness and substantially improves the statistical power of our score tests.

Although SPA<sub>SQR</sub> is a quantile GWAS method, LOCO PGS computed using *linear* GWAS software works very well for our purpose. We therefore outsource PGS construction to either [**LDAK-KVIK**](https://dougspeed.com/ldak-kvik/) or [**REGENIE**](https://rgcgithub.github.io/regenie/).


## Data that you will need

Throughout this page we use the same set of files:

```
geno.{bed,bim,fam}    PLINK fileset of variants to test
pheno.txt             phenotype file, two traits Y1 and Y2
covar.txt             covariate file with columns covar1 and covar2
```

For concreteness, `pheno.txt` has the following first few rows

```
$ head pheno.txt
FID    IID         Y1     Y2
fam1   sample1   11.8    9.1
fam1   sample2    3.7    0.5
fam2   sample3    9.5    2.2
fam3   sample4    8.4   12.0
```

and `covar.txt` similarly:

```
$ head covar.txt
FID    IID       covar1   covar2
fam1   sample1    0.31     1.04
fam1   sample2   -0.85    -0.22
fam2   sample3    0.07     0.91
fam3   sample4   -1.12     0.45
```

The FID and IID columns of `pheno.txt` and `covar.txt` should match those in `geno.fam`. These simulated data are available in our Github Repository so users can try to replicate this tutorial. Note that GRAB automatically adds an intercept to the covariates, so there is no need to include an intercept column in `covar.txt`.

To utilize the LOCO PGS, we also need a small text file called a **prediction list** — a two-column table pairing each phenotype name (column 1) with the absolute path to its LOCO PGS file (column 2). The prediction list is what GRAB reads (via the `--pred-list` argument) so that it knows which PGS file is paired with which trait. This file format is a REGENIE convention that we follow here.

## INT pre-transformation

Before computing the PGS we recommend applying a **rank-based inverse normal transformation (INT)** to each trait's non-missing values, because in our UK Biobank analysis we find that doing so generally yields more associations than skipping it; this is likely because the LOCO PGS is more accurate when computed on an INT-transformed trait. GRAB offers a utility for applying INT on the phenotype file:

```bash
./grab --int-pheno --pheno pheno.txt --out pheno_int
```

The output `pheno_int.txt` contains the normalized phenotypes, with missing values unchanged:

```
$ head pheno_int.txt
FID    IID         Y1           Y2
fam1   sample1     1.049         0.842
fam1   sample2    -0.299        -1.064
fam2   sample3     0.674         0.299
fam3   sample4     0.522         1.281
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
ldak_step1.step1.pheno1.loco.prs        (LOCO PGS for Y1)
ldak_step1.step1.pheno2.loco.prs        (LOCO PGS for Y2)
```

Within each LOCO PGS file, every row is one subject's LOCO PGS across the 22 autosomes:

```
$ head -3 ldak_step1.step1.pheno1.loco.prs
FID    IID         Chr1     Chr2     Chr3   ...   Chr22
fam1   sample1     0.124   -0.083    0.211         0.196
fam1   sample2    -0.241    0.018   -0.135        -0.057
```

Unlike REGENIE, LDAK-KVIK does not emit a prediction list pairing each phenotype with its LOCO PGS file. One can create `ldak_pred_list.txt` manually with the following shell snippet:

```bash
cat > ldak_pred_list.txt <<EOF
Y1	$(pwd)/ldak_step1.step1.pheno1.loco.prs
Y2	$(pwd)/ldak_step1.step1.pheno2.loco.prs
EOF
```

`$(pwd)` expands to the current working directory so each entry ends up as an absolute path.

GRAB also offers a utility that synthesizes `ldak_pred_list.txt` by scanning the current working directory for files ending in `.loco.prs` and starting with `ldak_step1` and matching them positionally to the columns of `pheno_int.txt`:

```bash
./grab --make-ldak-predlist --prefix ldak_step1 --pheno pheno_int.txt --out ldak_pred_list
```

Provided there are no irrelevant files that also start with `ldak_step1` and end with `.loco.prs`,  the above command writes `ldak_pred_list.txt` with the following content:

```
$ cat ldak_pred_list.txt
Y1	/abs/path/to/ldak_step1.step1.pheno1.loco.prs
Y2	/abs/path/to/ldak_step1.step1.pheno2.loco.prs
```


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
regenie_step1_1.loco          (LOCO PGS for Y1)
regenie_step1_2.loco          (LOCO PGS for Y2)
regenie_step1_pred.list       (pairs Y1 and Y2 with their .loco files)
```

REGENIE's `.loco` format is the **transpose** of LDAK-KVIK's: each row is one chromosome and each column is one subject (with FID and IID joined into a single `FID_IID` token):

```
$ head -3 regenie_step1_1.loco
FID_IID       fam1_sample1   fam1_sample2   fam2_sample3   ...
1             0.0497         -0.0624         -0.0156
2             0.0502         -0.0558         -0.0152
```
Unlike with LDAK-KVIK, REGENIE produces `regenie_step1_pred.list` which can be passed directly to GRAB's `--pred-list` option:

```
$ cat regenie_step1_pred.list
Y1    /abs/path/to/regenie_step1_1.loco
Y2    /abs/path/to/regenie_step1_2.loco
```

GRAB auto-detects the format from the file header and can distinguish whether the LOCO PGS files are produced by LDAK-KVIK or REGENIE.


## Running association testing with GRAB

Once the prediction list is ready, null model fitting and association testing is performed via a single `grab` call:

```bash
./grab --method SPAsqr \
    --bfile geno \
    --pheno pheno_int.txt \
    --covar covar.txt \
    --pred-list ldak_pred_list.txt \
    --out spasqr_results
```

The flags split into a small set you almost always set, plus a wider set you reach for as needed.

**Required:**

| Flag | What it does |
| --- | --- |
| `--method SPAsqr` | Selects the SPA<sub>SQR</sub> method (this is also what triggers all of the `--spasqr-*` options below). |
| `--bfile geno` | PLINK 1 genotype fileset (`geno.{bed,bim,fam}`). PLINK 2 (`--pfile PREFIX`), VCF (`--vcf FILE`), and BGEN (`--bgen FILE`) are also accepted — exactly one of the four is required. |
| `--pheno pheno_int.txt` | Phenotype file. Starts with `FID` and `IID`. |
| `--covar covar.txt` | Covariate file. Starts with `FID` and `IID`. GRAB adds an intercept automatically — do not include one in the file. |
| `--out spasqr_results` | Output prefix. GRAB appends `.<phenoname>.SPAsqr` so each trait gets its own tab-delimited result file. |

**Optional:**

| Flag | Default | What it does |
| --- | --- | --- |
| `--pred-list ldak_pred_list.txt` | — | Prediction list for LDAK-KVIK (or REGENIE's `regenie_step1_pred.list`). Omit to run with no LOCO offset — valid but much less powerful. |
| `--pheno-transform int` | `int` | `int` / `standardize`. **Must match the transform used during PGS construction.** With `pheno_int.txt` and the INT workflow, leave it at the default. With raw `pheno.txt` fed to LDAK-KVIK or REGENIE, set this to `standardize`. |
| `--pheno-name Y1,Y2` | all Y columns | Comma-separated list of trait columns to test. Omit to test every `Y` column found in `--pheno`. |
| `--covar-name covar1,covar2` | all covar columns | Comma-separated list of covariate columns to use. |
| `--spasqr-taus 0.1,0.3,0.5,0.7,0.9` | `0.1,0.3,0.5,0.7,0.9` | Quantile levels at which to test (max 20 levels). |
| `--spasqr-h-scale 3` | `3` (score mode) | Bandwidth divisor: $h = \mathrm{IQR}(\tilde Y - \hat Y_{-c}) / \text{scale}$. Larger $k$ → less smoothing. |

**SNP filters and runtime:**

| Flag | Default | What it does |
| --- | --- | --- |
| `--maf 1e-5` | `1e-5` | Minimum minor allele frequency. |
| `--mac 10` | `10` | Minimum minor allele count. |
| `--geno 0.1` | `0.1` | Maximum per-variant missingness fraction. |
| `--extract snps.txt` | — | Restrict testing to the variant IDs listed in `snps.txt` (one per line). |
| `--chr 1,2,5` | all autosomes | Comma-separated chromosomes to test. |
| `--threads 8` | `1` | Number of threads used for parallel computing. |

Variants that fail the above QC constraints are omitted from GWAS results with $p$-values and $Z$-scores filled by NA.

GRAB writes one output file per phenotype:

```
spasqr_results.Y1.SPAsqr
spasqr_results.Y2.SPAsqr
```

Each file lists per-variant statistics including per-quantile $p$-values, the Cauchy-combined $P_\mathrm{CCT}$, and per-$\tau$ $Z$-scores:

```
$ head -3 spasqr_results.Y1.SPAsqr
CHROM  POS      ID          REF  ALT  MISS_RATE   ALT_FREQ  MAC    HWE_P     P_CCT      P_tau0.1  P_tau0.3   P_tau0.5   P_tau0.7    P_tau0.9    Z_tau0.1    Z_tau0.3   Z_tau0.5  Z_tau0.7  Z_tau0.9
1      1171417  rs6603782   C    T    0.0137558   0.326478  25748  0.536365  0.0071169  0.923124  0.397414   0.0379922  0.00415006  0.00223711  -0.0964998  -0.846248  -2.07494  -2.86651  -3.0567
1      2236359  rs60363208  G    A    0.00435185  0.173843  13841  0.412157  0.227344   0.273908  0.0829882  0.152614   0.509465    0.70161     1.09411     1.73361    1.43036   0.65967   -0.383148
```

The first nine columns are the variant identifier and standard variant-level QC fields. `P_CCT` is the Cauchy-combined $p$-value across all $\tau$ levels and is the main genome-wide significance number to take to a Manhattan plot. The `P_tauX` and `Z_tauX` columns give the per-quantile $p$-values and $Z$-scores; comparing them across $\tau$ reveals whether a hit is driven by a mean shift (uniform signal across $\tau$), a tail effect (concentrated at extreme $\tau$), or a heteroskedastic / dispersion effect (sign change between low and high $\tau$).
### End-to-end recipes (with INT)

With LDAK-KVIK:

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

# 3. Build the pred-list (or write ldak_pred_list.txt by hand)
./grab --make-ldak-predlist --prefix ldak_step1 --pheno pheno_int.txt --out ldak_pred_list

# 4. Run SPAsqr; --pheno-transform int is the default and matches the INT-trained PGS
./grab --method SPAsqr \
    --bfile geno \
    --pheno pheno_int.txt \
    --covar covar.txt \
    --pred-list ldak_pred_list.txt \
    --out spasqr_results
```

With REGENIE (no separate `--make-ldak-predlist` step — REGENIE emits its own prediction list):

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

# 3. Run SPAsqr with REGENIE's native pred-list; --pheno-transform int is the default
./grab --method SPAsqr \
    --bfile geno \
    --pheno pheno_int.txt \
    --covar covar.txt \
    --pred-list regenie_step1_pred.list \
    --out spasqr_results
```

## Skipping the INT pre-transform

You may also skip the INT pre-transform and feed the raw `pheno.txt` directly to LDAK-KVIK or REGENIE. Before fitting the LOCO PGS, both backends internally regress the covariates out of the trait and then standardize the residuals to mean zero and unit variance, so the LOCO PGS still lives on a **standardized scale**, not on the scale of the raw `Y` column. In this case, pass `--pheno-transform standardize` to GRAB at association testing so that the trait GRAB constructs internally and the LOCO PGS it subtracts are on the same scale.

The complete LDAK-KVIK + SPA<sub>SQR</sub> workflow without INT:

```bash
# 1. Train the LOCO PGS on raw Y
./ldak6.2.linux \
    --kvik-step1 ldak_step1 \
    --bfile geno \
    --pheno pheno.txt --mpheno ALL \
    --covar covar.txt \
    --max-threads 8

# 2. Build the pred-list (or write ldak_pred_list.txt by hand)
./grab --make-ldak-predlist --prefix ldak_step1 --pheno pheno.txt --out ldak_pred_list

# 3. Run SPAsqr with --pheno-transform standardize
./grab --method SPAsqr \
    --bfile geno \
    --pheno pheno.txt \
    --covar covar.txt \
    --pred-list ldak_pred_list.txt \
    --pheno-transform standardize \
    --out spasqr_results
```

And the complete REGENIE + SPA<sub>SQR</sub> workflow without INT:

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
