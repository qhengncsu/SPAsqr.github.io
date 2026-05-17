---
layout: default
title: Step 0a — LOCO polygenic scores
nav_order: 3
description: "Train leave-one-chromosome-out polygenic scores with LDAK-KVIK or REGENIE, then pair them with GRAB via a pred-list."
has_children: false
---

# **Step 0a — LOCO polygenic scores**

A **leave-one-chromosome-out (LOCO) polygenic score (PGS)** is a per-subject prediction of the trait built from variants on every chromosome *except* the one currently being tested. SPA<sub>SQR</sub> uses LOCO PGS as an offset before fitting the null smoothed quantile regression model on each chromosome. This helps us better control for relatedness and also substantially improves the statistical power of our score tests. Although SPA<sub>SQR</sub> is a quantile GWAS method, LOCO PGS computed using linear GWAS software such as [**LDAK-KVIK**](https://dougspeed.com/ldak-kvik/) and [**REGENIE**](https://rgcgithub.github.io/regenie/) work very well for our purpose. Thus, SPA<sub>SQR</sub> outsources the complex task of PGS construction to existing software.

On this page we give a detailed tutorial on how to compute LOCO PGS using either LDAK-KVIK or REGENIE.

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

By the end of Step 0a we will have one LOCO PGS file per trait plus a small text file called a **pred-list** — a two-column table pairing each phenotype name (column 1) with the absolute path to its LOCO PGS file (column 2). The pred-list is what GRAB reads (via the `--pred-list` argument) to know which PGS to subtract from which trait. This two-column format is a REGENIE convention; GRAB adopts the same convention so that REGENIE's native output (`regenie_step1_pred.list`) can be passed in unchanged, while on the LDAK-KVIK side we build the same file ourselves and conventionally call it `ldak_pred_list.txt`.

## INT pre-transformation

Before computing the PGS we recommend applying a **rank-based inverse normal transformation (INT)** to each trait's non-missing values, because in our UK Biobank analysis we find that doing so generally yields more associations than skipping it; this is likely because the LOCO PGS is more accurate when computed on an INT-transformed trait. GRAB offers a utility for applying INT:

```bash
grab --int-pheno --pheno pheno.txt --out pheno_int
```

The output `pheno_int.txt` contains the normalized phenotypes:

```
$ head pheno_int.txt
FID    IID         Y1           Y2
fam1   sample1     1.049         0.842
fam1   sample2    -0.299        -1.064
fam2   sample3     0.674         0.299
fam3   sample4     0.522         1.281
```

## Computing the LOCO PGS with LDAK-KVIK

We can compute the LOCO PGS via LDAK-KVIK step 1 with:

```bash
ldak6.2.linux \
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

Unlike REGENIE, LDAK-KVIK does not emit a pred-list pairing each phenotype with its LOCO PGS file. One can create `ldak_pred_list.txt` manually with a short shell snippet:

```bash
cat > ldak_pred_list.txt <<EOF
Y1    $(pwd)/ldak_step1.step1.pheno1.loco.prs
Y2    $(pwd)/ldak_step1.step1.pheno2.loco.prs
EOF
```

`$(pwd)` expands to the current working directory so each entry ends up as an absolute path, which is what GRAB requires.

As a shortcut for the manual snippet above, GRAB also offers a utility that synthesizes `ldak_pred_list.txt` by scanning the current working directory for files ending in `.loco.prs` and matching them positionally to the columns of `pheno_int.txt`:

```bash
grab --make-ldak-predlist --pheno pheno_int.txt --out ldak_pred_list
```

Under ideal conditions, the command writes `ldak_pred_list.txt` with the following content:

```
$ cat ldak_pred_list.txt
Y1    /abs/path/to/ldak_step1.step1.pheno1.loco.prs
Y2    /abs/path/to/ldak_step1.step1.pheno2.loco.prs
```

**But `grab --make-ldak-predlist` may fail** — for example, when multiple LDAK runs share the same working directory and the match is ambiguous, or when not all expected `phenoN` files are present — the utility hard-errors with an explicit message, and it is recommended to fall back to the manual approach.

## Computing the LOCO PGS with REGENIE

We may also compute LOCO PGS via REGENIE Step 1:

```bash
regenie \
    --step 1 \
    --bed geno \
    --phenoFile pheno_int.txt \
    --covarFile covar.txt \
    --bsize 1000 \
    --out regenie_step1
```

This produces one LOCO PGS file per trait plus REGENIE's own pairing list:

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

This layout difference is transparent to the user: GRAB auto-detects the format from the file header and parses both LDAK-KVIK and REGENIE outputs uniformly.

Unlike with LDAK-KVIK, the `regenie_step1_pred.list` REGENIE produces can be passed directly to GRAB's `--pred-list` option — one row per phenotype, with the phenotype name in column 1 and the absolute path to the corresponding `.loco` file in column 2:

```
$ cat regenie_step1_pred.list
Y1    /abs/path/to/regenie_step1_1.loco
Y2    /abs/path/to/regenie_step1_2.loco
```

This two-column "name + absolute path" pred-list format is a REGENIE convention. GRAB adopts the same convention so REGENIE Step 1 output can be piped straight in without rewriting, and `grab --make-ldak-predlist` produces a file in the same format on the LDAK-KVIK side.

## Skipping the INT pre-transform

If you prefer to skip the INT transform, we can instead feed the raw `pheno.txt` directly to either backend. For example, for LDAK-KVIK, we now run

```bash
ldak6.2.linux \
      --kvik-step1 ldak_step1 \
      --bfile geno \
      --pheno pheno.txt --mpheno ALL \
      --covar covar.txt \
      --max-threads 8
grab --make-ldak-predlist --pheno pheno.txt --out ldak_pred_list
```

Before fitting the LOCO PGS, both LDAK-KVIK and REGENIE internally regress the covariates out of the trait and then standardize the residuals to mean zero and unit variance. The resulting LOCO PGS therefore lives on a **standardized scale**, not on the scale of the raw `Y` column.

To keep the LOCO PGS offset on the same scale as the trait it is subtracted from, pass `--pheno-transform standardize` to GRAB in Step 1–2; GRAB then centres and unit-scales the values in `pheno.txt` internally before subtracting the LOCO PGS.

> **Note on the `--pheno-transform` default.** The default value of `--pheno-transform` is `int`, matching the recommended INT workflow above. If you fed `pheno_int.txt` to LDAK-KVIK or REGENIE, you can omit `--pheno-transform` entirely when invoking GRAB. You only need to set `--pheno-transform standardize` (or `raw`) when you deliberately deviate from the default INT workflow.
