---
layout: default
title: Step 0a — LOCO polygenic scores
nav_order: 3
description: "Train leave-one-chromosome-out polygenic scores with LDAK-KVIK or REGENIE, then pair them with GRAB via a pred-list."
has_children: false
---

# **Step 0a — LOCO polygenic scores**

SPA<sub>SQR</sub> subtracts a **leave-one-chromosome-out (LOCO)
polygenic score** from the phenotype before fitting the null smoothed
quantile regression on each chromosome. The offset absorbs the
trans-chromosomal polygenic component of $Y$, sharpens the rank-score
residual, and translates directly into a power gain that grows with
trait heritability. The whole step is *optional* — when `--pred-list`
is omitted, SPA<sub>SQR</sub> fits the null model with no offset — but
in any cohort large enough to train a useful PGS (roughly $n \gtrsim
10^4$) we recommend running it.

Two backends are supported, both producing a per-subject table whose
columns are chromosome-indexed leave-one-chromosome-out predictions:

[**LDAK-KVIK**](https://dougspeed.com/kvik/) — recommended; fast,
GRM-aware, no R dependency.
[**REGENIE**](https://rgcgithub.github.io/regenie/) — alternative;
widely used, well-documented.

The two file formats differ in layout (LDAK writes one row per subject
with chromosome columns; REGENIE writes one row per chromosome with
subject columns), and GRAB auto-detects which format it is parsing
from the first header token.

## A running example

Throughout this page we use the same set of files:

```
geno.{bed,bim,fam}    PLINK fileset of variants to test
pheno.txt             phenotype file, two traits Y1 and Y2
covar.txt             covariate file with columns covar1 and covar2
```

For concreteness `pheno.txt` begins

```
$ head pheno.txt
FID    IID         Y1     Y2
fam1   sample1   11.8    9.1
fam1   sample2    3.7    0.5
fam2   sample3    9.5    2.2
fam3   sample4    8.4   12.0
```

and `covar.txt` similarly. By the end of Step 0a we will have one LOCO
PGS file per trait plus a small text file pairing each trait name with
the file it should use; that pairing is what SPA<sub>SQR</sub> reads
back via `--pred-list` in Step 1–2.

## INT pre-transformation

Before training the PGS we apply a **rank-based inverse normal
transformation (INT)** to each trait's non-missing values. This step
removes scale and skew from $Y$ without altering its ordering, so
quantile-regression score statistics on the post-transform residual
have well-behaved tails — which matters most for the rare-variant end
of the genome. GRAB ships a one-shot utility:

```bash
grab --int-pheno --pheno pheno.txt --out pheno_int
```

The output `pheno_int.txt` preserves the original row order and column
names; only the numeric entries change, with each $Y$ column INT'd
independently against its own non-missing values:

```
$ head pheno_int.txt
FID    IID         Y1           Y2
fam1   sample1     1.049         0.842
fam1   sample2    -0.299        -1.064
fam2   sample3     0.674         0.299
fam3   sample4     0.522         1.281
```

We recommend passing `pheno_int.txt` to whichever PGS backend you
choose; the rationale is that LDAK-KVIK and REGENIE both internally
*standardize* (centre and unit-scale) the trait they receive, but
neither performs an INT. Pre-transforming with INT yields a LOCO PGS
that lives on the INT scale, so the SPA<sub>SQR</sub> offset is on
the same scale as the response GRAB constructs internally.

## Training the LOCO PGS with LDAK-KVIK

LDAK-KVIK Step 1 fits a sparse elastic-net polygenic model and emits a
leave-one-chromosome-out prediction for every autosome:

```bash
ldak6.2.linux \
    --kvik-step1 ldak_step1 \
    --bfile geno \
    --pheno pheno_int.txt --mpheno ALL \
    --covar covar.txt \
    --max-threads 8
```

With `--mpheno ALL`, LDAK fits each $Y$ column independently and writes
one LOCO PGS file per trait, indexed by the trait's **position** in
the phenotype file rather than its column name:

```
ldak_step1.step1.pheno1.loco.prs        (paired with Y1)
ldak_step1.step1.pheno2.loco.prs        (paired with Y2)
```

Each row corresponds to one subject and each chromosome column is the
score predicted using all variants *except* those on that chromosome:

```
$ head -3 ldak_step1.step1.pheno1.loco.prs
FID    IID         Chr1     Chr2     Chr3   ...   Chr22
fam1   sample1     0.124   -0.083    0.211         0.196
fam1   sample2    -0.241    0.018   -0.135        -0.057
```

We caution that `--max-threads` should not exceed the number of
physical CPU cores: oversubscription causes thread contention and can
slow Step 1 down by an order of magnitude. We also caution that
LDAK-KVIK silently mean-imputes missing covariates before training. If
that does not match how you intend to handle missingness, impute the
covariate file yourself before invoking LDAK.

### Building the pred-list automatically

Because LDAK names its outputs positionally (`pheno1`, `pheno2`, …)
rather than by trait name, you cannot hand `--pred-list` a directory
listing — GRAB needs to know that `pheno1.loco.prs` corresponds to
column `Y1`. The GRAB binary ships a second utility for this:

```bash
grab --make-ldak-predlist --pheno pheno_int.txt --out grab_predlist
```

The command reads the column names from the supplied phenotype file,
scans the current directory for `*.step1.loco.prs` and
`*.step1.phenoN.loco.prs` files, and writes the pairing to
`grab_predlist.txt`:

```
$ cat grab_predlist.txt
Y1    /abs/path/to/ldak_step1.step1.pheno1.loco.prs
Y2    /abs/path/to/ldak_step1.step1.pheno2.loco.prs
```

The utility is strict: it refuses to write a pred-list whenever the
match is ambiguous (two LDAK runs side-by-side covering the same trait
count), incomplete (some `phenoN` files missing for the requested $K$
traits), or syntactically dangerous (the phenotype header contains
duplicate column names, or a candidate LOCO path contains whitespace
that the downstream pred-list parser would split on). The user is
told explicitly what went wrong rather than silently inheriting a
broken pairing.

## Training the LOCO PGS with REGENIE

REGENIE Step 1 fits a ridge-regression LOCO PGS via a two-level stacked
regression:

```bash
regenie \
    --step 1 \
    --bed geno \
    --phenoFile pheno_int.txt --phenoColList Y1,Y2 \
    --covarFile covar.txt --covarColList covar1,covar2 \
    --bsize 1000 \
    --lowmem --lowmem-prefix tmp_rg \
    --out regenie_step1
```

This produces one per-trait LOCO file plus REGENIE's own pairing list:

```
regenie_step1_1.loco         (paired with Y1; row = chromosome, column = subject)
regenie_step1_2.loco         (paired with Y2)
regenie_step1_pred.list      (Y1 → file 1, Y2 → file 2)
```

The `_pred.list` file REGENIE writes is already in the format GRAB
expects, so no auxiliary command is needed — you pass it directly to
`grab --method SPAsqr --pred-list regenie_step1_pred.list` in the next
step. (REGENIE files are read out by chromosome rather than by
subject, but GRAB detects the layout from the file header and parses
either format transparently.)

## Skipping the INT pre-transform

Both LDAK-KVIK and REGENIE re-centre and unit-scale every trait they
receive, regardless of its incoming distribution. As a consequence the
LOCO PGS each tool emits always lives on a standardized scale, and a
viable alternative to Workflow A is to skip `grab --int-pheno`
entirely and rely on the backend's internal standardization:

```bash
ldak6.2.linux                                  \
    --kvik-step1 ldak_step1 --bfile geno       \
    --pheno pheno.txt --mpheno ALL             \
    --covar covar.txt --max-threads 8
grab --make-ldak-predlist --pheno pheno.txt --out grab_predlist
```

or, equivalently, the analogous REGENIE invocation against
`pheno.txt`. In Step 1–2 you then pass `--pheno-transform standardize`
so that GRAB applies the same centring-and-scaling internally before
subtracting the offset. Explicitly pre-standardizing `pheno.txt` is a
no-op for the PGS — both backends would re-standardize the
already-standardized values — and so we omit that step.

The two workflows are statistically distinct: INT changes the *shape*
of the trait distribution (skewed traits become approximately
Gaussian) and is what we recommend for genome-wide screening,
particularly when rare variants are tested. Standardize alone changes
only location and scale and leaves the long tail intact. The remainder
of the documentation assumes INT unless stated otherwise.
