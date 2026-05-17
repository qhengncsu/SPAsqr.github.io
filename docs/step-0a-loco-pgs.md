---
layout: default
title: Step 0a — LOCO polygenic scores
nav_order: 3
description: "Train leave-one-chromosome-out polygenic scores with LDAK-KVIK or REGENIE."
has_children: false
---

# **Step 0a — LOCO polygenic scores**

A **leave-one-chromosome-out (LOCO) polygenic score** soaks up
trans-chromosomal polygenic signal so SPA<sub>SQR</sub>'s null fit on
chromosome *c* doesn't have to model it. Passing the LOCO PGS as an
offset via `--pred-list` noticeably boosts power. It is **optional** —
SPA<sub>SQR</sub> runs without it — but we recommend always including
it.

Two PGS backends are supported, both producing a per-subject table with
one column per autosome that GRAB reads via `--pred-list`:

- [**LDAK-KVIK**](https://dougspeed.com/kvik/) — recommended; fast,
  GRM-aware, no R dependency.
- [**REGENIE**](https://rgcgithub.github.io/regenie/) — alternative;
  widely used.

## The running example

We assume four files in your working directory:

`geno.{bed,bim,fam}` — PLINK fileset of variants to test.

`pheno.txt` — two-trait phenotype file:

```
$ head pheno.txt
FID    IID         Y1     Y2
fam1   sample1   11.8    9.1
fam1   sample2    3.7    0.5
fam2   sample3    9.5    2.2
fam3   sample4    8.4   12.0
...
```

`covar.txt` — covariate file (two covariates plus the implicit
intercept added by GRAB):

```
$ head covar.txt
FID    IID       covar1   covar2
fam1   sample1   0.31     1.04
fam1   sample2  -0.85    -0.22
fam2   sample3   0.07     0.91
fam3   sample4  -1.12     0.45
...
```

By the end of Step 0a you will have one LOCO PGS file per trait and a
small text file pairing each pheno name with the file it should use
when SPA<sub>SQR</sub> reads them back.

---

## Workflow A — INT transform (recommended)

### 1. INT-transform the phenotypes

The GRAB binary ships a utility that inverse-normal-transforms each
non-missing $Y$ column (Blom plotting position, average-rank ties):

```bash
grab --int-pheno --pheno pheno.txt --out pheno_int
```

This writes `pheno_int.txt` with the same FID/IID rows and the same
column names; only the numeric $Y$ entries change:

```
$ head pheno_int.txt
FID    IID         Y1           Y2
fam1   sample1     1.049         0.842
fam1   sample2    -0.299        -1.064
fam2   sample3     0.674         0.299
fam3   sample4     0.522         1.281
...
```

The two traits are transformed independently.

### 2A. LDAK-KVIK + auto-built pred-list

Run KVIK Step 1 on the INT-transformed file:

```bash
ldak6.2.linux \
    --kvik-step1 ldak_step1 \
    --bfile geno \
    --pheno pheno_int.txt --mpheno ALL \
    --covar covar.txt \
    --max-threads 8
```

LDAK writes one LOCO file per trait, indexed by position in
`pheno_int.txt`:

```
ldak_step1.step1.pheno1.loco.prs       # for Y1
ldak_step1.step1.pheno2.loco.prs       # for Y2
```

Each file has one row per subject:

```
$ head -3 ldak_step1.step1.pheno1.loco.prs
FID    IID         Chr1     Chr2     Chr3   ...   Chr22
fam1   sample1     0.124   -0.083    0.211         0.196
fam1   sample2    -0.241    0.018   -0.135        -0.057
...
```

`ChrK` is the polygenic score using all variants **except** chromosome
*K*. Note that the LDAK output uses `pheno1` / `pheno2` (positional),
not `Y1` / `Y2`, so you can't pass the files directly to
`--pred-list` — you need to map each column name to the right file.
GRAB ships a utility that does this for you:

```bash
grab --make-ldak-predlist --pheno pheno_int.txt --out grab_predlist
```

It reads the column names from `pheno_int.txt`, scans the current
directory for LDAK Step 1 outputs matching the trait count, and writes
`grab_predlist.txt`:

```
$ cat grab_predlist.txt
Y1    /abs/path/to/ldak_step1.step1.pheno1.loco.prs
Y2    /abs/path/to/ldak_step1.step1.pheno2.loco.prs
```

If anything is off — no LDAK output in the directory, two ambiguous
LDAK runs side-by-side, duplicate Y column names in the pheno header,
or a path containing whitespace — the utility hard-errors with a
specific message instead of writing a broken file.

### 2B. REGENIE + native pred-list

Run REGENIE Step 1 on the same INT-transformed file:

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

REGENIE produces a per-trait `.loco` file plus its own native pred-list:

```
regenie_step1_1.loco         # for Y1
regenie_step1_2.loco         # for Y2
regenie_step1_pred.list      # pairs each name with its .loco path
```

`regenie_step1_pred.list` is already in the format GRAB wants:

```
$ cat regenie_step1_pred.list
Y1    /abs/path/to/regenie_step1_1.loco
Y2    /abs/path/to/regenie_step1_2.loco
```

so you can pass it directly to `grab --method SPAsqr --pred-list ...` —
no rewriting needed.

### 3. Match `--pheno-transform`

INT is GRAB's default, so when you reach
[Step 1–2]({{ site.baseurl }}/docs/running-spasqr.html) you don't have
to pass `--pheno-transform` explicitly. You may also feed
`pheno.txt` (raw) instead of `pheno_int.txt` to SPA<sub>SQR</sub> — GRAB
will INT internally and produce identical results.

---

## Workflow B — standardize (skip the INT step)

If you don't want to pre-INT, you can give either LDAK-KVIK or REGENIE
your **raw** phenotype directly and ask GRAB to standardize on the same
scale. Both backends *already* standardize internally before building
the PGS, so the LOCO PGS comes out on the standardized scale of your
raw $Y$ either way.

### 1. Build the LOCO PGS on raw `pheno.txt`

LDAK-KVIK:

```bash
ldak6.2.linux \
    --kvik-step1 ldak_step1 \
    --bfile geno \
    --pheno pheno.txt --mpheno ALL \
    --covar covar.txt \
    --max-threads 8
grab --make-ldak-predlist --pheno pheno.txt --out grab_predlist
```

REGENIE:

```bash
regenie \
    --step 1 \
    --bed geno \
    --phenoFile pheno.txt --phenoColList Y1,Y2 \
    --covarFile covar.txt --covarColList covar1,covar2 \
    --bsize 1000 \
    --lowmem --lowmem-prefix tmp_rg \
    --out regenie_step1
```

### 2. Match `--pheno-transform`

In Step 1–2 you pass `--pheno-transform standardize`. GRAB then
centres and unit-scales the supplied $Y$ to match the scale of the
LOCO PGS the backend produced.

> **Why this works without explicit pre-standardization.** Both
> LDAK-KVIK and REGENIE re-centre and unit-scale the trait you give
> them before fitting the PGS, regardless of its incoming distribution.
> Pre-standardizing `pheno.txt` would be a no-op for the PGS; the
> useful pre-transform is INT (Workflow A), which changes the
> *shape* of the distribution, not just its scale.

---

## When to skip Step 0a

If your cohort is too small to train a useful PGS (say,
*n* ≲ 10⁴), or if you just want a quick first pass, omit `--pred-list`
when calling `grab --method SPAsqr`. SPA<sub>SQR</sub> will fit the
null SQR model with no offset — you lose the polygenic power gain but
everything else still works, and the `--pheno-transform` consistency
rule no longer constrains you.

> **Threading caveat for LDAK-KVIK.** `--max-threads` should not exceed
> the number of physical CPU cores. Oversubscription causes thread
> contention and can slow Step 1 down by an order of magnitude.

You'll feed the resulting `grab_predlist.txt` (or
`regenie_step1_pred.list`) to SPA<sub>SQR</sub> via `--pred-list` in
[Step 1–2]({{ site.baseurl }}/docs/running-spasqr.html).
