---
layout: default
title: Step 0a — LOCO polygenic scores
nav_order: 3
description: "Train leave-one-chromosome-out polygenic scores with LDAK-KVIK or REGENIE."
has_children: false
---

# **Step 0a — LOCO polygenic scores**

SPA<sub>SQR</sub> takes a chromosome-specific **leave-one-chromosome-out
(LOCO) polygenic score** $\hat Y_{-c}$ as an *offset* in the null SQR
fit on chromosome $c$. The LOCO PGS soaks up trans-chromosomal
polygenic signal, sharpens the rank-score residual, and noticeably
boosts power. It is **optional** — SPA<sub>SQR</sub> runs without it
— but we recommend always including it.

We support two PGS backends, both producing a one-row-per-subject
table with chromosome-indexed columns that GRAB reads via
`--pred-list`:

- [**LDAK-KVIK**](https://dougspeed.com/kvik/) — recommended; fast,
  GRM-aware, no external R dependencies.
- [**REGENIE**](https://rgcgithub.github.io/regenie/) — alternative;
  widely used, well-documented.

## Phenotype preprocessing — INT transform

Before training the PGS, we apply a **rank-based inverse normal
transformation (INT)** to the non-missing values of each trait. The
GRAB binary ships a utility for this:

```bash
grab --int-pheno --pheno pheno.txt --out pheno_int
```

The input `pheno.txt` is whitespace-separated with a header line
`FID IID Y1 Y2 …`; missing entries may be `NA`, `.`, or blank. The
output `pheno_int.txt` has the same columns and row order, with each
$Y$ column independently INT-transformed (Blom plotting position,
average-rank ties) on its own non-missing scope. Missing entries stay
missing.

```
$ head pheno.txt
FID    IID        Y1     Y2
fam1   sample1   11.8   9.1
fam1   sample2    3.7   0.5
fam2   sample3     NA    2.2
...
$ head pheno_int.txt
FID    IID        Y1           Y2
fam1   sample1    1.049         0.842
fam1   sample2   -0.299        -1.064
fam2   sample3    NA            0.299
...
```

> **Why INT-transform first?** LDAK-KVIK and REGENIE internally
> *standardize* the trait before PGS construction; they do not INT.
> Pre-transforming yields PGS that live on the INT scale, which (i)
> stabilizes rare-variant tail behavior in Step 0a and (ii) lets you
> later pass `--pheno-transform int` (the GRAB default) so the SPA<sub>SQR</sub>
> offset is on the same scale as the response. See the
> [`--pheno-transform` discussion]({{ site.baseurl }}/docs/running-spasqr.html#pheno-transform)
> for the full consistency rule.

## Option A — LDAK-KVIK (recommended)

Download the LDAK 6.2 binary from
[dougspeed.com/kvik](https://dougspeed.com/kvik/). Then run KVIK
Step 1 on the INT-transformed phenotype:

```bash
ldak6.2.linux \
    --kvik-step1 ldak_step1 \
    --bfile geno \
    --pheno pheno_int.txt --mpheno ALL \
    --covar covar.txt \
    --max-threads 8
```

The command produces one LOCO PGS file per trait:

```
ldak_step1.step1.Y1.loco.prs
ldak_step1.step1.Y2.loco.prs
```

Each file has one row per subject with the schema

```
FID    IID    Chr1    Chr2    ...    Chr22
```

where `ChrK` is the polygenic score predicted using all variants
*except* those on chromosome `K`. GRAB reads these columns directly
via `--pred-list`.

> **Threading caveat.** `--max-threads` should not exceed the number
> of physical CPU cores. Oversubscription causes thread contention and
> can slow Step 1 down by an order of magnitude.

## Option B — REGENIE

Download REGENIE from
[rgcgithub.github.io/regenie](https://rgcgithub.github.io/regenie/).
Then run REGENIE Step 1 on the INT-transformed phenotype:

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

The command produces `regenie_step1.list` plus per-trait LOCO blup
files
```
regenie_step1_1.loco
regenie_step1_2.loco
```
each with one row per subject and chromosome columns analogous to the
LDAK output. REGENIE expects an integer trait index in
`regenie_step1.list`, not a name; GRAB resolves this when you provide
the `--pred-list` (see below).

## The `--pred-list` file

Whichever backend you used, you write a small text file telling GRAB
which phenotype is paired with which LOCO PGS file (one row per
phenotype):

```
Y1    /path/to/ldak_step1.step1.Y1.loco.prs
Y2    /path/to/ldak_step1.step1.Y2.loco.prs
```

The pheno name in column 1 **must** match the column name in your
phenotype file (and the `--pheno-name` value you pass to
`grab --method SPAsqr`). Column 2 is the absolute path to the LOCO
PGS file. The same format works for both LDAK-KVIK and REGENIE
outputs.

You'll pass this file to SPA<sub>SQR</sub> via `--pred-list` in
[Step 1–2]({{ site.baseurl }}/docs/running-spasqr.html).

## When to skip Step 0a

If your cohort is too small to train a useful PGS (say,
$n \lesssim 10^4$), or if you only want a quick first-pass run,
**omit** `--pred-list` and SPA<sub>SQR</sub> will fit the null SQR
model with no offset. You forfeit the polygenic power gain but
everything else still works.

> **Note**
> - The same `pheno_int.txt` should be passed to both Step 0a *and*
>   `grab --method SPAsqr`. GRAB will apply
>   `--pheno-transform int` internally, which is idempotent on an
>   already-INT'd column — so passing either `pheno.txt` or
>   `pheno_int.txt` to SPA<sub>SQR</sub> gives the same result, provided the
>   chosen transform matches what the PGS was trained on.
