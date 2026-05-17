---
layout: default
title: Step 0b — Sparse GRM
nav_order: 4
description: "Construct a sparse genetic relationship matrix for variance calibration."
has_children: false
---

# **Step 0b — Sparse GRM**

If the cohort contains a non-trivial number of related individuals,
SPA<sub>SQR</sub> calibrates the variance of its score statistic using
a sparse genetic relationship matrix $\Phi$:

$$
\widehat{\mathrm{Var}}(S_j) \;=\;
\widehat{\sigma}_g^{\,2}(G_j)\, R^{\!\top}\,\Phi\, R.
$$

When the data are essentially unrelated, $\Phi = I_n$ and the formula
reduces to $\widehat\sigma_g^{\,2}(G_j)\, R^{\!\top} R$ — in which
case the GRM input is **not required**. Supplying a sparse GRM only
matters when off-diagonal kinship is non-negligible (close relatives,
clinical pedigrees, founder populations).

## Recommended path — PLINK 2 `--make-grm-sparse`

Since late 2025, [PLINK 2](https://www.cog-genomics.org/plink/2.0/)
supports sparse GRM construction natively. It is substantially faster
than GCTA and ships as a single binary:

```bash
plink2 \
    --bfile geno \
    --maf 0.01 \
    --make-grm-sparse 0.05 \
    --threads 8 \
    --out geno_sparse
```

- `--maf 0.01`: restrict GRM computation to variants with minor allele
  frequency $\geq 0.01$. Rare-variant noise inflates GRM elements; the
  MAF filter is essential.
- `--make-grm-sparse 0.05`: retain genetic correlation coefficients
  above $0.05$. Anything smaller is zeroed out (sparsified). The
  cutoff trades off coverage of distant cousins (smaller cutoff →
  denser matrix, more cousins captured) against computational cost.
  $0.05$ corresponds roughly to 3rd-degree relatives and is the
  default we recommend.

The command produces two files:

```
geno_sparse.grm.sp   # three columns: 0-based i, 0-based j, GRM[i,j]
geno_sparse.grm.id   # one row per subject: FID  IID
```

Pass `geno_sparse.grm.sp` to GRAB via `--sp-grm-plink2`. GRAB
auto-detects the companion `.grm.id` file based on the prefix.

> **Pre-pruning.** For populations with strong LD or admixed cohorts,
> consider LD-pruning the input PLINK fileset before
> `--make-grm-sparse` to reduce GRM bias:
>
> ```bash
> plink2 --bfile geno --indep-pairwise 500 50 0.2 --out prune
> plink2 --bfile geno --extract prune.prune.in --maf 0.01 \
>        --make-grm-sparse 0.05 --threads 8 --out geno_sparse
> ```

## Legacy path — GCTA

GCTA's sparse GRM works the same way:

```bash
gcta64 --bfile geno --maf 0.01 --make-grm --threads 8 --out geno_grm
gcta64 --grm geno_grm --make-bK-sparse 0.05 --out geno_sparse
```

GCTA's `.grm.sp` output has the **same three-column schema** as PLINK
2's; you can pass it to GRAB unchanged via
`--sp-grm-plink2 geno_sparse.grm.sp`. GCTA is significantly slower than
PLINK 2 on large cohorts, so we recommend the PLINK 2 path unless you
already have a GCTA-built GRM.

## When to skip Step 0b

You can omit the sparse GRM in any of these situations:

- The cohort is genuinely unrelated (e.g. independent biobank samples
  with no close relatives, KING relatedness < 1st-degree).
- A first-pass / sanity-check run while you iterate on phenotype QC.

When `--sp-grm-plink2` is omitted, SPA<sub>SQR</sub> falls back to the
unrelated-samples variance $\widehat\sigma_g^{\,2}(G_j)\, R^{\!\top}R$.
The score test stays valid, but tail $p$-values for rare variants in
related samples will be miscalibrated.

## What the `.grm.sp` file looks like

```
$ head geno_sparse.grm.sp
0     0    1.0024
1     1    0.9981
2     2    1.0107
3     0    0.5012      # half-sibling of subject 0
4     3    0.2503      # cousin of subject 3
...
```

Sample $i$/$j$ indices are **0-based** and correspond to row order in
the companion `geno_sparse.grm.id`. Diagonal entries are the
self-kinship $1 + F$ (inbreeding coefficient $F$ is typically near 0).

## Note on Y / X chromosomes

`--make-grm-sparse` works on autosomes by default. If you want chrX or
chrY contributions in the GRM, restrict the input to those
chromosomes with `--chr 23` / `--chr 24` and build a separate matrix.
Most analyses use autosomal GRMs only.

> **Note**
> - The sparse GRM is a one-time cost per cohort — you do not rebuild
>   it for every trait.
> - The same `.grm.sp` is reused across `--pred-list`,
>   `--pheno-transform`, and `--spasqr-mode` choices.
