---
layout: default
title: "Workflow 2: LOCO PGS + GRM + SPA<sub>SQR</sub>"
nav_order: 4
description: "End-to-end recipe for running SPAsqr with a LOCO polygenic score and a sparse GRM."
has_children: false
---

# **Workflow 2: LOCO PGS + GRM + SPA<sub>SQR</sub>**

In highly related cohorts the LOCO PGS offset alone may not calibrate the tests: REGENIE's PGS may leave residual inflation, and LDAK-KVIK's may deflate. A **sparse genetic relationship matrix (GRM)** calibrates the null variance of the score statistic under strong relatedness. This page adds one to [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html).

## Complete pipeline

Steps 1–3 are identical to Workflow 1; steps 4–5 are new.

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

# 4. Build the sparse GRM with PLINK 2
./plink2 \
    --bfile simu_geno \
    --maf 0.01 \
    --make-grm-sparse 0.05 \
    --threads 8 \
    --out simu_geno

# 5. Run SPAsqr with the LOCO PGS and the sparse GRM
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --sp-grm-plink2 simu_geno.grm.sp \
    --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
    --pheno-transform int \
    --threads 8 \
    --out spasqr_results
```

## Step by step

### 4. Build the sparse GRM

```bash
./plink2 \
    --bfile simu_geno \
    --maf 0.01 \
    --make-grm-sparse 0.05 \
    --threads 8 \
    --out simu_geno
```

- `--maf 0.01` uses only common variants.
- `--make-grm-sparse 0.05` keeps relatedness coefficients above 0.05 and zeroes the rest. (Available in PLINK 2 since late 2025; GCTA computes the same GRM but more slowly.)

Outputs `simu_geno.grm.sp` and its companion `simu_geno.grm.id`. Indices are **0-based** rows of `.grm.id`:

```
$ head simu_geno.grm.sp
0   0   1.0024
1   1   0.9981
2   2   1.0107
3   0   0.5012      # half-sibling of subject 0
4   3   0.2503      # first cousin of subject 3
```

### 5. Run SPA<sub>SQR</sub> with the GRM

```bash
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --sp-grm-plink2 simu_geno.grm.sp \
    --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
    --pheno-transform int \
    --threads 8 \
    --out spasqr_results
```

The only change from Workflow 1 is `--sp-grm-plink2 simu_geno.grm.sp`. GRAB finds `simu_geno.grm.id` from the same prefix. Output format is the same as Workflow 1.

## Other options

### REGENIE instead of LDAK-KVIK

Train the PGS with REGENIE as in Workflow 1, then pass `--pred-list simu_geno_regenie_pred.list` in step 5. Everything else is unchanged.

### GRM from other tools: `--sp-grm-grab`

For a GRM not computed by PLINK 2, write it as a single IID-keyed text file and pass it with `--sp-grm-grab` instead of `--sp-grm-plink2`:

```
$ head simu_geno.grm.grab
IID1    IID2     VALUE
IID_0   IID_0   1.0024
IID_1   IID_1   0.9981
IID_3   IID_0   0.5012
IID_4   IID_3   0.2503
```

- Tab-delimited, header `IID1 IID2 VALUE`; IIDs match the `.fam` file. No `.grm.id` needed.
- One row per related pair (GRAB symmetrizes) plus the diagonal. Unlisted pairs are zero.

```bash
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --sp-grm-grab simu_geno.grm.grab \
    --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
    --pheno-transform int \
    --threads 8 \
    --out spasqr_results
```

### When to omit the GRM

The GRM only changes the reference distribution, not the score statistic. For cohorts with low relatedness, omit it (Workflow 1) and the results will be very similar.

A GCTA-style GRM is unreliable for admixed or multi-ancestry cohorts, where population structure produces far too many entries above 0.05.

- **Admixed, low relatedness:** omit the GRM (Workflow 1).
- **Admixed and highly related** (e.g. the Mexico City Prospective Study): compute an ancestry-aware sparse GRM with [FastSparseGRM](https://github.com/rounakdey/FastSparseGRM), write it in the `--sp-grm-grab` format above, and pass it with `--sp-grm-grab`.
