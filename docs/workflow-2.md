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

Steps 1–3 are identical to Workflow 1. Steps 4 and 5 are new.

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

We compute the sparse GRM with PLINK 2, which has supported this since late 2025 (GCTA computes the same GRM, but more slowly). `--maf 0.01` restricts the computation to common variants, and `--make-grm-sparse 0.05` keeps only relatedness coefficients above 0.05 and zeroes the rest.

PLINK 2 writes `simu_geno.grm.sp` together with the companion ID file `simu_geno.grm.id`. Each row of `.grm.sp` holds a pair of **0-based** subject indices and their relatedness coefficient, where index $i$ refers to the $(i+1)$-th row of `.grm.id`:

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

The only change from Workflow 1 is the added `--sp-grm-plink2 simu_geno.grm.sp`. GRAB finds `simu_geno.grm.id` from the same prefix. The output format is the same as in Workflow 1.

## Other options

### REGENIE instead of LDAK-KVIK

Train the LOCO PGS with REGENIE as described in Workflow 1, then pass `--pred-list simu_geno_regenie_pred.list` in step 5. Everything else stays the same.

### GRM from other tools: `--sp-grm-grab`

If the GRM was computed by a tool other than PLINK 2, write it as a single IID-keyed text file and pass it with `--sp-grm-grab` instead of `--sp-grm-plink2`. The file looks like this:

```
$ head simu_geno.grm.grab
IID1    IID2     VALUE
IID_0   IID_0   1.0024
IID_1   IID_1   0.9981
IID_3   IID_0   0.5012
IID_4   IID_3   0.2503
```

The file is tab-delimited with the header `IID1 IID2 VALUE`, where the IIDs match the `.fam` file, so no companion `.grm.id` is needed. It lists each related pair once (GRAB symmetrizes the matrix) plus the diagonal entries; any pair not listed is treated as zero.

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
