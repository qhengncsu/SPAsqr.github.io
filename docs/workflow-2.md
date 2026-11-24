---
layout: default
title: "Workflow 2: LOCO PGS + GRM + SPA<sub>SQR</sub>"
nav_order: 4
description: "End-to-end recipe for running SPAsqr with a LOCO PGS and a sparse GRM."
has_children: false
---

# **Workflow 2: LOCO PGS + GRM + SPA<sub>SQR</sub>**

In cohorts with strong relatedness, a LOCO PGS offset alone may leave the tests poorly calibrated. REGENIE's LOCO PGS is usually less predictive than LDAK-KVIK's, so relatedness can still inflate the type I error rate. With LDAK-KVIK's LOCO PGS, the adjustment can go too far and deflate the test statistics.

In this workflow, we add a **sparse genetic relationship matrix (GRM)** to account for correlation between related participants. SPA<sub>SQR</sub> uses the GRM only when estimating the variance of the score statistic. This improves calibration under strong relatedness. For admixed cohorts, we do not recommend adding the GRM because population structure can confound estimates of relatedness.

The example below builds a sparse GRM with PLINK 2 and adds it to [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html).

## Complete pipeline

Follow steps 1–3 of Workflow 1 to prepare the LOCO PGS. In steps 4 and 5, build the sparse GRM and include it in the association test.

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

PLINK 2 has supported sparse GRMs since late 2025. GCTA can compute the same matrix, though it takes longer. We use `--maf 0.01` to build the GRM from common variants. The option `--make-grm-sparse 0.05` retains relatedness coefficients above 0.05 and sets the rest to zero.

PLINK 2 writes the matrix to `simu_geno.grm.sp` and the participant IDs to `simu_geno.grm.id`. Each row of `.grm.sp` contains two participant indices and their relatedness coefficient. The indices start at zero, so index $i$ refers to row $i+1$ in `.grm.id`.

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
    --threads 8 \
    --out spasqr_results
```

Add `--sp-grm-plink2 simu_geno.grm.sp` to the command from Workflow 1. GRAB uses the shared file prefix to find `simu_geno.grm.id` and writes results in the same format as in Workflow 1.

## Other options

### REGENIE instead of LDAK-KVIK

If you train the LOCO PGS with REGENIE as described in Workflow 1, use `--pred-list simu_geno_regenie_pred.list` in step 5. This is the only change you need to make.

### GRM from other tools: `--sp-grm-grab`

If you used another tool to build the GRM, you can supply a single text file with participant IDs using `--sp-grm-grab` in place of `--sp-grm-plink2`.

```
$ head simu_geno.grm.grab
IID1    IID2     VALUE
IID_0   IID_0   1.0024
IID_1   IID_1   0.9981
IID_3   IID_0   0.5012
IID_4   IID_3   0.2503
```

Use a tab-delimited file with the header `IID1 IID2 VALUE`. The IDs must match those in the `.fam` file, so a separate `.grm.id` file is unnecessary. Include the diagonal entries and list each related pair once. GRAB fills in the symmetric entries and treats unlisted pairs as zero.

```bash
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --sp-grm-grab simu_geno.grm.grab \
    --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
    --threads 8 \
    --out spasqr_results
```
