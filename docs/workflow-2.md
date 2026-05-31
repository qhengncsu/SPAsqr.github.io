---
layout: default
title: "Workflow 2: LOCO PGS + GRM + SPA<sub>SQR</sub>"
nav_order: 4
description: "End-to-end recipe for running SPAsqr with a LOCO polygenic score and a sparse GRM."
has_children: false
---

# **Workflow 2: LOCO PGS + GRM + SPA<sub>SQR</sub>**

When the study cohort has a relatively high degree of relatedness, incorporating the LOCO PGS as an offset alone does not ensure test calibration. More specifically, REGENIE LOCO PGS is typically less predictive of the phenotype than LDAK-KVIK LOCO PGS: using REGENIE LOCO PGS as an offset may fall short of completely eliminating type-I error inflation due to relatedness, while using LDAK-KVIK LOCO PGS may result in deflated test statistics.

SPA<sub>SQR</sub> calibrates the score variance with a **sparse genetic relationship matrix (GRM)** $\Phi \in \mathbb{R}^{n \times n}$, whose entries measure pairwise genetic similarity between subjects:

$$
\Phi_{ij} \;=\; \frac{1}{m} \sum_{k=1}^{m} \frac{(G_{ik} - 2 \mu_k)(G_{jk} - 2 \mu_k)}{2\, \mu_k (1 - \mu_k)},
$$

where $G_{ik}$ is the genotype of variant $k$ in subject $i$, $\mu_k$ is the alternate allele frequency of variant $k$, and $m$ is the number of variants used to construct the GRM. $\Phi_{ij} \approx 0.5$ for parent-offspring or full-sibling pairs, $\approx 0.25$ for half-siblings and grandparent-grandchild pairs, $\approx 0.125$ for first cousins, and so on, decaying toward zero for unrelated pairs. A **sparse** GRM zeroes out all off-diagonal entries below a chosen cutoff (e.g. $0.05$, corresponding roughly to 3rd-degree relatives), substantially reducing the memory footprint of the GRM.

SPA<sub>SQR</sub> plugs the sparse GRM into the score variance:

$$
\widehat{\mathrm{Var}}(S_j) \;=\; \widehat{\sigma}_g^{\,2}(G_j)\, R^{\top}\,\Phi\, R.
$$

For genuinely unrelated cohorts $\Phi \approx I_n$ and the formula reduces to the unrelated variance $\widehat{\sigma}_g^{\,2}(G_j)\, R^{\top}\, R$.

In addition to the files from [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html) — `simu_geno.{bed,bim,fam}`, `simu_geno.pheno`, and the LOCO PGS prediction list `simu_geno_ldak_pred.list` (or `simu_geno_regenie_pred.list`) — Workflow 2 requires two additional files:

```
simu_geno.grm.sp     sparse GRM, three columns: 0-based i, 0-based j, correlation
simu_geno.grm.id     companion ID file: one row per subject, FID  IID
```

The sparse GRM is a **one-time cost per cohort**: once built, the same `.grm.sp` is reused across every trait, every chromosome, and every quantile. The companion `.grm.id` file lists the subject IDs in the same order as the 0-based indices in `.grm.sp`; GRAB auto-detects it from the `.grm.sp` prefix.

A pre-built `simu_geno.grm.sp` / `simu_geno.grm.id` pair, together with all the inputs from Workflow 1, is available in the [`data/`](https://github.com/qhengncsu/SPAsqr.github.io/tree/main/data) folder of this documentation repository for users to replicate this tutorial verbatim.


## Computing the sparse GRM with PLINK 2

The GRM is a concept popularized by the [GCTA](https://yanglab.westlake.edu.cn/software/gcta/) software. However, using GCTA to compute the GRM has historically been relatively time-consuming. Since late 2025, [PLINK 2](https://www.cog-genomics.org/plink/2.0/) also supports sparse GRM construction, which is very simple and efficient to use:

```bash
./plink2 \
    --bfile simu_geno \
    --maf 0.01 \
    --make-grm-sparse 0.05 \
    --threads 8 \
    --out simu_geno
```

- `--maf 0.01` restricts GRM computation to common variants with minor allele frequency $\geq 0.01$.
- `--make-grm-sparse 0.05` retains genetic correlation coefficients above $0.05$ and zeroes out the rest.

The command produces `simu_geno.grm.sp` along with the companion `simu_geno.grm.id`. The `.grm.sp` file is a three-column text file:

```
$ head simu_geno.grm.sp
0   0   1.0024
1   1   0.9981
2   2   1.0107
3   0   0.5012      # half-sibling of subject 0
4   3   0.2503      # first cousin of subject 3
```

Sample $i$, $j$ indices are **0-based** and correspond to the row order in `simu_geno.grm.id`.


## SPA<sub>SQR</sub> association testing with GRM-aware variance

With the LOCO PGS prediction list from [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html) and the sparse GRM in hand, we only need to add the `--sp-grm-plink2` flag to use GRM-aware variance in our association testing procedure:

```bash
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --sp-grm-plink2 simu_geno.grm.sp \
    --out spasqr_results
```


## End-to-end recipes (with INT)

LDAK-KVIK LOCO PGS + sparse GRM, with INT:

```bash
# 1. INT-transform the selected traits
./grab2 --int-pheno --pheno simu_geno.pheno --pheno-name Quantitative1,Quantitative2 --out simu_geno_int

# 2. Train the LOCO PGS on the INT-transformed Y
./ldak6.2.linux \
    --kvik-step1 ldak_step1 \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --mpheno ALL \
    --covar simu_geno.pheno   --covar-names MALE,PC1,PC2,PC3,PC4 \
    --max-threads 8

# 3. Build the LDAK pred-list
cat > simu_geno_ldak_pred.list <<EOF
Quantitative1   $(pwd)/ldak_step1.step1.pheno1.loco.prs
Quantitative2   $(pwd)/ldak_step1.step1.pheno2.loco.prs
EOF

# 4. Build the sparse GRM (one-time cost)
./plink2 \
    --bfile simu_geno \
    --maf 0.01 \
    --make-grm-sparse 0.05 \
    --threads 8 \
    --out simu_geno

# 5. Run SPAsqr with both the LOCO PGS and the sparse GRM
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --sp-grm-plink2 simu_geno.grm.sp \
    --out spasqr_results
```

REGENIE LOCO PGS + sparse GRM, with INT:

```bash
# 1. INT-transform the selected traits
./grab2 --int-pheno --pheno simu_geno.pheno --pheno-name Quantitative1,Quantitative2 --out simu_geno_int

# 2. Train the LOCO PGS on the INT-transformed Y
./regenie \
    --step 1 \
    --bed simu_geno \
    --phenoFile simu_geno_int.txt --phenoColList Quantitative1,Quantitative2 \
    --covarFile simu_geno.pheno   --covarColList MALE,PC1,PC2,PC3,PC4 \
    --bsize 1000 --threads 8 \
    --out simu_geno_regenie

# 3. Build the sparse GRM (one-time cost)
./plink2 \
    --bfile simu_geno \
    --maf 0.01 \
    --make-grm-sparse 0.05 \
    --threads 8 \
    --out simu_geno

# 4. Run SPAsqr with REGENIE's native pred-list and the sparse GRM
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_regenie_pred.list \
    --sp-grm-plink2 simu_geno.grm.sp \
    --out spasqr_results
```

For analyses spanning dozens of phenotypes on a multi-core node, see [All traits at once vs one trait at a time]({{ site.baseurl }}/docs/all-vs-per-trait.html), which wraps LDAK-KVIK Step 1 and SPA<sub>SQR</sub> in per-trait scripts and a bash master launcher.


## When to skip the sparse GRM

The sparse GRM is purely a variance-calibration device: omitting it does not change the score statistic itself, only its reference distribution. It can be omitted when the study cohort has an objectively low degree of relatedness, in which case $R^{\top}\, R \approx R^{\top}\,\Phi\, R$ and the GWAS results will be mostly similar. Additionally, it is well-known that the GRM can be highly inaccurate when computed from genotypes spanning multiple ancestries; we therefore also advise against using this feature for multi-ancestry data.
