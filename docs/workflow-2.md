---
layout: default
title: "Workflow 2: LOCO PGS + GRM + SPA<sub>SQR</sub>"
nav_order: 4
description: "End-to-end recipe for running SPAsqr with a LOCO polygenic score and a sparse GRM."
has_children: false
---

# **Workflow 2: LOCO PGS + GRM + SPA<sub>SQR</sub>**

When the study cohort has a relatively high degree of relatedness, incorporating the LOCO PGS as an offset alone does not ensure test calibration. More specifically, REGENIE LOCO PGS is typically less predictive of the phenotype than LDAK-KVIK LOCO PGS: using REGENIE LOCO PGS as an offset falls short of completely eliminating type-I error inflation due to relatedness, while using LDAK-KVIK LOCO PGS may result in deflated test statistics.

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

In addition to the files from [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html) — `1kg.{bed,bim,fam}`, `1kg.pheno`, and the LOCO PGS prediction list `ldak_pred_list.txt` (or `regenie_step1_pred.list`) — Workflow 2 requires two additional files:

```
1kg.grm.sp     sparse GRM, three columns: 0-based i, 0-based j, correlation
1kg.grm.id     companion ID file: one row per subject, FID  IID
```

The sparse GRM is a **one-time cost per cohort**: once built, the same `.grm.sp` is reused across every trait, every chromosome, and every quantile. The companion `.grm.id` file lists the subject IDs in the same order as the 0-based indices in `.grm.sp`; GRAB auto-detects it from the `.grm.sp` prefix.

A pre-built `1kg.grm.sp` / `1kg.grm.id` pair, together with all the inputs from Workflow 1, is available in the [`examples/`](https://github.com/GeneticAnalysisinBiobanks/GRAB/tree/main/examples) folder of the [GRAB GitHub repository](https://github.com/GeneticAnalysisinBiobanks/GRAB) (also linked at the top-right of this page) for users to replicate this tutorial verbatim.


## Computing the sparse GRM with PLINK 2

The GRM is a concept popularized by the [GCTA](https://yanglab.westlake.edu.cn/software/gcta/) software. However, using GCTA to compute the GRM has historically been relatively time-consuming. Since late 2025, [PLINK 2](https://www.cog-genomics.org/plink/2.0/) also supports sparse GRM construction, which is very simple and efficient to use:

```bash
./plink2 \
    --bfile 1kg \
    --maf 0.01 \
    --make-grm-sparse 0.05 \
    --threads 8 \
    --out 1kg
```

- `--maf 0.01` restricts GRM computation to common variants with minor allele frequency $\geq 0.01$.
- `--make-grm-sparse 0.05` retains genetic correlation coefficients above $0.05$ and zeroes out the rest.

The command produces `1kg.grm.sp` along with the companion `1kg.grm.id`. The `.grm.sp` file is a three-column text file (`i`, `j`, $\Phi_{ij}$):

```
$ head 1kg.grm.sp
0	0	0.65971415
1	0	0.064281474
1	1	0.50382732
2	0	0.064634687
2	1	0.32203796
2	2	0.95093186
3	1	0.13422458
3	2	0.2433665
3	3	0.444169
4	0	0.10572714
```

Sample $i$, $j$ indices are **0-based** and correspond to the row order in `1kg.grm.id`. Diagonal entries are each subject's self-similarity (≈ 1 in expectation); off-diagonal entries above the `--make-grm-sparse 0.05` cutoff are retained, the rest are zeroed.


## SPA<sub>SQR</sub> association testing with GRM-aware variance

With the LOCO PGS prediction list from [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html) and the sparse GRM in hand, we only need to add the `--sp-grm-plink2` flag to use GRM-aware variance in our association testing procedure:

```bash
./grab2 --method SPAsqr \
    --bfile 1kg \
    --pheno 1kg_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar 1kg.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list ldak_pred_list.txt \
    --sp-grm-plink2 1kg.grm.sp \
    --out spasqr_results
```


## End-to-end recipes (with INT)

LDAK-KVIK LOCO PGS + sparse GRM, with INT:

```bash
# 1. INT-transform the selected traits
./grab2 --int-pheno --pheno 1kg.pheno --pheno-name Quantitative1,Quantitative2 --out 1kg_int

# 2. Train the LOCO PGS on the INT-transformed Y
./ldak6.2.linux \
    --kvik-step1 ldak_step1 \
    --bfile 1kg \
    --pheno 1kg_int.txt --mpheno ALL \
    --covar 1kg.pheno   --covar-names MALE,PC1,PC2,PC3,PC4 \
    --max-threads 8

# 3. Build the LDAK pred-list
cat > ldak_pred_list.txt <<EOF
Quantitative1   $(pwd)/ldak_step1.step1.pheno1.loco.prs
Quantitative2   $(pwd)/ldak_step1.step1.pheno2.loco.prs
EOF

# 4. Build the sparse GRM (one-time cost)
./plink2 \
    --bfile 1kg \
    --maf 0.01 \
    --make-grm-sparse 0.05 \
    --threads 8 \
    --out 1kg

# 5. Run SPAsqr with both the LOCO PGS and the sparse GRM
./grab2 --method SPAsqr \
    --bfile 1kg \
    --pheno 1kg_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar 1kg.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list ldak_pred_list.txt \
    --sp-grm-plink2 1kg.grm.sp \
    --out spasqr_results
```

REGENIE LOCO PGS + sparse GRM, with INT:

```bash
# 1. INT-transform the selected traits
./grab2 --int-pheno --pheno 1kg.pheno --pheno-name Quantitative1,Quantitative2 --out 1kg_int

# 2. Train the LOCO PGS on the INT-transformed Y
./regenie \
    --step 1 \
    --bed 1kg \
    --phenoFile 1kg_int.txt --phenoColList Quantitative1,Quantitative2 \
    --covarFile 1kg.pheno   --covarColList MALE,PC1,PC2,PC3,PC4 \
    --bsize 1000 --threads 8 \
    --out regenie_step1

# 3. Build the sparse GRM (one-time cost)
./plink2 \
    --bfile 1kg \
    --maf 0.01 \
    --make-grm-sparse 0.05 \
    --threads 8 \
    --out 1kg

# 4. Run SPAsqr with REGENIE's native pred-list and the sparse GRM
./grab2 --method SPAsqr \
    --bfile 1kg \
    --pheno 1kg_int.txt --pheno-name Quantitative1,Quantitative2 \
    --covar 1kg.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list regenie_step1_pred.list \
    --sp-grm-plink2 1kg.grm.sp \
    --out spasqr_results
```


## When to skip the sparse GRM

The sparse GRM is purely a variance-calibration device: omitting it does not bias the score statistic itself, only its reference distribution. It can be omitted when the study cohort has an objectively low degree of relatedness, in which case $R^{\top}\, R \approx R^{\top}\,\Phi\, R$ and the GWAS results are mostly similar. Additionally, it is well-known that the GRM can be highly inaccurate when computed from genotypes spanning multiple ancestries; we therefore also advise against using this feature for multi-ancestry data.
