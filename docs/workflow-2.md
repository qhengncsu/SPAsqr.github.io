---
layout: default
title: "Workflow 2: LOCO PGS + GRM + SPA<sub>SQR</sub>"
nav_order: 4
description: "End-to-end recipe for running SPAsqr with a LOCO polygenic score and a sparse GRM."
has_children: false
---

# **Workflow 2: LOCO PGS + GRM + SPA<sub>SQR</sub>**

When the study cohort contains a non-trivial number of related subjects — close relatives, clinical pedigrees, founder populations — the unrelated-samples variance $\widehat{\sigma}_g^{\,2}(G_j)\,R^{\!\top}R$ used in [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html) is anti-conservative: subjects in the same family share genotype and rank-score, so the score statistic carries more variance than the diagonal form predicts and the resulting $p$-values are inflated, especially in the tails. SPA<sub>SQR</sub> calibrates the score variance with a **sparse genetic relationship matrix (GRM)** $\Phi \in \mathbb{R}^{n \times n}$, whose entries measure pairwise genome-wide genetic similarity between subjects:

$$
\Phi_{ij} \;=\; \frac{1}{m} \sum_{k=1}^{m} \frac{(G_{ik} - 2 p_k)(G_{jk} - 2 p_k)}{2\, p_k (1 - p_k)},
$$

where $G_{ik} \in \{0, 1, 2\}$ is the alternate-allele dosage of variant $k$ in subject $i$, $p_k$ is the alternate allele frequency of variant $k$, and $m$ is the number of variants used to construct the GRM. $\Phi_{ij}$ is approximately $1 + F_i$ on the diagonal (with inbreeding coefficient $F_i$ usually near $0$), $\approx 0.5$ for parent-offspring or full-sibling pairs, $\approx 0.25$ for half-siblings and grandparent-grandchild pairs, $\approx 0.125$ for first cousins, and so on, decaying toward zero for unrelated pairs. The **sparse** GRM zeroes out all off-diagonal entries below a chosen cutoff (e.g. $0.05$, corresponding roughly to 3rd-degree relatives), keeping the matrix linear in the number of related pairs.

GRAB plugs the sparse GRM into the score variance:

$$
\widehat{\mathrm{Var}}(S_j) \;=\; \widehat{\sigma}_g^{\,2}(G_j)\, R^{\!\top}\,\Phi\, R.
$$

For genuinely unrelated cohorts $\Phi \approx I_n$ and the formula reduces to the Workflow 1 expression. This page extends the LOCO-PGS pipeline of [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html) by constructing a sparse GRM and passing it to GRAB alongside the LOCO PGS.

## Data that you will need

In addition to the files from [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html) — `geno.{bed,bim,fam}`, `pheno.txt`, `covar.txt`, and the LOCO PGS prediction list `ldak_pred_list.txt` (or `regenie_step1_pred.list`) — Workflow 2 introduces one additional artifact:

```
geno_sparse.grm.sp    sparse GRM, three columns: 0-based i, 0-based j, kinship
geno_sparse.grm.id    companion ID file: one row per subject, FID  IID
```

The sparse GRM is a **one-time cost per cohort**: once built, the same `.grm.sp` is reused across every trait, every PGS choice, and every `--pheno-transform` setting. The companion `.grm.id` file lists the subject IDs in the same order as the 0-based indices in `.grm.sp`; GRAB auto-detects it from the `.grm.sp` prefix and aligns its rows to the `--pheno` analysis set.

## Computing the sparse GRM with PLINK 2

Since late 2025, [PLINK 2](https://www.cog-genomics.org/plink/2.0/) supports sparse GRM construction natively and is substantially faster than the legacy GCTA path on biobank-scale cohorts. We recommend the PLINK 2 path:

```bash
./plink2 \
    --bfile geno \
    --maf 0.01 \
    --make-grm-sparse 0.05 \
    --threads 8 \
    --out geno_sparse
```

- `--maf 0.01` restricts GRM computation to common variants with minor allele frequency $\geq 0.01$. Rare-variant noise inflates GRM entries, so the MAF filter is essential.
- `--make-grm-sparse 0.05` retains genetic correlation coefficients above $0.05$ and zeroes out the rest. Smaller cutoffs capture more distant cousins at the cost of a denser matrix; $0.05$ corresponds roughly to 3rd-degree relatives and is the value we recommend.

The command produces `geno_sparse.grm.sp` along with the companion `geno_sparse.grm.id`. The `.grm.sp` file is a three-column text file:

```
$ head geno_sparse.grm.sp
0   0   1.0024
1   1   0.9981
2   2   1.0107
3   0   0.5012      # half-sibling of subject 0
4   3   0.2503      # first cousin of subject 3
```

Sample $i$/$j$ indices are **0-based** and correspond to the row order in `geno_sparse.grm.id`. Diagonal entries are the self-kinship $1 + F$ (the inbreeding coefficient $F$ is typically near $0$). Only nonzero off-diagonals above the sparsification cutoff are stored, which keeps the file linear in the number of related pairs rather than quadratic in cohort size.

## Running association testing with GRAB

With the LOCO PGS prediction list from [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html) and the sparse GRM in hand, both null model fitting (SPA<sub>SQR</sub> step 1) and association testing (SPA<sub>SQR</sub> step 2) are performed via a single `grab` call. The only addition compared to Workflow 1 is the `--sp-grm-plink2` flag, which points GRAB at the sparse GRM file (the companion `.grm.id` is auto-detected from the prefix):

```bash
./grab --method SPAsqr \
    --bfile geno \
    --pheno pheno_int.txt \
    --covar covar.txt \
    --pred-list ldak_pred_list.txt \
    --sp-grm-plink2 geno_sparse.grm.sp \
    --out spasqr_results
```

All other required and optional flags — `--method`, `--bfile`, `--pheno`, `--covar`, `--out`, `--pheno-transform`, `--spasqr-taus`, `--spasqr-h-scale`, `--threads`, SNP-level filters, and so on — behave exactly as documented in [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html#running-association-testing-with-grab). The output schema is also unchanged: one tab-delimited result file per trait with per-$\tau$ $p$-values, per-$\tau$ $Z$-scores, and a Cauchy-combined $P_\mathrm{CCT}$. What does change is the value entering each per-$\tau$ $p$-value — the score variance is now the GRM-aware quantity $\widehat\sigma_g^{\,2}(G_j)\,R^{\!\top}\Phi R$ rather than its diagonal counterpart, which restores tail calibration in the presence of relatedness.

### End-to-end recipes (with INT)

LDAK-KVIK LOCO PGS + sparse GRM, with INT:

```bash
# 1. INT-transform the phenotype
./grab --int-pheno --pheno pheno.txt --out pheno_int

# 2. Train the LOCO PGS on the INT-transformed Y
./ldak6.2.linux \
    --kvik-step1 ldak_step1 \
    --bfile geno \
    --pheno pheno_int.txt --mpheno ALL \
    --covar covar.txt \
    --max-threads 8

# 3. Build the LDAK pred-list
./grab --make-ldak-predlist --prefix ldak_step1 --pheno pheno_int.txt --out ldak_pred_list

# 4. Build the sparse GRM (one-time cost)
./plink2 \
    --bfile geno \
    --maf 0.01 \
    --make-grm-sparse 0.05 \
    --threads 8 \
    --out geno_sparse

# 5. Run SPAsqr with both the LOCO PGS and the sparse GRM
./grab --method SPAsqr \
    --bfile geno \
    --pheno pheno_int.txt \
    --covar covar.txt \
    --pred-list ldak_pred_list.txt \
    --sp-grm-plink2 geno_sparse.grm.sp \
    --out spasqr_results
```

REGENIE LOCO PGS + sparse GRM, with INT:

```bash
# 1. INT-transform the phenotype
./grab --int-pheno --pheno pheno.txt --out pheno_int

# 2. Train the LOCO PGS on the INT-transformed Y
./regenie \
    --step 1 \
    --bed geno \
    --phenoFile pheno_int.txt \
    --covarFile covar.txt \
    --bsize 1000 \
    --out regenie_step1

# 3. Build the sparse GRM (one-time cost)
./plink2 \
    --bfile geno \
    --maf 0.01 \
    --make-grm-sparse 0.05 \
    --threads 8 \
    --out geno_sparse

# 4. Run SPAsqr with REGENIE's native pred-list and the sparse GRM
./grab --method SPAsqr \
    --bfile geno \
    --pheno pheno_int.txt \
    --covar covar.txt \
    --pred-list regenie_step1_pred.list \
    --sp-grm-plink2 geno_sparse.grm.sp \
    --out spasqr_results
```

## When the sparse GRM can be skipped

The sparse GRM is purely a variance-calibration device: omitting it does not bias the score statistic itself, only its reference distribution. We recommend skipping the GRM only when the cohort is genuinely unrelated — e.g. independent biobank samples filtered to KING relatedness below the 3rd-degree threshold — in which case $\Phi \approx I_n$ and the GRM-aware variance is numerically indistinguishable from the diagonal one. In any setting that retains first- or second-degree relatives, tail $p$-values for rare variants are systematically inflated without the GRM, and genome-wide significance thresholds are no longer interpretable at face value.
