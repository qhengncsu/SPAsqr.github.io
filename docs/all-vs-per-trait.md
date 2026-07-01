---
layout: default
title: "LDAK-KVIK Step 1 runs faster one-trait-at-a-time"
nav_order: 5
description: "For many phenotypes: build LOCO PGS one trait at a time in parallel, then run SPAsqr over all traits at once."
has_children: false
---

# **LDAK-KVIK Step 1 runs faster one-trait-at-a-time**

[Workflow 1]({{ site.baseurl }}/docs/workflow-1.html) and [Workflow 2]({{ site.baseurl }}/docs/workflow-2.html) analyze all traits in a single command. However, we have found that **LDAK-KVIK Step 1** is faster run **one trait per process**: $N$ single-threaded processes in parallel beat one $N$-threaded call over all $N$ traits. **SPA<sub>SQR</sub>** is the opposite — one multi-trait call is usually much faster, mainly from **shared I/O**: genotypes are read once and reused across traits.

Here we illustrate computing the LDAK-KVIK LOCO PGS using a separate process for each trait, then assemble them into one prediction list to supply to GRAB.

## Prepare the shared inputs

Again, we first generate the INT-transformed phenotype `simu_geno_int.txt`:

```bash
# INT-transform every trait once (all columns in a single call)
./grab2 --int-pheno --pheno simu_geno.pheno --pheno-name Quantitative1,Quantitative2 \
        --out simu_geno_int
```

## Build the LOCO PGS for every trait in parallel

This snippet reads the trait list from the phenotype header and launches a single-threaded LDAK-KVIK Step 1 for each trait in the background, then assembles the prediction list once all finish.

```bash
#!/usr/bin/env bash
set -e
PHENO=simu_geno_int.txt

# Trait list = every column after the FID/IID keys in the phenotype header
TRAITS=$(head -n1 "$PHENO" | tr -s ' \t' ' ' | cut -d' ' -f3-)

# One single-threaded LDAK-KVIK Step 1 per trait, all launched in the background
for trait in $TRAITS; do
    ./ldak6.2.linux --kvik-step1 ldak_${trait} --bfile simu_geno \
        --pheno "$PHENO" --pheno-name ${trait} \
        --covar simu_geno.pheno --covar-names MALE,PC1,PC2,PC3,PC4 \
        --max-threads 1 > log.ldak.${trait} 2>&1 &
done
wait                                      # all LOCO PGS finished

# Build one prediction list covering every trait
: > simu_geno_ldak_pred.list
for trait in $TRAITS; do
    echo "${trait}    $(pwd)/ldak_${trait}.step1.loco.prs" >> simu_geno_ldak_pred.list
done
```

Every trait runs as its own single-threaded process, and the operating-system scheduler spreads those processes across the machine's CPU cores. The trailing `wait` ensures every LOCO PGS is on disk before the prediction list is assembled.

## Run SPA<sub>SQR</sub> for all traits

```bash
./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno_int.txt \
    --covar simu_geno.pheno --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list simu_geno_ldak_pred.list \
    --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
    --pheno-transform int \
    --threads 8 \
    --out spasqr_results
```

Omitting `--pheno-name` makes SPA<sub>SQR</sub> analyze every trait column in `--pheno`; `--threads` should not exceed the number of available cores on the machine. Thread contention can slow down the program significantly.
