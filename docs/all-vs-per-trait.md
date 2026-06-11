---
layout: default
title: "All-traits-at-once vs one-trait-at-a-time"
nav_order: 5
description: "The fastest recipe for many phenotypes: build LOCO PGS one trait at a time in parallel, then run SPAsqr over all traits at once."
has_children: false
---

# **All-traits-at-once vs one-trait-at-a-time**

[Workflow 1]({{ site.baseurl }}/docs/workflow-1.html) and [Workflow 2]({{ site.baseurl }}/docs/workflow-2.html) analyze multiple traits in a single LDAK-KVIK / REGENIE / GRAB call. This **all-traits-at-once** style is the simplest possible recipe to use our software, and is good enough when we are analyzing just a handful of phenotypes.

When analyzing dozens or even hundreds of phenotypes, the fastest recipe is a **hybrid** that treats the two stages differently:

- **LOCO PGS construction** (LDAK-KVIK Step 1) is fastest run **one trait at a time**, as many single-threaded processes in parallel.
- **SPA<sub>SQR</sub> association testing** is fastest run **all traits at once**, in a single multi-threaded call.

So we build the LOCO PGS per trait in parallel, assemble one combined prediction list, and then run SPA<sub>SQR</sub> once over all traits.

## LDAK-KVIK Step 1 performs better one-trait-at-a-time

LDAK-KVIK Step 1 is more efficient run one trait at a time: $N$ single-threaded processes in parallel beat one $N$-threaded multi-trait call. SPA<sub>SQR</sub> is the opposite — a single all-traits call streams the genotypes only once and shares that work across every trait, so it beats launching a separate SPA<sub>SQR</sub> process per trait.

## Prepare the shared inputs

A single one-time step produces the file the launcher assumes — the INT-transformed phenotype `simu_geno_int.txt`:

```bash
# INT-transform every trait once (all columns in a single call)
./grab2 --int-pheno --pheno simu_geno.pheno --pheno-name Quantitative1,Quantitative2 \
        --out simu_geno_int
```

For related cohorts you can additionally supply a sparse GRM via `--sp-grm-plink2`; see [Workflow 2]({{ site.baseurl }}/docs/workflow-2.html).

## Build the LOCO PGS for every trait in parallel

This snippet reads the trait list from the phenotype-file header, launches one single-threaded LDAK-KVIK Step 1 per trait **all at once** in the background, waits for them to finish, then assembles a single prediction list:

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

Every trait runs as its own single-threaded process, and the operating-system scheduler spreads those processes across the node's cores; each failure lands in its own `log.ldak.${trait}`. The trailing `wait` ensures every LOCO PGS is on disk before the prediction list is assembled.

## Run SPA<sub>SQR</sub> over all traits at once

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

Omitting `--pheno-name` makes SPA<sub>SQR</sub> analyze every trait column in `--pheno`; give it as many `--threads` as the node has cores.

## When to prefer which mode

| Setting | Recommended recipe |
| ------- | ---------------- |
| A handful of traits | **All-traits-at-once** for both stages ([Workflow 1]({{ site.baseurl }}/docs/workflow-1.html), [Workflow 2]({{ site.baseurl }}/docs/workflow-2.html)) — simplest, no script plumbing. |
| Dozens to hundreds of traits | **Per-trait LOCO PGS in parallel + one all-traits SPA<sub>SQR</sub> call**, as on this page. |
