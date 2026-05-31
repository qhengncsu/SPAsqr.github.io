---
layout: default
title: "All traits at once vs one trait at a time"
nav_order: 5
description: "Two ways to fan out across many traits: all-trait calls vs per-trait runners with a bash master launcher."
has_children: false
---

# **All traits at once vs one trait at a time**

[Workflow 1]({{ site.baseurl }}/docs/workflow-1.html) and [Workflow 2]({{ site.baseurl }}/docs/workflow-2.html) batch every trait of interest into a single LDAK / REGENIE / GRAB invocation — `--mpheno ALL`, `--phenoColList Quantitative1,Quantitative2`, `--pheno-name Quantitative1,Quantitative2`. That **all-traits-at-once** style is the simplest possible recipe, fits the tutorial, and is what we recommend when you are running, say, fewer than a handful of phenotypes on a multi-core node.

At biobank scale (UK Biobank, Biobank Japan, etc.) the picture changes: you typically have dozens to a hundred-plus phenotypes to analyze, the per-trait runtime is non-trivial (tens of minutes to a couple of hours), and a single multi-core node has more cores than any one GRAB / LDAK process can saturate. In that regime the natural decomposition is the reverse: **one trait per process**, each process pinned to a single thread, and a thin `bash` master launcher fans the per-trait jobs out across the node's cores in parallel. This page documents that pattern for LDAK-KVIK + SPA<sub>SQR</sub>; REGENIE itself already handles multi-trait Step 1 efficiently in a single invocation, so the per-trait pattern is only useful for the SPA<sub>SQR</sub> stage downstream of REGENIE.

The recipe below uses the same `simu_geno` tutorial fixture as the rest of the documentation, with `Quantitative1` and `Quantitative2` playing the role of a real trait list of length $N$.

## Shared preparatory steps

These two are run once and reused across every trait:

```bash
# INT-transform every trait you plan to analyze (one call, all columns)
./grab2 --int-pheno --pheno simu_geno.pheno --pheno-name Quantitative1,Quantitative2 \
        --out simu_geno_int

# Build the sparse GRM (one-time cost per cohort)
./plink2 --bfile simu_geno --maf 0.01 --make-grm-sparse 0.05 \
         --threads 8 --out simu_geno
```

## Per-trait runners

`run_ldak_step1.sh` — fit the LDAK Step 1 null model for one trait, single-threaded:

```bash
#!/usr/bin/env bash
set -e
trait=$1
./ldak6.2.linux \
    --kvik-step1 ldak_${trait} \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name ${trait} \
    --covar simu_geno.pheno   --covar-names MALE,PC1,PC2,PC3,PC4 \
    --max-threads 1
```

`run_spasqr.sh` — assemble a single-line pred-list and run SPA<sub>SQR</sub> for one trait, single-threaded:

```bash
#!/usr/bin/env bash
set -e
trait=$1
echo "${trait}    $(pwd)/ldak_${trait}.step1.loco.prs" > .pred_${trait}.list

./grab2 --method SPAsqr \
    --bfile simu_geno \
    --pheno simu_geno_int.txt --pheno-name ${trait} \
    --covar simu_geno.pheno   --covar-name MALE,PC1,PC2,PC3,PC4 \
    --pred-list .pred_${trait}.list \
    --sp-grm-plink2 simu_geno.grm.sp \
    --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
    --pheno-transform int \
    --threads 1 \
    --out spasqr_results.${trait}

rm -f .pred_${trait}.list
```

Each script takes the trait name as its single argument and is intentionally as small as possible: no `time -v` wrappers, no embedded logic, just one invocation pinned to one thread.

## Master launcher

The master loops over the trait list and backgrounds each per-trait call, holding the live job count at `MAX_PARALLEL` via `wait -n` — a pure `bash` semaphore that returns as soon as **any** background job finishes, no GNU parallel or `xargs -P` required:

```bash
#!/usr/bin/env bash
TRAITS="Quantitative1 Quantitative2"
MAX_PARALLEL=8                           # set to (#cores) on the node

# Phase 1: LDAK Step 1 for every trait
for trait in $TRAITS; do
    bash run_ldak_step1.sh $trait > log.ldak.${trait} 2>&1 &
    (( $(jobs -r | wc -l) >= MAX_PARALLEL )) && wait -n
done
wait                                      # let all Step 1 jobs finish before SPAsqr

# Phase 2: SPAsqr for every trait
for trait in $TRAITS; do
    bash run_spasqr.sh    $trait > log.spasqr.${trait} 2>&1 &
    (( $(jobs -r | wc -l) >= MAX_PARALLEL )) && wait -n
done
wait
```

The two phases are run sequentially because Phase 2 consumes the `.loco.prs` files produced by Phase 1. Within each phase every trait runs independently, so on an $n$-core node with $N$ traits and per-trait wall time $T$, the all-phase wall time is approximately $T \cdot \lceil N / n \rceil$, regardless of how many traits there are.

## When to prefer which mode

| Setting | Recommended mode |
| ------- | ---------------- |
| Small fixture or a handful of traits, single tutorial node | **All-traits-at-once** ([Workflow 1]({{ site.baseurl }}/docs/workflow-1.html), [Workflow 2]({{ site.baseurl }}/docs/workflow-2.html)) — minimal recipe, no script plumbing. |
| Biobank-scale analysis with many traits on one multi-core node | **One-trait-at-a-time** + bash master launcher, as above — saturates the node and isolates per-trait failures. |
| REGENIE Step 1 alone | Always **all-traits-at-once**; REGENIE Step 1 handles a multi-trait `--phenoColList` efficiently in a single invocation. Only the downstream SPA<sub>SQR</sub> stage benefits from per-trait fan-out. |
| Many nodes (SLURM, PBS, k8s) | The per-trait runner scripts above drop straight into a SLURM array job (`--array=1-${N}`, dispatching trait $i$ from a manifest); no further changes are needed beyond replacing the bash master with the cluster scheduler. |
