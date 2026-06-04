---
layout: default
title: "All-traits-at-once vs one-trait-at-a-time"
nav_order: 5
description: "Two ways to run many traits: a single all-trait call vs one process per trait with a bash master launcher."
has_children: false
---

# **All-traits-at-once vs one-trait-at-a-time**

[Workflow 1]({{ site.baseurl }}/docs/workflow-1.html) and [Workflow 2]({{ site.baseurl }}/docs/workflow-2.html) analyze multiple traits in a single LDAK-KVIK / REGENIE / GRAB call. This **all-traits-at-once** style is the simplest possible recipe to use our software. This style is good enough when we are analyzing just a handful of phenotypes.

When it is desired to analyze dozens or even hundreds of phenotypes, running one trait per process/command is the better practise, for the following three reasons:

- **Memory.** A multi-trait LDAK or SPA<sub>SQR</sub> process keeps every trait's working information in RAM at the same time, and peak memory scales roughly with the number of traits — at biobank scale this can easily run a node out of RAM. If we only analyze one trait at each call, total memory is more manageable.
- **Debugging and restart.** When one trait in a hundred fails to complete the analysis due to an error, a multi-trait run forces you to figure out which trait broke and re-run the whole thing. Per-trait runs put every failure in its own log file (`log.spasqr.${trait}`) and let you debug and restart just that one trait.
- **Speed.** Within a single LDAK-KVIK or SPA<sub>SQR</sub> process, using $N$ threads via `---max-threads N` or `---threads N` doesn't give you a clean $N$× speedup. Running $N$ single-threaded processes in parallel typically *beats* analyzing $N$ traits using $N$ threads in computation efficiency.

This page gives a tutorial on how to run LDAK-KVIK and SPA<sub>SQR</sub> one-trait-at-a-time. The recipe below uses the same `simu_geno` file sets and phenotypes as the previous workflows.

## Per-trait runner scripts

`run_ldak_step1.sh` runs LDAK-KVIK Step 1 for one trait, single-threaded:

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

`run_spasqr.sh` creates an one-line prediction list and run SPA<sub>SQR</sub> for one trait, single-threaded:

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

Each script takes the trait name as its single argument and is deliberately as small as possible — one call, one thread, one trait — so that any parallelism lives entirely outside of it. We may use the following bash script to run several such one-thread processes concurrently. 

## Master launcher

The master derives the trait list automatically from the phenotype-file header — `head -n1` reads the header, `tr -s ' \t' ' '` collapses runs of spaces/tabs to a single space, and `cut -d' ' -f3-` keeps every column after the `FID IID` keys (adjust `-f3-` to `-f2-` if your file has a single `IID` key column). It then loops over that list and sends each per-trait call to run in the background via `&`, holding the live job count at `MAX_PARALLEL` via `wait -n` — a semaphore that returns as soon as **any** background job finishes:

```bash
#!/usr/bin/env bash
PHENO=simu_geno_int.txt
MAX_PARALLEL=8                           # number of concurrent jobs

# Trait list = every column after the FID/IID keys in the phenotype header
TRAITS=$(head -n1 "$PHENO" | tr -s ' \t' ' ' | cut -d' ' -f3-)

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

The line `(( $(jobs -r | wc -l) >= MAX_PARALLEL )) && wait -n` is a `bash` semaphore. `jobs -r` lists currently *running* background jobs in this shell, `wc -l` counts them, and the `(( ... ))` test compares that count against `MAX_PARALLEL`. If the cap has been reached, `wait -n` blocks the `for` loop until **any one** background job finishes, freeing a slot for the next iteration to launch a new single-trait call.  The net effect is that no more than `MAX_PARALLEL` per-trait processes are alive at any moment. The trailing `wait` ensures that all step 1 jobs are finished before we start step 2 association testing.

The two phases run sequentially because Phase 2 needs the `.loco.prs` files produced by Phase 1. Within each phase, every trait runs as its own **single-threaded process**, and the **operating system scheduler** spreads those processes across the node's CPU cores. This is the key difference from the all-traits-at-once style: parallelism comes from the OS scheduling $N$ independent processes on different cores, not from multi-threading inside one big LDAK-KVIK or SPA<sub>SQR</sub> call. 
## When to prefer which mode

| Setting | Recommended mode |
| ------- | ---------------- |
| A handful of traits | **All-traits-at-once** ([Workflow 1]({{ site.baseurl }}/docs/workflow-1.html), [Workflow 2]({{ site.baseurl }}/docs/workflow-2.html)) — simple and easy to use, no bash script plumbing. |
| Dozens to hundreds of traits | **One-trait-at-a-time** + bash master launcher, as on this page, typically runs faster than **All-traits-at-once** and isolates per-trait failures. |
