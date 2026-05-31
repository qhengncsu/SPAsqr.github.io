---
layout: default
title: "All traits at once vs one trait at a time"
nav_order: 5
description: "Two ways to fan out across many traits: all-trait calls vs per-trait runners with a bash master launcher."
has_children: false
---

# **All traits at once vs one trait at a time**

[Workflow 1]({{ site.baseurl }}/docs/workflow-1.html) and [Workflow 2]({{ site.baseurl }}/docs/workflow-2.html) batch every trait of interest into a single LDAK / REGENIE / GRAB call — `--mpheno ALL`, `--phenoColList Quantitative1,Quantitative2`, `--pheno-name Quantitative1,Quantitative2`. This **all-traits-at-once** style is the simplest possible recipe and is what we first show in the tutorial. It is also the right choice when you are running just a handful of phenotypes.

When the trait list is long (dozens to hundreds of phenotypes), running one trait per process is the more practical pattern, for three reasons:

- **Memory.** An all-trait LDAK or SPA<sub>SQR</sub> process keeps every trait's working buffers in RAM at the same time, and peak memory scales roughly with the number of traits — on a biobank with hundreds of phenotypes this can run a node out of RAM. A per-trait process only holds one trait's working set, so total memory is capped at `MAX_PARALLEL × per-trait-memory`, which you control directly.
- **Debugging and restart.** When one trait in a hundred fails to converge or hits a numerical edge case, an all-trait run forces you to figure out which trait broke and re-run the whole batch. Per-trait runs put every failure in its own log file (`log.spasqr.${trait}`) and let you restart just that one trait.
- **Speed.** Within a single LDAK or SPA<sub>SQR</sub> process, splitting work across traits doesn't give you a clean $N$× speedup. Running $N$ single-threaded processes in parallel typically *beats* one $N$-trait multi-threaded run in speed.

This page gives a tutorial on how to run LDAK-KVIK and SPA<sub>SQR</sub> one trait at a time. The recipe below uses the same `simu_geno` file sets and phenotypes as the previous workflows.

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

`run_spasqr.sh` — assemble a one-line pred-list and run SPA<sub>SQR</sub> for one trait, single-threaded:

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

Each script takes the trait name as its single argument and is deliberately as small as possible — one invocation, one thread, one trait — so that the parallelism lives entirely in the launcher above it.

## Master launcher

The master loops over the trait list and backgrounds each per-trait call, holding the live job count at `MAX_PARALLEL` via `wait -n` — a pure `bash` semaphore that returns as soon as **any** background job finishes, with no dependency on GNU parallel or `xargs -P`:

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

The line `(( $(jobs -r | wc -l) >= MAX_PARALLEL )) && wait -n` is a `bash` semaphore. `jobs -r` lists currently *running* background jobs in this shell, `wc -l` counts them, and the `(( ... ))` test compares that count against `MAX_PARALLEL`. If the cap has been reached, `wait -n` blocks the `for` loop until **any one** background job finishes, freeing a slot for the next iteration to launch a new trait into. The trailing `wait` after the loop then drains the remaining jobs before the phase ends. The net effect is that no more than `MAX_PARALLEL` per-trait processes are alive at any moment — pure bash, no GNU parallel, no `xargs -P`, no external scheduler.

The two phases run sequentially because Phase 2 needs the `.loco.prs` files produced by Phase 1. Within each phase, every trait runs as its own **single-threaded process**, and the **operating system scheduler** spreads those processes across the node's CPU cores. This is the key difference from the all-traits-at-once style: parallelism comes from the OS scheduling $N$ independent processes on different cores, not from multi-threading inside one big LDAK or SPA<sub>SQR</sub> call. On an $n$-core node with $N$ traits and per-trait wall time $T$, the per-phase wall time is approximately $T \cdot \lceil N / n \rceil$ — linear in the trait count, inversely linear in the core count, with no internal-threading overhead to dilute it.

The same per-trait runner scripts drop straight into a SLURM array job (`#SBATCH --array=1-${N}`, with trait $i$ pulled from a manifest file) when you want to scale beyond a single node; the only thing the cluster scheduler replaces is the bash master launcher itself.
