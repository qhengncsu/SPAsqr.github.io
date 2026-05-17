---
layout: default
title: Running SPAsqr
nav_order: 5
description: "Run the SPAsqr association test via the GRAB CLI."
has_children: false
---

# **Running SPA<sub>SQR</sub>**

A single `./grab --method SPAsqr` call:

1. Reads the phenotype + covariates and applies the chosen
   `--pheno-transform`.
2. Loops over chromosomes. For each chromosome $c$:
   - Subtracts the LOCO PGS column $\hat Y_{-c}$ as an offset.
   - Fits the null **smoothed quantile regression** model at every
     requested $\tau$ with bandwidth
     $h = \mathrm{IQR}(\tilde Y - \hat Y_{-c}) / k$.
   - For every variant on chromosome $c$, computes the rank-score
     statistic $S_j = G_j^{\!\top} R_\tau$, the GRM-aware variance
     $\widehat\sigma_g^{\,2}(G_j)\, R_\tau^{\!\top}\Phi R_\tau$, and a
     $Z$-score $Z_{j,\tau} = S_j / \widehat{\mathrm{sd}}(S_j)$.
   - Applies saddlepoint approximation when $|Z_{j,\tau}|$ exceeds
     `--spa-z-threshold` (default 2.0).
3. Combines per-$\tau$ $p$-values into one $P_\mathrm{CCT}$ via the
   Cauchy combination test.
4. Writes one tab-delimited result file per phenotype.

## Worked example

We assume the artifacts produced by the previous two steps:

- `geno.{bed,bim,fam}` — PLINK fileset of variants to test.
- `pheno.txt` (or `pheno_int.txt`) — phenotype file, two traits
  `Y1, Y2`.
- `covar.txt` — covariate file with columns `covar1, covar2`.
- `geno_sparse.grm.sp` — sparse GRM (optional but recommended).
- `grab_predlist.txt` — prediction list pairing each phenotype with its
  LOCO PGS file.

The prediction list looks like:

```
Y1    /path/to/ldak_step1.step1.Y1.loco.prs
Y2    /path/to/ldak_step1.step1.Y2.loco.prs
```

The full command:

```bash
./grab --method SPAsqr \
     --bfile geno \
     --pheno pheno_int.txt --pheno-name Y1,Y2 \
     --covar covar.txt --covar-name covar1,covar2 \
     --sp-grm-plink2 geno_sparse.grm.sp \
     --pred-list grab_predlist.txt \
     --spasqr-h-scale 3 \
     --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
     --pheno-transform int \
     --threads 8 \
     --out spasqr_results
```

This produces

```
spasqr_results.Y1.SPAsqr
spasqr_results.Y2.SPAsqr
```

each with the schema described in [Output format](#output-format).

---

## Required flags

| Flag | Value | Purpose |
| ---- | ----- | ------- |
| `--method SPAsqr` | — | Selects the SPA<sub>SQR</sub> method. |
| `--bfile` / `--pfile` / `--vcf` / `--bgen` | PREFIX or FILE | Genotype input. Exactly **one** of these is required. PLINK 1 (`bed/bim/fam`), PLINK 2 (`pgen/pvar/psam`), VCF/BCF, and BGEN are all supported. |
| `--pheno` | FILE | Phenotype file (`FID IID Y1 Y2 …` header). |
| `--out` | PREFIX | Output prefix — GRAB appends `.<phenoname>.SPAsqr` per trait. |

## Frequently used optional flags

| Flag | Default | Purpose |
| ---- | ------- | ------- |
| `--pheno-name` | (all) | Comma-separated list of trait columns to test. Omitting tests every Y column. |
| `--covar` + `--covar-name` | — | Covariate file + columns. An intercept is added automatically. |
| `--pred-list` | — | LOCO PGS file pairing. See [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html). |
| `--sp-grm-plink2` | — | Sparse GRM (`.grm.sp`). See [Step 0b]({{ site.baseurl }}/docs/step-0b-sparse-grm.html). |
| `--spasqr-taus` | `0.1,0.3,0.5,0.7,0.9` | Quantile grid. Max 20 values. |
| `--spasqr-h-scale` | `3` (score) / `10` (wald) | Bandwidth divisor: $h = \mathrm{IQR}(Y)/k$. Smaller $k$ ⇒ larger $h$ ⇒ smoother. |
| `--spasqr-h` | (auto) | Explicit bandwidth. Mutually exclusive with `--spasqr-h-scale`. |
| `--spasqr-solver` | `qmme` | Null-model solver: `qmme` (faster, default) or `conquer`. |
| `--spasqr-mode` | `score` | `score` for genome-wide screening; `wald` for [effect-size estimation]({{ site.baseurl }}/docs/effect-size-estimation.html). |
| `--pheno-transform` | `int` | `raw` / `int` / `standardize`. See note below. |
| `--threads` | `1` | Worker threads. Scales linearly up to # physical cores. |
| `--maf` / `--mac` / `--geno` / `--hwe` | `1e-5` / `10` / `0.1` / `0` | Standard variant QC filters. |
| `--chr` | (all) | Restrict to a subset of chromosomes (e.g. `1-4,22`). |
| `--extract` / `--exclude` | — | SNP include / exclude lists. |
| `--keep` / `--remove` | — | Subject include / exclude lists. |

Run `grab --help SPAsqr` for the full reference.

---

## `--pheno-transform` and consistency with the PGS scale {#pheno-transform}

`--pheno-transform` controls the transform GRAB applies **internally**
to the non-missing trait values before subtracting the LOCO PGS and
fitting the null SQR model:

- `raw` — no transform; SPA<sub>SQR</sub> fits $Y$ as supplied.
- `int` (default) — inverse normal transform (Blom plotting position,
  average-rank ties).
- `standardize` — centre and scale to unit variance.

Because the transform is applied internally, **you may pass either the
raw `pheno.txt` or the pre-INT'd `pheno_int.txt`** as the
`--pheno` input — INT of an already-INT'd column reproduces the same
column (ranks are preserved), so the two are equivalent under
`--pheno-transform int`.

**The transform must match the one used during PGS construction**, so
that the offset and the response live on the same scale:

| If the LOCO PGS was trained on … | … pass to SPA<sub>SQR</sub> |
| -------------------------------- | ---------------- |
| INT-transformed $Y$ (`pheno_int.txt`) | `--pheno-transform int` |
| Raw $Y$ (`pheno.txt`)                 | `--pheno-transform standardize` — LDAK-KVIK / REGENIE internally standardize. |
| Pre-standardized $Y$                  | `--pheno-transform standardize` |

We recommend the **INT path** documented in
[Workflow 1]({{ site.baseurl }}/docs/workflow-1.html): pre-INT with
`grab --int-pheno`, feed `pheno_int.txt` to LDAK-KVIK, and leave
`--pheno-transform` at its default.

---

## Bandwidth: `--spasqr-h-scale` vs `--spasqr-h`

The null-model SQR fit smooths the check loss with a Gaussian kernel
of bandwidth $h$. SPA<sub>SQR</sub> picks $h$ automatically as
$h = \mathrm{IQR}(\tilde Y - \hat Y_{-c}) / k$ where $\tilde Y$ is the
transformed response and $k$ is `--spasqr-h-scale`:

- $k = 3$ (default in **score** mode) — matches our main analyses,
  approximately equal to Silverman's rule for unit-variance Gaussian.
- $k = 10$ (default in **wald** mode) — narrower bandwidth, less
  smoothing bias, slower convergence.
- Pass `--spasqr-h FLOAT` to override the auto rule with an explicit
  bandwidth (advanced — useful for sensitivity analyses).

`--spasqr-h` and `--spasqr-h-scale` are mutually exclusive.

---

## Genotype QC

The defaults match common biobank-GWAS conventions and are usually
fine:

- `--maf 1e-5` — drop singletons.
- `--mac 10` — drop variants with fewer than 10 minor alleles.
- `--geno 0.1` — drop variants with >10% missingness.
- `--hwe 0` — HWE filter disabled (set to e.g. `1e-10` to enable).

For a quick test run on a single chromosome, add `--chr 22`.

---

## Output format {#output-format}

GRAB writes one tab-delimited file per phenotype:

```
spasqr_results.Y1.SPAsqr
spasqr_results.Y2.SPAsqr
```

**Score mode** (default): one row per variant, with per-$\tau$
$p$-values, Z-scores, and a CCT-combined $p$-value.

```
CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P
P_CCT  P_tau0.1 ... P_tau0.9  Z_tau0.1 ... Z_tau0.9
```

Field meanings:

- `CHROM POS ID REF ALT` — variant identifiers from the input genotype
  file.
- `MISS_RATE` — per-marker missing rate (post-`--keep` / `--remove`).
- `ALT_FREQ` — alternative-allele frequency on the analysis subset.
- `MAC` — minor allele count.
- `HWE_P` — Hardy-Weinberg equilibrium $p$-value.
- `P_CCT` — Cauchy-combined $p$-value across the $\tau$ grid. **This
  is the headline GWAS $p$-value for the variant.**
- `P_tau{val}` — per-quantile $p$-value (one column per
  `--spasqr-taus` entry).
- `Z_tau{val}` — per-quantile signed $Z$-score
  ($Z = S / \widehat{\mathrm{sd}}(S)$, post-SPA).

Suggested Manhattan / QQ inputs: plot $-\log_{10}(\texttt{P\_CCT})$
against `POS`. For tau-specific follow-up at a hit, inspect
`Z_tau0.1, …, Z_tau0.9` — a monotone trend across $\tau$ suggests a
mean-only effect, while non-monotone or U-shaped patterns suggest
heteroskedastic or dispersion effects.

**Wald mode** (`--spasqr-mode wald`): one row per (variant, $\tau$),
with $\hat\beta_G$ and SE instead of $p$-values only. See
[Effect-size estimation]({{ site.baseurl }}/docs/effect-size-estimation.html).

---

## Threading and resource use

- **Threads.** SPA<sub>SQR</sub> uses a chunk-level work-stealing pool;
  scaling is roughly linear up to the number of physical cores. The
  null-model fit (per chromosome × per $\tau$) is parallelized
  internally; the per-variant scoring loop is parallel across chunks.
- **Memory.** Dominated by the genotype matrix chunk in RAM
  (`--chunk-size`, default 8192 markers). At $n = 500\text{K}$ that
  is about 4 GB per chunk; reduce `--chunk-size` if RAM-constrained.
- **Disk.** Output is plain text; gz/zst compression can be enabled
  via `--compression gz` / `--compression zst`.

> **Note**
> - Both `--sp-grm-plink2` and `--pred-list` are optional. You can run
>   SPA<sub>SQR</sub> with either, both, or neither.
> - Omitting the GRM falls back to $\Phi = I_n$
>   ($\widehat{\mathrm{Var}}(S_j) =
>   \widehat\sigma_g^{\,2}(G_j)\, R^{\!\top}R$).
> - Omitting the PGS skips the polygenic offset and forfeits the
>   associated power gain.
> - We recommend **always including the LOCO PGS**, and additionally
>   supplying a sparse GRM when the cohort has non-trivial relatedness.
