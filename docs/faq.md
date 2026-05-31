---
layout: default
title: FAQ
nav_order: 8
description: "Frequently asked questions about SPAsqr."
has_children: false
---

# **Frequently asked questions**

## Methodology

### How is SPA<sub>SQR</sub> different from classical mean-based GWAS (BOLT-LMM / REGENIE)?

Classical methods test $\beta_G$ in
$\mathbb{E}[Y \mid G, X] = X^\top \alpha + G\, \beta_G$. SPA<sub>SQR</sub>
tests $\beta_{G,\tau}$ in
$Q_\tau(Y \mid G, X) = X^\top \alpha_\tau + G\, \beta_{G,\tau}$ at a
grid of $\tau$ values and combines the per-$\tau$ $p$-values via CCT.
On purely mean-shift traits SPA<sub>SQR</sub> matches mean-based GWAS;
on heteroskedastic / dispersion / heavy-tailed traits it wins.

### Why "smoothed" quantile regression?

Conventional QR's check loss is non-differentiable at zero, producing a
discrete rank score $\psi_i = \tau - \mathbb{1}(u_i < 0)$ with
variance $\tau(1-\tau)$. **Smoothing** the check loss with a Gaussian
kernel of bandwidth $h$ replaces the indicator with the Gaussian CDF
$\psi_i^h = \tau - \Phi(-u_i / h)$, which has strictly smaller
variance. That smaller variance directly translates to a smaller
denominator in the score statistic — i.e. higher power — without
biasing the numerator at the population $\beta_G = 0$ null.

### Why CCT for the multi-$\tau$ combination?

Cauchy combination is calibrated under *arbitrary dependence*, which is
exactly the situation across $\tau$ levels for a given variant (the
per-$\tau$ score statistics share the same null residuals and are
highly correlated). Bonferroni would be too conservative; Fisher
requires independence; Stouffer requires known correlation.

### What's the difference between `int` and `irn`?

They are the same algorithm, with two naming conventions in the
literature:

- **IRN** = "inverse rank normal" transform.
- **INT** = "inverse normal transform" (rank-based).

GRAB / SPA<sub>SQR</sub> uses the `int` spelling on the command line.
Older versions accepted `irn`; current versions accept only `int`.

---

## Workflow

### Do I have to run the LOCO-PGS / sparse-GRM steps?

Both are optional. Recommendation:

- **LOCO PGS** — yes, almost always. The power gain is large and the
  cost is one Workflow 1 run per trait.
- **Sparse GRM** — yes if the cohort contains close relatives.
  Skipping when relatives are present causes inflation in the tail.

### Can I pass `pheno.txt` directly, or do I have to pre-INT with `--int-pheno`?

You can pass either. GRAB applies the same INT internally via
`--pheno-transform int` (the default), and INT of an already-INT'd
column is idempotent. The reason to pre-INT with `grab --int-pheno`
is so that LDAK-KVIK / REGENIE Step 1 also see the INT $Y$ — those
external tools don't INT themselves.

### Which `--pheno-transform` should I use?

It has to match the transform used during PGS construction:

| PGS trained on … | `--pheno-transform` |
| ---------------- | ------------------- |
| INT $Y$ (recommended) | `int` (default) |
| Raw $Y$ — LDAK-KVIK / REGENIE internally standardize | `standardize` |
| No PGS (omitting `--pred-list`) | any — `int` is the sensible default |

### What `--spasqr-taus` should I pick?

Start with the default `0.1,0.3,0.5,0.7,0.9` (5 values). If you suspect
heteroskedastic effects, go denser:
`0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9` (9 values). The cost is
sub-linear in `n_taus` (the null-model fit is amortized across all
variants), and CCT remains calibrated.

### What `--spasqr-h-scale` should I pick?

Default `3` in score mode is the sweet spot for $n \gtrsim 10^4$ in
our simulation studies. Don't change unless you have a specific reason
— you'd be picking a different smoothing bias / variance trade-off.

### How many threads should I give it?

Up to the number of **physical** cores. Hyperthreaded logical cores
help marginally; over-subscribing hurts. On a 16-core node, use
`--threads 16`.

---

## Inputs and outputs

### What phenotype-file format does SPA<sub>SQR</sub> expect?

Whitespace-separated, with a header line `FID IID Y1 Y2 ...`. Missing
entries may be `NA`, `na`, `NaN`, `nan`, `.`, `-`, or blank. Rows are
matched to genotype FIDs/IIDs.

### What does SPA<sub>SQR</sub> do with missing phenotype values?

For each trait independently, SPA<sub>SQR</sub> drops the subjects
with missing $Y$ and refits the null model on the non-missing scope.
Multi-trait jobs are not penalized — a subject missing trait 1 still
contributes to trait 2.

### How are multi-allelic variants handled?

Each non-reference allele is split into a separate biallelic test (the
PLINK 2 default). SPA<sub>SQR</sub> does not perform a joint
multi-allelic test.

### Where do per-quantile effect sizes go?

Score mode reports per-$\tau$ $Z$-scores but **not** $\hat\beta_G$
itself — the score statistic is variance-stabilized but not
back-transformed to the $\beta_G$ scale. For effect-size estimates,
re-run on a candidate list with `--spasqr-mode wald`; see
[Effect-size estimation]({{ site.baseurl }}/docs/effect-size-estimation.html).

---

## Errors and edge cases

### `Error: --pheno-transform must be 'raw', 'int', or 'standardize'`

You passed an unrecognized value. Older docs / scripts may use `irn`
(the previous spelling for the inverse normal transform); replace with
`int`.

### `Error: phenotype 'YYY' not found in --pred-list`

The pheno name in the first column of the prediction list must match the
trait column name in your phenotype file. Check spelling and case.

### Very different `Z_tau` signs across $\tau$ at a hit

This is **not an error**. It signals a heteroskedastic / dispersion
effect — the variant changes the spread of $Y$ without shifting the
centre. SPA<sub>SQR</sub>'s CCT combination is designed to catch
exactly these. Inspect $\hat\beta_{G,\tau}$ in Wald mode for a clearer
view.

### Why are my rare-variant $P_\mathrm{CCT}$ different from REGENIE-QRS / NormalQR?

REGENIE-QRS and NormalQR rely on a normal approximation to the score
statistic, which under-reports tail $p$-values for rare variants in
heavy-tailed traits. SPA<sub>SQR</sub> applies saddlepoint
approximation when $|Z| >$ `--spa-z-threshold` (default 2.0), which is
the calibrated answer in the tail. Trust SPA<sub>SQR</sub>.

### Output file is empty

Check the log for `nothing to test (empty marker list)` — your
`--extract` / `--chr` filters likely zeroed out every variant. Also
verify that `--maf` / `--mac` aren't too strict for your cohort.

---

## Performance

### How long does a UK Biobank-scale SPA<sub>SQR</sub> run take?

Approximate, single trait, ~500 K samples, ~20 M imputed variants,
score mode, 9 tau levels, 16 threads, with LOCO PGS + sparse GRM:

- Null-model fit (per chromosome × per $\tau$): a few minutes.
- Per-variant score loop: a few hours.
- Total wall time: 4–8 hours on a single biobank-spec node.

Wald mode on the same scan is roughly 9× slower; you'd restrict it to
a candidate list of a few hundred variants.

### Memory footprint?

Dominated by the genotype chunk in RAM (`--chunk-size 8192` markers).
For $n = 500\text{K}$ samples that is about 4 GB per chunk; with
default settings GRAB peaks at ~6–8 GB resident.

---

## Where to ask for help

- Open an issue at
  [github.com/GeneticAnalysisinBiobanks/GRAB](https://github.com/GeneticAnalysisinBiobanks/GRAB/issues)
  with a minimal repro and the relevant `--log` output.
- Email the corresponding author of the SPA<sub>SQR</sub> manuscript
  (address on the paper).
