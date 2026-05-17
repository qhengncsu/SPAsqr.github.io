---
layout: default
title: Effect-size estimation
nav_order: 6
description: "Per-marker per-tau β̂_G and SE via SPAsqr Wald mode."
has_children: false
---

# **Effect-size estimation — `--spasqr-mode wald`**

The default `--spasqr-mode score` is optimized for **genome-wide
screening**: it tests $H_0\!: \beta_G = 0$ against
$H_1\!: \beta_G \neq 0$ at every $\tau$ using a rank-score statistic
computed from a single null-model fit per chromosome. Score mode
returns calibrated $p$-values and signed $Z$-scores per marker, but
**not** the effect-size estimate $\hat\beta_G$ itself.

For follow-up analyses — characterizing the magnitude of a top hit,
plotting per-quantile effect curves, or feeding effects into a
downstream meta-analysis — switch to **Wald mode** with
`--spasqr-mode wald`. Wald mode fits the *full* model
$Q_\tau(Y \mid X, G_j) = X^\top \beta_\tau + G_j\, \beta_{G,\tau}$
once **per (marker, $\tau$)** by smoothed M-estimation, and returns

- $\hat\beta_{G,\tau}$ — the per-quantile genetic effect on the
  transformed $Y$ scale,
- $\widehat{\mathrm{SE}}(\hat\beta_{G,\tau})$ — sandwich-form standard
  error using the bandwidth-aware Hessian and the empirical score
  meat,
- $Z_{j,\tau} = \hat\beta_{G,\tau} / \widehat{\mathrm{SE}}$ — Wald
  statistic,
- $P_{j,\tau}$ — two-sided $p$-value from the standard normal.

## When to use Wald mode

| Use case | Recommended mode |
| -------- | ---------------- |
| Genome-wide screening (millions of variants) | **`score`** — much faster, gives calibrated $P_\mathrm{CCT}$. |
| Effect-size estimation on top hits (≤ a few hundred variants) | **`wald`** — gives $\hat\beta_G$ + SE per $\tau$. |
| Generating summary statistics for meta-analysis | **`wald`** restricted via `--extract`. |
| Sensitivity analyses (alternative bandwidth, tau-specific) | **`wald`** — finer control over per-marker fit. |

> Wald mode refits the full model per (marker, $\tau$), so total cost
> scales as $\mathcal O(\text{n\_markers} \times \text{n\_taus})$
> versus score mode's $\mathcal O(\text{n\_markers})$. For a
> genome-wide scan with 9 tau levels this is roughly 9× slower than
> score mode — fine for a candidate list, infeasible for the full
> genome.

## Example

Restrict to a set of GWAS hits (one ID per line in `hits.txt`) and
re-run in Wald mode:

```bash
./grab --method SPAsqr \
     --spasqr-mode wald \
     --bfile geno \
     --pheno pheno_int.txt --pheno-name Y1 \
     --covar covar.txt --covar-name covar1,covar2 \
     --pred-list grab_predlist.txt \
     --spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9 \
     --extract hits.txt \
     --spasqr-h-scale 10 \
     --threads 8 \
     --out spasqr_wald
```

This produces `spasqr_wald.Y1.SPAsqr` with one row per
(marker, $\tau$):

```
CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  TAU  BETA  SE  Z  P
```

Column meanings:

- `CHROM POS ID REF ALT MISS_RATE ALT_FREQ MAC HWE_P` — same as score
  mode.
- `TAU` — the quantile level for this row.
- `BETA` — Wald estimate $\hat\beta_{G,\tau}$ on the
  pheno-transform scale (e.g. INT scale if
  `--pheno-transform int`).
- `SE` — sandwich-variance standard error.
- `Z` — Wald statistic $\hat\beta_{G,\tau} / \widehat{\mathrm{SE}}$.
- `P` — two-sided $p$-value, $P = 2 \cdot \Phi(-|Z|)$.

## Interpreting per-quantile effects

For a single hit, you typically have 5–9 rows (one per $\tau$). The
**effect-quantile curve** $\tau \mapsto \hat\beta_{G,\tau}$ tells you
where in the conditional distribution of $Y$ the variant acts:

- **Flat curve** — mean-only effect; classical OLS / linear GWAS
  would have caught it. SPA<sub>SQR</sub> still gives the same hit
  with slightly conservative SE.
- **Monotone in $\tau$** — location-scale effect; the variant shifts
  both the centre and the tails.
- **U-shaped or non-monotone** — heteroskedastic / dispersion effect;
  the variant changes the *spread* of $Y$ without shifting the
  centre. These are the hits where SPA<sub>SQR</sub> typically beats
  classical mean-based GWAS.

Plot $\hat\beta_{G,\tau}$ against $\tau$ with $\pm 1.96 \cdot$ SE
error bars; cross-check the sign pattern against the per-$\tau$
$Z$-scores from score mode.

## Bandwidth and solver notes

- Wald mode defaults to `--spasqr-h-scale 10` (narrower bandwidth than
  score mode's 3). This reduces smoothing bias in $\hat\beta_G$ at the
  cost of slower convergence; it is the right default when you care
  about the *estimate*, not just the test.
- The same `--spasqr-solver` choice applies (`qmme` default,
  `conquer` as alternative). `qmme` is the recommended solver for both
  modes.
- `--spasqr-tol` controls convergence tolerance for the M-estimation
  iterations. The default `1e-7` is usually fine; tighten to `1e-9`
  for very rare variants or extreme tail $\tau$.

## When NOT to use Wald mode

- For **genome-wide screening**, always use score mode first. Score
  mode's $P_\mathrm{CCT}$ is the calibrated GWAS $p$-value; Wald
  per-$\tau$ Z-scores are not multiplicity-adjusted.
- For **rare-variant testing** at extreme $\tau$ (e.g. $\tau \in
  \{0.05, 0.95\}$) with sparse data, Wald sandwich SEs can be
  unstable. Use score mode instead, which leverages SPA in the tails.

> **Note**
> - Wald mode honors the same `--pred-list`, `--sp-grm-plink2`,
>   `--pheno-transform` and `--spasqr-taus` flags as score mode.
> - The output schema differs: Wald is **one row per (marker, $\tau$)**
>   with `TAU BETA SE Z P` columns; score is **one row per marker**
>   with `P_CCT P_tau{val}... Z_tau{val}...` columns.
