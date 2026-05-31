---
layout: default
title: Strategies for improving statistical power
nav_order: 8
description: "Power-improving knobs in the SPAsqr pipeline."
has_children: false
---

# **Strategies for improving statistical power**

SPA<sub>SQR</sub> has several power-improving knobs beyond the
default settings. The table below summarizes what to consider; the
sections that follow give detail.

| Strategy | When it helps | Cost |
| -------- | ------------- | ---- |
| LOCO polygenic score offset | Always — substantial gain on polygenic traits | One Workflow 1 run per trait |
| Sparse GRM for variance calibration | Cohorts with non-trivial relatedness | One sparse-GRM run per cohort |
| INT pre-transform on Y | Heavy-tailed / skewed traits, rare variants | Negligible |
| Wider $\tau$ grid | Heteroskedastic / dispersion effects | Linear in `n_taus`, score mode is fast |
| Smaller bandwidth (`--spasqr-h-scale` ↑) | Sharper tail discrimination at small $n$ | More iterations to converge |
| Multi-trait CCT post-hoc combination | Correlated traits with shared causal variants | Negligible (offline) |

---

## 1. LOCO polygenic score offset

By far the largest power gain in the pipeline. For a polygenic
quantitative trait with heritability $h^2 \approx 0.3$, the
chromosome-specific LOCO PGS soaks up the trans-chromosomal genetic
variance, sharpening the rank-score residual on the tested chromosome
by roughly a factor of $\sqrt{1 - h^2(1-1/22)}$. For trait
architectures with $h^2 > 0.2$ this typically corresponds to a
1.2–1.5× boost in effective sample size.

**Action.** Always run [Workflow 1]({{ site.baseurl }}/docs/workflow-1.html)
and pass `--pred-list` to SPA<sub>SQR</sub> unless your cohort is too
small ($n \lesssim 10^4$) to train a useful PGS.

---

## 2. Sparse GRM for variance calibration

Under unrelated samples, $\Phi = I_n$ and the SPA<sub>SQR</sub>
variance reduces to $\widehat\sigma_g^{\,2}(G)\, R^\top R$. With close
relatives this is anticonservative — score variances are
underestimated and inflation creeps in. The sparse GRM
$\widehat{\mathrm{Var}}(S) = \widehat\sigma_g^{\,2}(G)\, R^\top \Phi R$
restores calibration.

**Action.** Follow [Workflow 2]({{ site.baseurl }}/docs/workflow-2.html)
and pass `--sp-grm-plink2` when:

- King-robust kinship analysis shows non-trivial 1st–3rd degree
  relatives (typical UK Biobank fraction: ~5%).
- The cohort is a founder population or contains clinical pedigrees.

It is **safe** to always supply the GRM — if the cohort is genuinely
unrelated, $\Phi \approx I_n$ and the result is unchanged.

---

## 3. INT pre-transform on Y

For heavy-tailed, skewed, or zero-inflated quantitative traits, the
**inverse normal transform** (Blom plotting position, average-rank
ties) stabilizes the rank-score variance and improves the
saddlepoint-approximation accuracy at rare variants. The transform is
applied **per phenotype, per non-missing scope**:

```bash
./grab --int-pheno --pheno pheno.txt --out pheno_int
```

Then feed `pheno_int.txt` to **both** Workflow 1 and Step 1–2 (and leave
`--pheno-transform` at the default `int`). Because INT of an
already-INT'd column is idempotent, passing the raw `pheno.txt` to
SPA<sub>SQR</sub> with `--pheno-transform int` gives an identical
result — provided the PGS was also trained on INT $Y$.

This is the recommended default. The only cases where you might want
to skip INT are:

- The trait has scientific meaning on its native scale (e.g.
  blood-pressure mmHg) and you want effect sizes you can interpret.
  In that case use `--pheno-transform standardize` and live with
  slightly less robust tail behaviour.
- You are comparing SPA<sub>SQR</sub> against a benchmark method that
  uses raw $Y$; pass `--pheno-transform raw` to both for a fair
  comparison.

---

## 4. Wider $\tau$ grid

The default `--spasqr-taus 0.1,0.3,0.5,0.7,0.9` covers location and
tail behaviour with five quantiles. For heteroskedastic / dispersion
effects, a denser grid catches more signal:

```bash
--spasqr-taus 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9
```

Score mode cost is $\mathcal O(n_\tau)$ per null-model fit, which is
amortized across millions of variants — going from 5 to 9 $\tau$
levels typically adds < 30% to total runtime. The Cauchy combination
test is well-calibrated even when the per-$\tau$ tests are correlated,
so adding more $\tau$ rarely hurts.

For extreme-tail traits (e.g. event-time, dosage of an
event-triggering molecule), include $\tau \in \{0.05, 0.95\}$ in the
grid.

---

## 5. Bandwidth choice

The null-SQR bandwidth $h = \mathrm{IQR}(\tilde Y - \hat Y_{-c}) / k$
trades off bias and variance:

- **Smaller $h$ (larger `--spasqr-h-scale`)** ⇒ less smoothing bias,
  more variance, slower convergence.
- **Larger $h$ (smaller `--spasqr-h-scale`)** ⇒ more smoothing bias,
  smaller score-statistic variance, faster convergence.

The defaults — $k = 3$ in score mode, $k = 10$ in Wald mode — are
based on simulation studies in the SPA<sub>SQR</sub> manuscript and
cover the bias-variance sweet spot for $n \gtrsim 10^4$.

For small studies ($n \lesssim 5000$), try a smaller $k$ (more
smoothing) to reduce solver variance. For very large studies
($n > 10^6$) the choice matters less; stick with the default.

---

## 6. Multi-trait combination

If you analyse multiple correlated traits and expect shared causal
variants, you can run SPA<sub>SQR</sub> separately per trait and then
combine $P_\mathrm{CCT}$ values across traits via the Cauchy
combination test (offline, in R / Python):

```r
P_combined <- function(P_vec) {
  T <- mean(tan((0.5 - P_vec) * pi))
  0.5 - atan(T) / pi
}
```

This is valid under arbitrary dependence among the per-trait
$p$-values and often picks up pleiotropic variants that any single
trait misses. The same idea underlies the per-$\tau$ combination
inside SPA<sub>SQR</sub>.

---

## 7. Solver and tolerance

The default `--spasqr-solver qmme` is a quantile MM-estimator that is
both faster and more numerically stable than the
[Conquer](https://github.com/XiaoouPan/conquer)-style first-order
solver on typical biobank workloads. Switch to `--spasqr-solver
conquer` only as a sanity check.

`--spasqr-tol 1e-7` is the default; tighten to `1e-9` for
analyses where you care about effect-size accuracy at rare variants.
The QMME solver internally clips the tolerance to at most `1e-9`, so
tightening past that point has no effect.

---

> **Note**
> - Strategies 1, 2, 3 are independent — apply all three by default.
> - Strategies 4, 5, 6, 7 are tuning knobs to pull out in
>   problem-specific contexts.
> - When in doubt, run with defaults; SPA<sub>SQR</sub> is calibrated
>   out of the box.
