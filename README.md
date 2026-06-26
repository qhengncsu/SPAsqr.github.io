# Documentation for SPA<sub>SQR</sub>

**SPA<sub>SQR</sub> is implemented in the [GRAB](https://github.com/GeneticAnalysisinBiobanks/GRAB) command-line binary.**

**The full documentation is available at [SPA<sub>SQR</sub> documentation](https://qhengncsu.github.io/SPAsqr.github.io/).**

SPA<sub>SQR</sub> is a saddlepoint-approximated, smoothed quantile
regression (SQR) framework for genome-wide association studies on
quantitative traits. It tests association across multiple quantile
levels simultaneously, accommodates LOCO polygenic scores as offsets
and an optional sparse GRM for variance calibration, and applies
saddlepoint approximation in the tails so that rare variants stay
well-calibrated even under heavy-tailed or skewed phenotypes.

The pipeline is:

1. **Workflow 1 — LOCO PGS + SPAsqr** via LDAK-KVIK or REGENIE.

2. **Workflow 2 — LOCO PGS + GRM + SPAsqr** (optional GRM) via PLINK 2 `--make-grm-sparse`
   (preferred since late 2025) or GCTA.

3. **Steps 1–2 — Run `grab --method SPAsqr`.** Reads phenotype +
   covariates, subtracts the chromosome-specific LOCO PGS as an offset,
   fits the null SQR model per chromosome with bandwidth
   `h = IQR(Y - PGS) / k` (default `k = 3`; `--spasqr-h-scale`), applies SPA in the tails, and writes
   per-phenotype results with per-tau p-values, Z-scores, and a
   Cauchy-combined P_CCT.

4. **(Optional) Effect-size estimation** via `--spasqr-mode wald`.

See the documentation site for installation, command-line reference,
input/output formats, and FAQ.
