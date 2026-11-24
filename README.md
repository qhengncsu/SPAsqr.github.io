# Documentation for SPA<sub>SQR</sub>

SPA<sub>SQR</sub> is a method for genome-wide association studies of quantitative traits. It uses smoothed quantile regression to test for association at several quantiles and combines the results into a single test for each variant. A saddlepoint approximation in the tails keeps rare-variant tests well calibrated even when the trait distribution is skewed or heavy-tailed.

SPA<sub>SQR</sub> is available in the [GRAB](https://github.com/GeneticAnalysisinBiobanks/GRAB) command-line program. The [documentation](https://qhengncsu.github.io/SPAsqr.github.io/) covers installation, commands, file formats, and frequently asked questions.

1. **Prepare the polygenic scores.** Both workflows use leave-one-chromosome-out (LOCO) polygenic scores from LDAK-KVIK or REGENIE. Workflow 1 uses these scores alone. Workflow 2 also uses a sparse genetic relationship matrix (GRM) to calibrate the variance of the score statistic. To build the GRM, we recommend PLINK 2 with `--make-grm-sparse`, available since late 2025. GCTA is another option.
2. **Test for association.** Run `grab --method SPAsqr` with your phenotype and covariate data. For each chromosome, GRAB subtracts the LOCO score from the phenotype and fits the null smoothed quantile regression model with bandwidth `h = IQR(Y - PGS) / k`. The default is `k = 3`, which you can change with `--spasqr-h-scale`. Each trait has its own output file with quantile-specific −log₁₀ p-values (`LOG10P_tau`), Z-scores (`Z_Norm_tau`), and the Cauchy-combined result (`LOG10P_CCT`).
3. **Estimate effects at selected variants (optional).** Use `--spasqr-mode wald` to estimate effects at each quantile and test whether they differ across quantiles.
