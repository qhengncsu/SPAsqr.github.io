---
layout: default
title: FAQ
nav_order: 6
description: "Frequently asked questions about SPAsqr."
has_children: false
---

# **Frequently asked questions**

## Data input and missing values

### How does GRAB handle missing phenotypes and covariates?
{: #how-does-grab-represent-and-handle-missing-phenotype-and-covariate-values }

GRAB treats a blank cell or any of the following values as missing: `NA`, `na`, `NaN`, `nan`, `.`, or `-`.

- **Phenotypes.** GRAB analyzes each trait separately. If a participant is missing one trait, GRAB excludes them from that analysis only. It does not impute missing phenotypes.
- **Covariates.** GRAB replaces a missing value with the mean of that covariate, so the participant stays in the analysis. If the entire column is missing, GRAB reports an error.

---

## Workflow

### Which `--pheno-transform` should I use?

In most cases, you can omit this option. Supply the phenotype file and covariates you used to train the LOCO PGS. GRAB detects the prediction file format and puts the trait on the appropriate scale.

| LOCO PGS from | Transformation GRAB applies |
| ---------------- | ------------------- |
| LDAK-KVIK (PGS on the training trait scale) | GRAB uses the trait as supplied (`raw`) unless you specify a transformation with `--pheno-transform`. |
| REGENIE (PGS on the scale of the covariate-adjusted, unit-variance trait) | GRAB applies `standardize`, then divides the trait by the standard deviation of its residuals from regression on the covariates. If you request a different transformation, GRAB overrides it and issues a warning. |
| No LOCO PGS (`--pred-list` omitted) | GRAB applies `int` by default. |

See [LOCO PGS scale]({{ site.baseurl }}/docs/workflow-1.html#loco-pgs-scale) for details.

