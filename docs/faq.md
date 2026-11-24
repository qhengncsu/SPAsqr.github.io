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

If you supply a LOCO PGS, use the same phenotype file and covariates you used for training and leave `--pheno-transform` unset. For example, if you trained on `simu_geno_int.txt` after RINT, pass that file to GRAB as well.

- The **LDAK-KVIK LOCO PGS** is already on the scale of the training phenotype. By default GRAB uses the file as supplied (`raw`). You can override this setting explicitly.
- The **REGENIE LOCO PGS** is on the scale of the covariate-adjusted, unit-variance trait. GRAB standardizes and rescales the phenotype automatically. If you specify a transformation other than `standardize`, GRAB overrides it and issues a warning.
- **Without a LOCO PGS**, GRAB applies RINT by default (`int`).

See [LOCO PGS scale]({{ site.baseurl }}/docs/workflow-1.html#loco-pgs-scale) for how GRAB matches the phenotype to the LOCO PGS from each program.
