---
layout: default
title: FAQ
nav_order: 6
description: "Frequently asked questions about SPAsqr."
has_children: false
---

# **Frequently asked questions**

## Data input and missing values

### How does GRAB represent and handle missing phenotype and covariate values?

- **Missing tokens.** A phenotype / covariate cell is read as missing if it
  is `NA`, `na`, `NaN`, `nan`, `.`, `-`, or blank.
- **Phenotype** — handled per trait: a subject missing a trait is dropped
  from *that* trait only (never imputed) and still contributes to others.
- **Covariate** — mean-imputed per column; no subject is dropped for a
  missing covariate (a fully-missing column is an error).

---

## Workflow

### Which `--pheno-transform` should I use?

In most cases, none. Pass GRAB the same phenotype file and covariates that were used to train the LOCO PGS; GRAB detects the LOCO file format and puts the trait on the scale of the LOCO PGS:

| LOCO PGS from | Transformation GRAB applies |
| ---------------- | ------------------- |
| LDAK-KVIK — PGS on the scale of the training trait | `raw` by default (trait used as supplied); an explicit `--pheno-transform` is honored |
| REGENIE — PGS on the scale of the covariate-adjusted, unit-variance trait | always `standardize`, followed by division by the standard deviation of the trait's residuals on the covariates; other values are overridden with a warning |
| No LOCO PGS (omitting `--pred-list`) | `int` by default |

See [LOCO PGS scale]({{ site.baseurl }}/docs/workflow-1.html#loco-pgs-scale) for details.

