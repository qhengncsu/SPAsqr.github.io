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

It has to put the trait on the same scale as the LOCO PGS:

| LOCO PGS trained on … | `--pheno-transform` |
| ---------------- | ------------------- |
| RINT $Y$ (recommended) | `int` (default) |
| Raw $Y$, LDAK-KVIK — PGS on the original trait scale | `raw` |
| Raw $Y$, REGENIE — PGS on a standardized scale | `standardize` |
| No LOCO PGS (omitting `--pred-list`) | `standardize` |

