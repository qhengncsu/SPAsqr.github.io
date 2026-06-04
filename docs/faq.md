---
layout: default
title: FAQ
nav_order: 8
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

It has to match the transform used during PGS construction:

| PGS trained on … | `--pheno-transform` |
| ---------------- | ------------------- |
| INT $Y$ (recommended) | `int` (default) |
| Raw $Y$ — LDAK-KVIK / REGENIE internally standardize | `standardize` |
| No PGS (omitting `--pred-list`) | `standardize` |
