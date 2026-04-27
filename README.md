# BRCA-Omics-Combination-Recommendation

A disease-specific, cross-algorithm and multi-metric framework for identifying optimal multi-omics combinations for breast cancer subtyping.

## Overview

This repository contains the code used in our study on breast cancer-specific multi-omics combination recommendation.  
The project systematically evaluates admissible multi-omics combinations across multiple integration algorithms and derives final recommendations using a cross-algorithm normalized ranking framework.

The workflow is designed to address the following question:

> Which omics combinations are most informative for breast cancer subtyping, beyond default combinations commonly used in previous studies?

## Study design

The analytical workflow includes the following major steps:

1. Multi-omics data preprocessing
2. Construction of admissible omics combinations
3. Multi-omics integration using representative algorithms
4. Evaluation using individual metrics
5. Composite scoring within each algorithm
6. Cross-algorithm normalized ranking
7. Final recommendation of BRCA-specific omics combinations
8. Omics contribution and sensitivity analyses

## Supported omics layers

The framework considers the following omics types:

- CNV: copy number variation
- MET: DNA methylation
- RNA: mRNA expression
- MIR: miRNA expression
- PRO: proteomics
- SNP: somatic mutation

