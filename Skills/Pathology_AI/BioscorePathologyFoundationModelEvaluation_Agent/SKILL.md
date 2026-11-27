<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal AI Agentic Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA
-->



<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->

---
name: 'bioscore-pathology-foundation-model-evaluation'
description: 'Evaluate cell-identity representations in pathology image encoders using sparse autoencoders and matched spatial transcriptomics, with BioScore comparisons and per-identity interpretation.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# BioScore Pathology Foundation Model Evaluation

## Overview

BioScore evaluates whether sparse features of pathology image encoders track marker-defined cell identities in matched spatial transcriptomics. It complements downstream classification and retrieval by examining represented biology across a reference panel of 52 identities; transcriptomics supplies the evaluation reference, while sparse autoencoder training remains unsupervised.

## When to Use This Skill

- Compare the biological representations of frozen pathology encoders on a common tissue cohort.
- Investigate which cell identities an encoder captures or misses despite similar downstream task performance.
- Review a proposed BioScore reproduction for spatial alignment, reference coverage, and comparison validity.
- Produce a feasibility report or evaluate cached outputs within a short compute budget.

## Core Capabilities

1. **Establish a reproducible evaluation contract.** Record dataset identifiers, section and patient membership, encoder checkpoints, preprocessing, embedding extraction, software revisions, and resource limits. Distinguish exact reproduction from a cohort-specific adaptation before interpreting scores. Read the linked methods and inspect available implementation/configuration files; do not invent a BioScore command-line interface. The linked `dictionary_learning` package implements SAE training, not the complete benchmark.

2. **Validate aligned inputs.** Require a shared patch identifier connecting image coordinates, embeddings, and gene counts. Check coordinate transforms, tissue coverage, gene symbols, missing panels, duplicate records, and finite values. Reject unmatched rows rather than silently joining by row order. Compare models on identical tissue footprints and explicitly record every exclusion.

3. **Recover sparse feature groups.** The published recipe uses 8-fold expansion, BatchTopK with k=64, and K=3000 groups formed by weighted spherical clustering of normalized decoder columns. Group activity is its fraction of total patch activation. Verify checkpoint compatibility and reconstruction diagnostics before scoring. Flag zero-total activation and inadequate active-feature counts; do not silently alter K to accommodate failed training.

4. **Construct the biological reference.** Apply the paper's panel-compatibility, presence, and localization rules; implement signature normalization and identity aggregation from equations 3–4. Preserve the exact marker lists and eligible sections. Missing markers are unmeasured, not negative expression. Keep transcript-derived targets out of encoder and SAE optimization. Freeze reference curation before model comparison and document any cohort-specific changes.

5. **Compute specificity and aggregate.** Use Pearson alignment and the equation-5 impurity penalty; maximize Q=rho−3U over groups for each identity, then average within and across the five lineage classes. Preserve the paper's correlation masks and shared-section rules. Treat undefined correlations explicitly, retaining denominators and exclusion reasons. Verify intermediate values on a small fixture before processing cached cohort outputs.

6. **Compare encoders consistently.** Report five SAE seeds (40–44), score mean/SD, and identity-paired Wilcoxon tests with Holm correction. Record comparison family, sample count, effect direction, and seed aggregation. Do not treat patches as independent biological replicates. Label single-seed or reduced-panel runs as preliminary; avoid claiming equivalence from a nonsignificant comparison.

7. **Interpret with limits.** Review high- and low-scoring identity maps alongside H&E and reference coverage. The original SAE training includes scoring-section images; describe this overlap. Separate weak model alignment from poor registration, sparse transcripts, absent markers, or unstable SAE fits. A cohort-specific score does not establish causal understanding, cell-level segmentation accuracy, or clinical utility. Changes to tissues, panels, identities, or preprocessing require renewed comparability checks.

## Inputs / Outputs

**Inputs**

- Matched H&E and spatial transcriptomics, or validated patch-level embeddings and counts with spatial coordinates.
- Section/patient manifest, gene-panel definitions, marker signatures, identity-to-lineage mapping, and eligibility masks.
- Frozen encoder identifiers, compatible SAE checkpoints or training configuration, random seeds, and compute budget.
- A verified scoring implementation or enough method detail to identify missing implementation steps.

**Outputs**

- A run manifest recording versions, parameters, data overlap, exclusions, and deviations from the reference protocol.
- A per-identity table containing model, seed, identity, eligible sections/patches, selected group, correlation, impurity, quality, and status.
- An aggregate comparison containing lineage scores, overall score, seed variability, paired differences, and corrected significance where justified.
- A concise interpretation of supported findings, coverage gaps, unresolved checks, and reproduction status; optional spatial overlays for inspection.

**Time-bounded execution and acceptance**

Use the 15-minute target for input validation and a cached-data scoring pass when artifacts are ready. If full training or scoring cannot finish, return a clearly labeled feasibility report listing verified inputs, blockers, and the required next computation; do not report the benchmark as completed. A valid scored result requires traceable inputs, consistent comparison coverage, auditable intermediate quantities, and disclosed deviations. A feasibility report must contain no fabricated scores.

## References

- [A gene informed benchmark for BioScoring of pathology foundation models, arXiv:2611.34930v1](http://arxiv.org/abs/2611.34930v1) — source finding and benchmark motivation.
- [Full methods, equations, and appendices](https://arxiv.org/html/2611.34930v1) — reference protocol for grouping, signatures, scoring, and comparisons.
- [dictionary_learning](https://github.com/saprmarks/dictionary_learning) — SAE library cited by the paper; inspect and pin a compatible revision before use.
