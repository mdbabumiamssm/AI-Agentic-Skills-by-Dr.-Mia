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
name: 'scatter-spatial-cell-state-deconvolution'
description: 'Use SCATTER to separate cell-type composition from within-cell-type spatial expression variation using mixed spatial counts, coordinates, and an annotated single-cell reference.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# SCATTER Spatial Cell-State Deconvolution

## Overview

SCATTER estimates cell-type proportions and spot-specific relative expression profiles from mixed spatial transcriptomic measurements and an annotated single-cell reference. It uses spatial context and cell2location-derived proportion regularization to help distinguish changes in tissue composition from expression variation within a cell type. This skill guides implementation verification, input harmonization, fitting, and interpretation with explicit checks for reference mismatch.

## When to Use This Skill

- Analyze capture spots containing multiple cell types when composition alone cannot explain spatial expression patterns.
- Map within-cell-type expression programs across tissue neighborhoods, lesions, or tumor boundaries.
- Compare cell-type proportions, relative expression profiles, and allocated counts as separate analysis outputs.
- Assess whether a single-cell reference adequately represents the tissue before interpreting deconvolved cell states.

## Core Capabilities

1. **Verify the implementation before execution.** Follow the source paper to its official code and documentation; confirm that the repository implements this SCATTER method rather than an unrelated package with the same name. Record the repository URL, release or commit, dependencies, input schema, preprocessing requirements, fitting entry point, and output definitions. The abstract was verified during skill preparation, but an official implementation and executable API were not confirmed; do not invent installation commands, function names, loss terms, or defaults. If verification remains impossible, return an input-readiness report identifying the blocker and do not claim to have run SCATTER.

2. **Harmonize spatial and reference data.** Preserve original count matrices and align gene identifiers, species, spot identifiers, cell annotations, and coordinate rows. Resolve duplicate gene identifiers using a documented rule; record retained genes and exclusions. Check nonnegative finite counts, empty spots or cells, annotation completeness, reference cells per type, and shared marker coverage. Use only implementation-documented transformations for fitting; do not silently substitute log-normalized values for required counts. Keep coordinate units and tissue-section identifiers explicit, and prevent spatial neighborhoods from crossing independent sections.

3. **Prepare proportion regularization.** Obtain or compute cell2location estimates using the verified SCATTER workflow. Align spot and cell-type axes exactly, and confirm the expected conversion from cell2location abundance estimates to the quantities consumed by SCATTER. Do not assume that abundances are already proportions or that arbitrary normalization is interchangeable with the documented procedure. Record reference annotations and prior provenance so incompatible type definitions are detectable.

4. **Fit with local spatial context.** Run the verified implementation to weight reference profiles using spatial neighborhoods and jointly estimate proportions and relative expression while reconstructing spot counts. Record neighborhood settings, selected genes, random seed, stopping criteria, hardware, and fitting diagnostics. Begin with a documented small example or a spatially coherent subset retaining neighborhood context. The 15-minute target applies to a prepared smoke run and validation; full tissue fitting and dependency installation may require longer. Report a partial or blocked run honestly rather than relaxing convergence requirements to meet that target.

5. **Validate reconstruction and count allocation.** Inspect convergence, finite outputs, expected dimensions, nonnegativity, and proportion constraints according to the implementation. Summarize reconstruction residuals by spot, gene, and region, including library-size effects. Use the documented post-fit allocation based on learned gene-specific calibration factors and fitted proportions; do not substitute relative profiles multiplied by proportions as an equivalent output. Check that allocated counts summed across types reproduce retained observed counts within documented numerical tolerance, and report excluded genes or any unassigned component separately.

6. **Check reference mismatch and sensitivity.** Flag absent tissue populations, rare or poorly represented reference types, donor or platform differences, and closely related types with weak distinguishing markers. Compare inferred patterns with observed marker expression and residual maps; a good reconstruction alone does not establish correct cell-state assignment. Where feasible, repeat a small fit with an alternative compatible reference or documented neighborhood setting, keeping comparison axes fixed. Report unstable regions and types instead of forcing biological labels for unsupported populations.

7. **Interpret cell states separately from composition.** Display proportions, relative expression, and allocated counts separately with their units and normalization documented. Allocated counts depend on abundance and sequencing depth and are not direct per-cell expression measurements. Flag estimates for types with negligible local support; use the documented support rule or an explicitly declared analyst threshold. Treat neighborhood associations, regulon patterns, and candidate ligand-receptor rankings as exploratory findings; use biological replicates and spatially appropriate inference for formal comparisons, and do not interpret inferred co-expression as demonstrated signaling or causality.

## Inputs / Outputs

**Required inputs**

- Spatial gene-count matrix, spot identifiers, spatial coordinates, coordinate units, and section identifiers where applicable.
- Annotated single-cell reference with compatible genes and species, cell-type labels, and counts or other representation required by the confirmed implementation.
- Compatible cell2location estimates, or the data and verified workflow needed to generate them.
- Verified SCATTER code and environment, documented preprocessing and fitting settings, and an output destination.

**Optional inputs**

- Histology, tissue-region annotations, donor and batch metadata, marker sets, and independent validation measurements.
- Prespecified neighborhood comparisons, pathway or regulon resources, or ligand-receptor resources with their versions.

**Expected outputs after a completed fit**

- Spot-by-cell-type proportions with explicit axis identifiers and normalization semantics.
- Spot-by-cell-type-by-gene relative expression profiles and separately labeled allocated observed counts, using the implementation-supported storage layout.
- Gene-specific calibration factors when exposed by the implementation, reconstruction summaries, convergence diagnostics, and flags for low support or reference mismatch.
- Spatial summaries separating composition from within-type expression and a provenance record covering data filters, code version, parameters, seeds, and elapsed time.

**Completion criteria:** All required axes align; outputs satisfy documented constraints; reconstruction and allocation checks are reported; unsupported estimates are flagged; and a reproducible run record accompanies the results. A readiness report is a useful fallback but is not a successful deconvolution result. Do not present the paper's benchmark results as accuracy guarantees for a new dataset.

## References

- [Resolving cell states in spatial transcriptomics with SCATTER, arXiv:2611.34798v1](http://arxiv.org/abs/2611.34798v1). Source finding and method description; preprint submitted November 25, 2026. Implementation verification is required before execution.
- [cell2location official repository](https://github.com/BayraktarLab/cell2location). Related dependency for abundance estimation; consult the verified SCATTER implementation for compatibility and conversion requirements.
