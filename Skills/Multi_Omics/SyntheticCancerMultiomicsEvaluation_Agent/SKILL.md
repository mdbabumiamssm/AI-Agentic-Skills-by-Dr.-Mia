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
name: 'synthetic-cancer-multiomics-evaluation'
description: 'Evaluate synthetic paired cancer RNA-expression and CNV cohorts using patient-disjoint comparisons, biological coherence, and training-data proximity diagnostics.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Synthetic Cancer Multi-Omics Evaluation

## Overview

Evaluate synthetic paired RNA-expression and copy-number variation (CNV) cohorts beyond marginal similarity, using biological relationships and independent patients as evaluation anchors. The source study applies this framework to longitudinal DECIDER high-grade serous ovarian cancer data and finds that stronger fidelity can coexist with stronger training-data proximity. Transfer to other cancers is a proposed application requiring validation, not a demonstrated result.

## When to Use This Skill

- Compare generators that produce complete RNA–CNV pairs for the same synthetic sample.
- Audit whether repeated biopsies or longitudinal samples leak across evaluation splits.
- Assess whether plausible individual modalities retain meaningful relationships when paired.
- Review a synthetic cohort before selecting it for a specified research use.
- Design an evaluation for another cancer while separating methodological adaptation from replication.

## Core Capabilities

1. **Verify the protocol before implementing metrics.**
   Read the versioned full text, especially Section 3 and Appendix C; inspect Appendix B for generation settings and Appendices J–K for controls and sensitivity analyses. Create a protocol ledger mapping each selected metric to its source location, mathematical definition, preprocessing, comparison direction, aggregation, dependencies, and unresolved parameters. Do not implement a metric from the abstract alone. If access or implementation details are missing, return a blocked metric entry with the exact missing information; never silently substitute a familiar formula.

   Verified metric families include CNV G-score agreement; RNA Wasserstein and Jensen–Shannon distances, gene means, variance ratios, PCA-space MMD, and classifier separability; co-expression correlation, edge, module, and hub agreement; cis RNA–CNV slopes; and nearest-neighbor proximity and membership AUC. Appendix C specifies matching and scaling choices. The paper's discriminator uses sample-stratified internal folds, so its score is descriptive. These details require verification together, not independent metric-name matching.

2. **Audit paired data and patient separation.**
   Require a sample-to-patient mapping for every real record. Verify unique sample keys, paired row correspondence, gene identifiers, genome build, missingness, RNA units, and CNV representation. Distinguish continuous ratios from absolute copy counts and discrete event calls; do not apply thresholds across incompatible units. Record any segment-to-gene mapping and purity/ploidy adjustment already performed to avoid double correction.

   Assert that patient sets are disjoint across training, validation, and held-out evaluation. Keep all related biopsies together and identify duplicate profiles. If a patient key is absent or overlaps remain, stop claims of patient-disjoint validity and report the repair needed. Record treatment phase, anatomical site, batch, and repeated-sample counts so cohort composition is visible.

3. **Freeze preprocessing and comparison design.**
   Fit learned transformations and feature selection within the relevant training partition; inner validation must not influence its own fitted preprocessing. Freeze the final configuration before examining held-out outcomes. Preserve paired sample alignment through every transformation, and document retained and excluded genes independently for each modality.

   The paper uses patient-grouped nested evaluation, training-only feature selection, and real training–held-out agreement as a reference. Treat that reference as empirical variability, not an acceptance threshold. For an adaptation, justify split counts and stratification from available patients rather than copying cohort-specific settings. Require generator provenance showing which real records were used for fitting or tuning.

4. **Evaluate biological relationships with explicit comparability checks.**
   Execute only ledger entries whose definitions and data requirements are resolved. Report modality-specific fidelity separately from within-pair coherence. Declare the gene intersection used for paired analyses and its coverage; do not generalize selected-gene results to the whole genome. Record constant-feature exclusions, undefined correlations, empty networks, and insufficient sample sizes as unavailable results with reasons, rather than zero scores.

   For another cancer, prespecify a relevant biomarker panel and rationale before inspecting synthetic results. Label new biological endpoints, alternative thresholds, or patient-grouped discriminator folds as extensions. Any proposed shuffle control must preserve marginal profiles while breaking pairing; verify that marginal diagnostics remain stable and that the intended paired diagnostic responds. Association preservation alone does not establish a causal mechanism.

5. **Interpret proximity alongside fidelity.**
   Keep proximity outputs separate from biological quality scores and inspect sensitivity to cohort size, modality scale, and duplicate records. The study's proximity diagnostics do not establish copying or guarantee privacy. Do not certify a dataset as anonymous, safe to release, or clinically useful from these evaluations. If such a conclusion is requested, identify the separate evidence needed for that specific use.

6. **Produce a reproducible, bounded evaluation report.**
   Start with a feasibility check for the 15-minute execution target. Prefer existing synthetic cohorts and frozen configurations; full generator training is outside this short evaluation unless already feasible. If necessary, run a clearly identified pilot with a prespecified subset and list the omitted analyses. Do not label a pilot, blocked run, or partial report as a completed full benchmark.

   Retain results by fold, generator, seed, metric, and comparison before aggregation. Keep seed variability distinct from between-patient uncertainty; use patient-level resampling if adding uncertainty estimates for repeated observations, and label this as an extension. Avoid selecting a universal winner from incompatible endpoints or treating descriptive differences as statistical significance. Record software versions, random seeds, configuration, elapsed time, and deviations so the run can be repeated.

## Inputs / Outputs

**Inputs**

- Real paired RNA and CNV matrices, sample identifiers, patient identifiers, and a split manifest.
- Synthetic paired matrices for each generator and seed, with training-data provenance and generation settings.
- Feature annotations, measurement scales, preprocessing history, and optional treatment/site/batch metadata.
- The versioned paper, any available implementation, selected endpoints, compute budget, and intended research use.

**Outputs**

- A pairing and leakage audit with explicit pass, fail, or unresolved status.
- A protocol ledger distinguishing verified source definitions, implementation choices, and proposed extensions.
- A tidy results table containing metric, modality, feature scope, fold, seed, comparison, value, units, direction of preference, and unavailable-value reasons.
- A concise report comparing synthetic–held-out results with the real-data reference, alongside training proximity, biological limitations, and reproducibility details.
- A completion status of complete, pilot, partial, or blocked, with remaining work stated explicitly.

A valid completed evaluation requires verified metric definitions, no patient overlap, aligned pairs, interpretable scales, and traceable results for every planned endpoint. A blocked report is a valid diagnostic deliverable, but it is not evidence that the synthetic cohort passed evaluation.

## References

- [Synthesis of Ovarian Cancer Multi-Omics Datasets: A Biology-Informed Evaluation Framework, arXiv:2611.08394v1](http://arxiv.org/abs/2611.08394v1), November 2026. Methodological basis and DECIDER-specific findings.
- [Versioned full text](https://arxiv.org/html/2611.08394v1): Section 3 and Appendix C for evaluation definitions; Appendix B for generation settings; Appendices J–K for controls and sensitivity analyses.
