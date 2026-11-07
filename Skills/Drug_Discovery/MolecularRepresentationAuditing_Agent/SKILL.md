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
name: 'molecular-representation-auditing'
description: 'Audit molecular fingerprints and learned embeddings for collisions, structural and 3D similarity, and endpoint relevance before clustering, dataset splitting, or virtual screening.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Molecular Representation Auditing

## Overview

Audit molecular representations against the information needed for a specific clustering, dataset splitting, or virtual screening task. The source study motivates examining representation degeneracy, graph structure, three-dimensional shape, and endpoint similarity beyond downstream prediction scores. This skill defines a reproducible local audit; its suggested checks are not an exact reproduction of the paper or universal acceptance criteria.

## When to Use This Skill

- Selecting fingerprints, physicochemical descriptors, or learned embeddings for a molecular workflow.
- Investigating identical vectors for distinct molecules, unexpected nearest neighbors, or unstable clusters.
- Checking whether train/test partitions separate structurally related molecules adequately for the intended evaluation.
- Assessing whether similarity rankings preserve shape, stereochemistry, or relevant ADMET measurements.
- Reviewing screening ties, representation changes, or claims that a learned embedding improves on classical baselines.

## Core Capabilities

1. **Define scope and evidence before computation.** Record the intended use, dataset version, available labels and conformers, compute budget, and task-specific acceptance criteria. Retrieve the versioned full paper and relevant Methods before attributing exact metrics, parameters, or conclusions to it; if retrieval fails, label evidence as abstract-only and continue with explicitly local checks. Do not transfer library-wide findings to the user's collection.

2. **Build a bounded, reproducible audit cohort.** Use a fixed seed and retain selected molecule IDs, input hashes, sampling rules, and exclusions. Within a 15-minute budget, prefer cached embeddings and conformers and a bounded set of molecules and pairs; record the actual size and coverage. Avoid unrestricted all-pairs calculations, model training, or large downloads. Reserve time for reporting, and mark incomplete axes as unassessed rather than passed.

3. **Separate identity from representation equality.** Parse structures and document salt handling, charge, aromaticity, tautomer policy, and stereochemistry retention. Preserve original IDs and structures alongside standardized identities. Count duplicate records separately from distinct molecules that share a vector; inspect stereoisomers and normalization-induced merges explicitly. Report invalid structures and representation failures, and compare methods on a shared valid cohort while retaining method-specific coverage counts.

4. **Construct fair representation baselines.** Compare candidate embeddings with at least one suitable fingerprint and a physicochemical descriptor baseline when available. Record toolkit versions, fingerprint family, radius or path settings, bit length, binary/count mode, and chirality flags. For embeddings, record checkpoint identifier or hash, tokenizer, pooling, vector dimension, inference settings, and known training overlap. Record descriptor columns, missing-value handling, and scaling; fit any learned preprocessing on training data only when evaluating a split. Flag constant, nonfinite, or zero-norm vectors before similarity calculations.

5. **Measure collisions and near-collisions.** Group exactly equal vectors and verify equality after any hash-based grouping. For N distinct standardized molecules with valid vectors and collision-group sizes g greater than one, report the number of groups, sum(g), sum(g)/N, and maximum group size. These count molecules participating in collisions, not only redundant entries. Treat floating-point near-equality as a separate analysis with an explicit norm, tolerance, and dtype; never silently round embeddings into exact collisions. Inspect representative pairs, distinguishing possible bit-folding collisions from information omitted by the encoding. A sample with no collisions does not establish library-wide uniqueness.

6. **Audit structural and three-dimensional neighborhoods.** Specify the similarity function for each representation, such as Tanimoto for binary fingerprints and cosine for nonzero embeddings, and retain its parameters and direction. Compare rankings on the same sampled pairs against an independent graph comparison, with explicit atom/bond matching rules, normalization, and timeout behavior. If using a maximum-common-substructure proxy, identify it as a proxy rather than the paper's exact metric. For 3D checks, document conformer generation, protonation, stereochemistry, alignment, shape score, and ensemble aggregation. Compare both broad pair samples and local nearest neighbors using rank association and top-k overlap, with k and tie handling fixed in advance. Report conformer failures, graph timeouts, pair counts, and whether distances were reversed before correlation; do not infer 3D fidelity from 2D agreement alone.

7. **Evaluate endpoint relevance and uncertainty.** Analyze each measured ADMET endpoint separately, preserving assay context, units, transformations, censoring, missingness, and replicate handling. Compare endpoint differences or class agreement among retrieved neighbors against the same baseline cohorts. Inspect activity cliffs and label disagreements within collision groups. Describe sampling variability using molecule-level resampling when feasible; shared-molecule pairs are not independent observations. Record known pretraining overlap and unknown provenance. Missing endpoint data permits no conclusion about endpoint fidelity.

8. **Translate findings into use-specific decisions.** For clustering, inspect size distributions, singleton prevalence, seed or threshold sensitivity, and independently measured structural neighbors split between clusters. For dataset splitting, check exact identities, collision groups, and close structural pairs across partitions; report affected test-molecule counts and denominators, without treating cluster separation as proof of independence. For virtual screening, retain all identities in tied groups, report ranking ambiguity, and assess known-active retrieval only where valid labels exist. Propose task-appropriate mitigations such as revised encoding parameters, alternative split constraints, or independent graph/3D reranking, and rerun the affected checks before claiming improvement.

9. **Issue an auditable disposition.** Assign each representation and intended use a status of suitable within tested scope, conditional, unsuitable, or insufficient evidence. Tie the decision to predeclared criteria, observed failures, coverage, and uncertainty; do not invent a universal score threshold or declare one representation universally best. If criteria are absent, provide a descriptive comparison and a conditional recommendation. A valid report contains the reproducibility manifest, numerical results with denominators, concrete failure examples where observed, and explicit unassessed work; successful workflow completion does not mean the representation passed.

## Inputs / Outputs

**Inputs**

- Required: molecule IDs with SMILES or molecular structures, intended use, and candidate representations or their generation configurations.
- Optional: existing vectors, checkpoints, conformers, endpoint measurements with assay metadata, cluster/split assignments, screening queries and relevance labels, and acceptance criteria.
- Execution context: available toolkit/model versions, hardware, runtime and memory limits, random seed, and permitted output location. If essential data are absent, produce an insufficient-evidence report identifying the missing inputs instead of fabricated results.

**Outputs**

- A concise audit report with scope, sampling coverage, identity policy, per-use disposition, limitations, and prioritized mitigations.
- A reproducibility manifest containing input and checkpoint identifiers, software versions, preprocessing, representation and distance parameters, seeds, sampled IDs, and commands or pseudocode sufficient to rerun the audit.
- A comparison table with representation, cohort size, failures, collision counts and rates, available structural/shape/endpoint checks, runtime, and assessment status. Use explicit unavailable or unassessed values rather than zeros for missing checks.
- An evidence table of collision members, discordant neighbor pairs, cross-partition neighbors, and screening ties as applicable, preserving molecule IDs and the relevant measurements.

Return these as structured report sections unless separate artifacts are requested. The 15-minute target applies to the bounded audit and honest reporting of incomplete work, not exhaustive validation of an ultra-large library.

## References

- [On the limitations of molecular representations, arXiv:2611.08420v1](http://arxiv.org/abs/2611.08420v1) — source finding; consult the [versioned full text](https://arxiv.org/html/2611.08420v1) for study-specific methods and scope.
