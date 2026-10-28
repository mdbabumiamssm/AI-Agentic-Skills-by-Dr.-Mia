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
name: 'decision-model-document-retrieval'
description: 'Navigate document hierarchies using shared-context decision-model judgments without token generation, and compare retrieval quality, latency, and cost with lexical and embedding retrieval.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Decision Model Document Retrieval

## Overview

Apply the DecisionRAG approach to retrieve passages by navigating document structure with a pretrained decision model. Jointly assess candidate branches in a shared context without generating navigation tokens, allowing retrieval decisions to use language-model capabilities beyond lexical or embedding similarity. Use the workflow below to build a bounded pilot and measure retrieval quality, latency, and cost independently.

## When to Use This Skill

- Relevant evidence is missed by lexical or embedding retrieval before answer generation begins.
- A corpus exposes useful topics, document titles, sections, or other navigable structure.
- A decision-model backend can return candidate judgments without autoregressive token generation.
- A retrieval experiment needs a fair comparison of hierarchical navigation, lexical retrieval, and dense retrieval.
- A paper-grounded feasibility assessment is needed before implementing a larger retrieval system.

## Core Capabilities

1. **Verify the method and execution prerequisites.** Read the versioned paper, especially its method and configuration appendices, and follow any official artifact links it supplies. Verify model access, scoring semantics, request limits, and artifact versions before using implementation-specific settings. Do not invent an SDK, checkpoint, repository, or scoring formula. If the required backend is unavailable, return a feasibility specification with explicit blockers rather than claiming an executed reproduction.

2. **Construct traceable candidates.** Map the available corpus structure into parent-child relationships with stable document, section, and passage identifiers. Preserve passage text and source offsets. Choose branch representations from titles, headings, or supported summaries; record how they were obtained. Treat any newly constructed hierarchy or summary strategy as an experimental choice unless verified in the paper, and account for its preprocessing separately.

3. **Judge candidates together.** Supply the query, relevant conversation context, and identifiable candidate representations in a shared request. Use the verified decision interface to obtain relevance probabilities or supported choice distributions in one forward pass per request. Validate that scores map to the intended candidate IDs. A generated JSON ranking is a different implementation and must be labeled as such; do not describe it as token-free decision scoring.

4. **Navigate under explicit budgets.** Start at the hierarchy root, score the eligible branches, retain a configured subset, and expand their children until passages are reached. Record branch limits, pruning rules, maximum depth, context limits, and stopping conditions. Keep a trace of retained and discarded IDs so missed evidence can be traced to pruning. If candidates exceed a request limit, use the verified batching policy; do not assume scores from separate contexts or hierarchy levels are globally calibrated. Declare any lexical or dense candidate prefilter as a hybrid variant.

5. **Return grounded retrieval results.** Rank the reached passages using the documented scoring or optional reranking procedure. Return up to the requested number of unique passages with source locations and selection traces. Report an empty or incomplete result when warranted, without fabricating evidence. Evaluate generated answers separately if answer generation is requested.

6. **Compare retrieval methods fairly.** Run lexical and embedding baselines on the same corpus snapshot, passage boundaries, query set, and relevance judgments. Record model versions, baseline configuration, candidate budgets, and reranker use. Keep development queries separate from evaluation queries; label any tuning or task-specific training. Inspect shared-context versus independent scoring and pruning-related misses when the available budget permits.

7. **Measure and close a bounded pilot.** Choose a workload that fits the 15-minute execution target with available resources; this target is not a claim about paper-scale reproduction. Report quality using clearly defined metrics appropriate to the labels, such as Recall@k or Hit@k. Without relevance labels, report inspected examples and mark quantitative quality as unavailable. Measure wall-clock retrieval latency with stated boundaries, sample count, hardware, concurrency, and warm-up policy; summarize median and tail latency when sample size supports them. Report request counts, input usage, accelerator time or billed cost separately, and distinguish offline preparation from query-time work. Stop at the time budget and identify incomplete measurements.

## Inputs / Outputs

**Inputs**

- Corpus snapshot, source locations, existing hierarchy metadata, and passage boundaries.
- Queries and optional dialogue history; optional relevance judgments and evaluation split.
- Accessible decision model with a documented non-generative scoring interface.
- Requested result count, navigation and context budgets, runtime budget, and available compute.
- Lexical and embedding baseline configurations, plus optional reranker settings and pricing information.

**Outputs**

- Execution status: completed pilot, partial pilot, or feasibility-only assessment with missing prerequisites.
- Ranked passages containing query ID, passage ID, document ID, source location, text, rank, and scoring context.
- Navigation trace containing candidate IDs, parent IDs, judgments, retained branches, and stopping reason.
- Comparison report separating retrieval quality, online latency, online cost, and offline preparation; mark unavailable values explicitly.
- Reproducibility record containing corpus and model versions, settings, query sample, evaluation definitions, deviations from the paper, and observed failures.

A valid completed pilot has resolvable passage identifiers, a trace consistent with returned evidence, and measurements tied to the actual workload. Paper-reported results are reference observations, never substitutes for local measurements or guarantees of performance.

## References

- [Retriever-Free Retrieval: Navigating Documents with Decision Models, arXiv:2610.35868v1](http://arxiv.org/abs/2610.35868v1) — source finding for DecisionRAG and its retrieval claims.
- [Versioned full text](https://arxiv.org/html/2610.35868v1) — consult Sections 3–4 and the configuration, latency, and cost appendices before implementing paper-specific details.
