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
name: 'gepa-curator-agent-improvement'
description: 'Curate agent experiences and retrieve reliable prior successes for reflective prompt optimization, with failure attribution and held-out evaluation.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# GEPA-Curator Agent Improvement

## Overview

Use GEPA-Curator to improve the evidence supplied to reflective prompt optimization. Its task curation and success retrieval help distinguish prompt weaknesses from unreliable execution feedback. This skill adds a bounded local audit and evaluation workflow; it does not promise improvement within the time limit.

## When to Use This Skill

- Prompt revisions repeatedly fail to improve an agent's task completion.
- Rollout logs mix tool errors, inconsistent judgments, and potentially correctable behavior.
- A team needs an auditable experiment comparing an existing prompt with a proposed revision.
- Historical traces are available but their suitability as learning examples is unclear.

## Core Capabilities

### Published method

The following summarizes Sections 3–4 and Appendix B.2 of the [paper](https://arxiv.org/pdf/2610.37651v1).

1. **History-based selection.** Maintain per-task traces and judgments. Run underexplored, uncertain, or prompt-sensitive tasks; skip only with strong evidence of low learning value. Reconsider skips each epoch.
2. **Failure attribution.** Distinguish mastery, missing capabilities, persistent unproductive failures, evaluator inconsistency, and infrastructure errors. Rejected mutations can demonstrate sensitivity without establishing current mastery.
3. **Success retrieval.** Jointly decide selection and retrieve one same-task success corroborated by another independent clean success. Use no anchor for skipped tasks or insufficient evidence.
4. **Reflective optimization.** Supply selected evidence as optimizer metadata, then append resulting traces and outcomes to history. The reported setup uses a greedy lineage and selects its final checkpoint by validation performance.

### Local execution workflow

These operational conventions are recommendations for applying the method, not claims about an official implementation or API.

1. **Establish the experiment contract.** Record the baseline prompt, actor and optimizer versions, harness revision, evaluator configuration, task split, reset procedure, rollout budget, and deadline. Keep development, validation, and final test task IDs disjoint. Set the minimum useful improvement and regression tolerance before viewing candidate results.
2. **Verify the integration.** Read the complete paper, especially Appendix B.2, and inspect the installed optimizer's actual interfaces. The linked GEPA repository is an upstream optimizer, not proof that GEPA-Curator is installed. Identify where task sampling, reflection metadata, and history persistence enter the local code. Do not invent package names, command flags, or callable signatures. If the integration is unavailable, complete an evidence audit and report execution as blocked.
3. **Normalize evidence.** Assign immutable task, rollout, and prompt IDs. For each attempt, retain the initial environment identity, tool observations, terminal state, evaluator output, exception information, and mutation acceptance status. Preserve raw evidence alongside summaries. Mark missing resets and configuration changes so comparisons do not silently cross incompatible conditions.
4. **Produce a reviewable curation record.** Apply the published method above and the full appendix instructions. For every decision, record task ID, run/skip choice, supporting trace IDs, attribution, uncertainty, revisit epoch, selected anchor ID or null, corroborating IDs, and a separate anchor justification. This is a local audit schema; adapt it explicitly to the verified implementation. Validate that all referenced IDs exist and belong to the declared task and training split.
5. **Run a bounded optimization trial.** Assemble the eligible batch and invoke the verified optimizer interface. Require each proposed edit to identify the observed behavior it is intended to change. Keep the baseline and candidate prompts separately addressable. Put infrastructure repair proposals in a separate record; changing the harness during a prompt comparison invalidates that comparison. If no eligible tasks remain or the budget expires, return the audit instead of looping indefinitely.
6. **Evaluate without leakage.** Compare baseline and candidate under matched environment resets and evaluator settings. Use validation to choose a candidate, then freeze it before final test evaluation. Never retrieve validation or test traces into optimization context. Report task counts, successes, failures, infrastructure errors, repeat count, resource use, and uncertainty. Keep infrastructure errors visible in end-to-end results; label any additional analysis that excludes them and show its denominator.
7. **Finish with an evidence-based decision.** Return retain-baseline, candidate-supported, inconclusive, or blocked, with the reason and supporting measurements. If the harness or judge changes, require a fresh baseline comparison. An apparent gain from a small or noisy sample is inconclusive unless it meets the declared decision rule.

## Inputs / Outputs

**Inputs**

- Baseline prompt and explicit optimization objective.
- Task definitions, stable IDs, and declared dataset splits.
- Rollout histories, evaluator evidence, and environment diagnostics.
- Local optimizer interface, versioned configuration, and available execution budget.

**Outputs**

- Curation ledger with evidence links and unresolved attribution questions.
- Selected training batch and validated anchor mapping, when available.
- Candidate prompt and explained diff, if optimization ran.
- Evaluation report with counts, denominators, uncertainty, configuration, and final decision.
- Explicit missing prerequisites and next actions when execution cannot be completed.

**Fifteen-minute completion contract:** Use existing logs and an already configured harness for a small trial. Reserve time to return the ledger and status report before the deadline. A valid audit has resolvable evidence IDs, explicit split boundaries, justified decisions, and honest execution status. If training or held-out evaluation cannot finish, return partial measurements and label improvement unverified; a completed audit is not a completed benchmark reproduction.

## References

- [GEPA-Curator: Curated Learning for Self Improving Agents — supplied source finding, arXiv:2610.37651v1](http://arxiv.org/abs/2610.37651v1).
- [Full paper PDF](https://arxiv.org/pdf/2610.37651v1) — Sections 3–4 and Appendix B.2; use the PDF for appendix prompt details omitted from HTML rendering.
- [GEPA upstream implementation](https://github.com/gepa-ai/gepa) — verify the installed revision and integration points before execution; this reference does not establish a released GEPA-Curator implementation.
