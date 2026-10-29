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
name: 'attention-gated-prompt-injection-defense'
description: 'Design and evaluate Heads of Security attention gates for prompt injection research, covering token trust boundaries, gate training, adaptive attacks, and utility regression checks.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Attention-Gated Prompt Injection Defense

## Overview

This skill turns Heads of Security (HoS) into a scoped research implementation and evaluation workflow. HoS learns per-head attention gates that separate instruction authority from external data while retaining access to useful content; the source finding evaluates Llama-3.3-70B-Instruct and GPT-OSS-20B. The workflow distinguishes model modifications from hosted API controls and pairs attack resistance checks with task utility measurements.

## When to Use This Skill

- Implementing or reviewing attention-based instruction–data separation in a model whose forward pass is accessible.
- Designing token-level trust annotations for retrieved documents, tool outputs, or mixed-trust prompts.
- Planning gate-only training and reproducible comparisons with an undefended model.
- Evaluating whether security gains survive adaptive attacks and preserve few-shot, multi-turn, and tool-use behavior.
- Assessing whether a hosted model interface supports the required architectural intervention.

## Core Capabilities

1. **Verify evidence and implementation access.** Read the versioned paper, especially Sections 3 and Appendix J, and follow its artifact links before choosing an implementation. As checked on 2026-10-29, the paper links a GitHub repository, but direct retrieval returned 404; released code and checkpoint compatibility were not verified. Recheck access and record exact revisions, licenses, model identifiers, and unresolved dependencies. Do not invent installation commands or checkpoint names.

   Classify the work as reproduction, independent prototype, or feasibility assessment. A hosted completion API without attention hooks cannot implement HoS through prompts, delimiters, or ordinary API fine-tuning. For that case, return an architectural access requirement and an evaluation plan; describe any application-level trust controls as separate measures.

2. **Define and validate trust boundaries.** Build channel labels from application provenance, before serializing model input. The paper uses trusted, untrusted, and template channels; generated tokens are trusted in its setup. Preserve external-data labels even when documents contain role markers or claims of higher authority. Distinguish trusted user instructions from quoted external material inside the same message.

   Verify one channel label per token after chat templating, truncation, batching, and padding. Inspect boundary tokens and ensure attacker-controlled text cannot manufacture template metadata. Carry labels consistently through cached decoding. Include cases with empty external spans, multiple tool messages, and role-like text embedded in data. Treat a label mismatch as an implementation failure, not an attack result.

3. **Specify the attention intervention.** HoS adds `log(sigmoid(g[layer, head, query_channel, key_channel]))` to attention logits before softmax. It freezes base weights and learns gates with preference optimization; GPT-OSS additionally gates its explicit attention sink. These are paper-derived mechanisms, not verified repository APIs. See the [method and appendices](https://arxiv.org/html/2610.37657v1).

   Before implementing, map gate dimensions to the actual attention heads, including grouped-query attention. Preserve causal and padding masks, positional handling, sliding windows, and the model's sink normalization semantics. Do not substitute output-vector scaling or uniform deletion of untrusted tokens for the described operation.

   Compare an equal-gate initialization with the base model on identical inputs within a declared numerical tolerance, including sink entries where applicable. Check finite outputs and gradients, cached versus uncached agreement, and that only intended parameters receive updates. Stop implementation validation if the backend silently bypasses gates or changes masking.

4. **Prepare gate training without overstating completion.** Use preference pairs that retain the legitimate task and external data while contrasting task-following and injection-following responses. Keep training, model selection, and final attack evaluation disjoint by task and attack family where feasible. Record tokenizer, template, data split, seed, optimizer, loss configuration, precision, and trainable parameter inventory.

   Recover model-specific optimization settings from the source before claiming reproduction. For a new architecture, label changes as experimental. If resources permit, run a small gradient and checkpoint round-trip check before scheduling full training. Confirm frozen parameters remain unchanged and save gate configuration together with the channel policy and base-model revision. Never interpret a decreasing training loss as evidence of injection robustness.

5. **Evaluate attacks and utility together.** Use an identical harness, task set, decoding configuration, and tool permissions for the base model and defended variant. Separate static attacks, transferred attacks, and attacks optimized against the defended model. Specify attacker access, writable spans, query or compute budget, stopping rule, and success adjudication before evaluation. Run tool-action tests in a sandbox with synthetic data and mocked consequential actions.

   Report attack success as an explicit numerator and denominator. Count malformed calls, refusals, timeouts, and task failures separately so disabling useful behavior cannot masquerade as security. For adaptive evaluations, report success versus budget and distinguish attempts from cumulative success across attempts. Do not infer adaptive robustness from static attack results.

   Measure clean task completion and task completion under attack. Include tasks that require reading untrusted content, zero-shot and few-shot prompts, later-turn instructions, and valid tool-call generation. Compare utility changes against a tolerance selected for the intended task before inspecting results. Add uncertainty estimates when sample sizes permit; mark small smoke checks as diagnostic rather than benchmark results. Published results are contextual evidence, never acceptance thresholds or locally measured outcomes.

6. **Deliver a bounded research result.** Within 15 minutes, complete source/access triage, an annotated example, an intervention specification, and an evaluation matrix. Run a small smoke check only if compatible local resources are already available; full model training and adaptive attack campaigns are outside this time-boxed outcome.

   Return `ready-for-experiment`, `needs-implementation`, or `blocked` with evidence. If artifacts, hardware, or attention access are missing, provide the exact dependency and a concrete next step. Label unexecuted checks `not run` and report no invented metrics. Passing this workflow means producing a usable research plan or validated smoke report, not certifying a secure deployment. Gate training does not replace application authorization or enforcement of tool permissions.

## Inputs / Outputs

**Inputs**

- Base-model and tokenizer identifiers, revisions, runtime/backend, hardware limits, and access mode.
- Legitimate tasks, external content, provenance-based trust policy, and serialization/template rules.
- Optional compatible gate artifacts or implementation, preference data, and attack/utility harnesses.
- Experiment budget, attack surface and attacker capabilities, utility tolerance, and desired reproduction scope.

**Outputs**

- Evidence record distinguishing verified paper details, accessible artifacts, assumptions, and missing dependencies.
- Token/channel annotation example and a model-specific gate integration specification.
- Training configuration or configuration requirements, with frozen/trainable parameter checks.
- Evaluation matrix covering base and gated models, clean and attacked utility, static and adaptive attacks, and resource budgets.
- Smoke-test results when executed, including sample counts and failures; otherwise explicit `not run` entries.
- A concise readiness decision, limitations, and next experiment. A valid result includes all applicable items without claiming that planned experiments were performed.

## References

- [Heads of Security: Learned Attention Gates Improve Prompt Injection Robustness, arXiv:2610.37657v1](http://arxiv.org/abs/2610.37657v1) — source finding, submitted 2026-10-27.
- [Versioned full text](https://arxiv.org/html/2610.37657v1) — method in Section 3, sink treatment in Appendix D, and training/channel assignment in Appendix J.
- [Paper-linked HoS repository](https://github.com/imperial-aisp/heads_of_security) — returned 404 during authoring on 2026-10-29; verify availability before relying on code or releases.
