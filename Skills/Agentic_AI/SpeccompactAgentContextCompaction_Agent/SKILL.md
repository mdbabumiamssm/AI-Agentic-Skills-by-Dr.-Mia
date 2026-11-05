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
name: 'speccompact-agent-context-compaction'
description: 'Reproduce and evaluate SpecCompact context compaction using a fine-tuned small drafter, parallel summary slots, and target-model semantic verification, measuring retention and latency against ordinary summarization.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# SpecCompact Agent Context Compaction

## Overview

SpecCompact accelerates agent-history summarization through parallel slot drafting, target-model semantic verification, and local repair. This skill checks implementation readiness and evaluates whether reduced compaction latency preserves the information needed to continue a task; semantic acceptance does not guarantee lossless summaries or unchanged task performance.

## When to Use This Skill

- Context compaction interrupts a long-running agent with noticeable delays.
- Reproducing SpecCompact on an existing agent history or supported benchmark.
- Comparing a fine-tuned drafter with ordinary target-only summarization.
- Investigating whether verification, repairs, or repeated context processing erase drafting savings.
- Checking retention before integrating a faster compactor into an agent.

## Core Capabilities

1. **Confirm implementation availability.** The paper-linked [repository](https://github.com/shawnyin128/Speculative-Compaction) was accessible during authoring and documents inference, training, and checkpoint links. Recheck availability at execution time; record the commit, checkpoint revision, license, dependency lock, and endpoint capabilities. Repository discovery is not evidence that a local run succeeds.

2. **Prepare reproducible inputs.** Freeze the history, compaction boundary, slot schema, token budgets, model identifiers, and decoding settings. Keep system instructions, current user input, and agent setup outside the replaceable history. Mark tool results and quoted content as data; preserve message roles and ordering. Reject oversized inputs explicitly rather than silently truncating them.

3. **Reproduce the semantic loop.** Use the paper's slots: Goal, Constraints and Preferences, Progress, Key Decisions, Next Steps, and Critical Context. A fine-tuned small drafter proposes lines from the history and each slot's accepted prefix. The target checks support and consistency, accepts a proposal or generates a replacement without the rejected line, then updates that prefix. Run slots concurrently; discard buffered continuations after a repair changes their prefix. Assemble completed slots in schema order. See [method details](https://arxiv.org/html/2611.06876v1).

4. **Use the released execution interface.** Inspect the pinned README before running commands. The documented client requires Python 3.12 or 3.13, separately running model servers, a draft completions endpoint, and target chat completions with token log probabilities. Install the inference group separately from training. Configure endpoint URLs, served model names, and tokenizer paths in `compact.yaml`; reference credentials through environment variables. The documented commands are:

   ```bash
   uv sync --locked --group inference
   uv run --locked --group inference speccompact compact --config compact.yaml --print-config
   uv run --locked --group inference speccompact compact --config compact.yaml --input context.txt --output outputs/compacted.json
   ```

   Run these from an authorized checkout using new output paths. Confirm the pinned version still supports these arguments. The JSON result contains `text`, `metadata`, and `configuration`; detailed call recording can expose sensitive history, so enable it only within the permitted data boundary. [CLI documentation](https://github.com/shawnyin128/Speculative-Compaction#readme).

5. **Separate smoke testing from training.** Prefer a released drafter paired with the selected target for the short workflow. If training is needed, use target-generated slot supervision, separate trajectories across training and evaluation, mask prompt tokens from the training loss, and record the recipe. Treat an unfine-tuned drafter or substituted verifier as an adaptation rather than an exact reproduction. Retrieve the paper's verifier probability rule and termination rules instead of inventing equivalent-looking defaults.

6. **Measure quality and execution cost together.** Use the paired evaluation protocol below. Shared-prefix caching and line-aware stopping are method optimizations; verify their actual server behavior. A high draft acceptance rate alone cannot establish retention, because an entirely omitted fact never reaches verification.

### Fifteen-minute workflow

This time bound covers a small smoke test with prepared models and accessible histories, not model downloads, fine-tuning, or full benchmark reproduction. The following schedule is an operational recommendation for this skill.

1. **Minutes 0–3: readiness.** Check the pinned implementation, model/checkpoint pairing, server health, tokenizer compatibility, and output location. Define a short held-out history and its critical facts. If prerequisites are missing, return a readiness report with `status: blocked`; do not claim successful execution.
2. **Minutes 3–6: configure.** Validate the effective configuration, including verifier method and threshold, concurrency, per-slot bounds, and total summary budget. Record which structural lines bypass verification. Set a timeout that leaves time for inspection.
3. **Minutes 6–12: paired run.** Compact the same history with SpecCompact and an ordinary target-only summary. Use the same target, preserved-message boundary, and output budget. Include a target-only structured-slot comparator when time permits to distinguish schema effects from speculation. Log failures and retries rather than dropping them.
4. **Minutes 12–15: inspect and report.** Check schema completeness, token budget, unresolved requests, exact identifiers, constraints, completed actions, and next steps. Compare both summaries against the original history. Report timing, retention, missing evidence, and whether a wider evaluation is warranted.

If a run exceeds its timeout, retain diagnostic output and report `status: incomplete`. A valid report can document a blocker, but it does not satisfy the successful-execution outcome. Never replace a live agent's context merely because the smoke test returned syntactically valid JSON.

### Paired evaluation protocol

Use these evaluation choices as local experimental design, not as claims about the published implementation.

- **Controls:** Pair histories and continuation tasks across methods. Record hardware, resource allocation, serving software, precision, concurrency, cache state, prompt sizes, output sizes, and compaction frequency. Include the drafter's resource cost. Keep cold-cache and warm-cache results separate.
- **Retention:** Before inspecting summaries, annotate task-relevant facts and their source spans. Compute supported required facts retained divided by required facts, alongside unsupported statements, contradictions, and critical omissions. Preserve exact paths, IDs, numerical values, and explicit user restrictions where continuation requires them. Use independent review for ambiguous semantic matches.
- **Verification overhead:** Count proposed, accepted, repaired, structurally bypassed, and invalid lines separately. Define semantic acceptance as accepted verified proposals divided by all semantically verified proposals; report unavailable denominators as `N/A`. Measure verifier and repair call latency, token usage, queueing, and repeated prefill where observable.
- **Compaction latency:** Measure elapsed wall time from request submission to validated assembled summary, including orchestration, retries, and final checks. Parallel role times overlap; do not sum them and label the result wall time. Compute paired speedup as baseline wall time divided by SpecCompact wall time, with the aggregation method stated.
- **End-to-end latency:** Resume the same held-out task from each summary in an isolated evaluation environment. Measure total time through task termination, including all compactions and any extra actions caused by lost context. Report task success separately from summary fidelity and compaction speed.
- **Ablations:** Where supported by the pinned implementation, compare serial versus parallel slots, caching disabled versus enabled, and buffering disabled versus enabled. Alter one factor at a time. Record any fallback that changes semantic behavior.
- **Decision rule:** Agree on acceptable critical omissions, retention loss, and task-success loss before comparing results. Require measured latency improvement within those limits before recommending adoption. Report sample count, per-case results, failures, and uncertainty; a small smoke test is not a benchmark claim.

## Inputs / Outputs

**Inputs**

- Conversation histories with message roles, task boundaries, and permission to process them on the selected endpoints.
- Target and compatible fine-tuned drafter endpoints, tokenizer revisions, and compute limits.
- Pinned implementation, effective configuration, summary budget, and baseline prompt.
- Held-out continuation tasks, fact annotations, timing budget, and quality acceptance criteria.

**Outputs**

- Structured compacted text with its effective configuration and execution metadata.
- Availability manifest identifying code, models, dependencies, and any reproduction deviations.
- Per-run comparison of retention, unsupported content, critical omissions, verification/repair cost, compaction wall time, and end-to-end task latency when measured.
- A concise decision report containing `status` (`success`, `blocked`, `incomplete`, or `failed`), evidence, limitations, and the next justified action. Use `not measured` for unavailable metrics; never substitute paper results for local measurements.

Successful smoke-test output requires both paired summaries, checked budgets and required facts, recorded wall times, and an explicit quality comparison. Full reproduction additionally requires verified configuration parity and the broader held-out evaluation; training and benchmark runs can extend beyond the short workflow.

## References

- [SpecCompact: Accelerating LLM Agent Context Compaction with Slot-Parallel Semantic Speculative Verification, arXiv:2611.06876v1](http://arxiv.org/abs/2611.06876v1). Source finding; reported acceleration is accompanied by model- and task-dependent downstream trade-offs.
- [SpecCompact full text](https://arxiv.org/html/2611.06876v1). Consult Sections 3–4 and the prompt appendix for exact method and experimental settings.
- [Speculative-Compaction implementation](https://github.com/shawnyin128/Speculative-Compaction). Paper-linked source, runtime interface, configuration documentation, training recipe, and checkpoint links; pin a revision before execution.
