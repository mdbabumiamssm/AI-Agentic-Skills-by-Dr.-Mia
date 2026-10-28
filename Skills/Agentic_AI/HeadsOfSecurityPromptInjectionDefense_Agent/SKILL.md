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
name: 'heads-of-security-prompt-injection-defense'
description: 'Plan, implement, and evaluate learned attention gates for prompt injection defense, with explicit trust boundaries, verified architecture support, adaptive attacks, and utility checks.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Heads of Security Prompt Injection Defense

## Overview

This skill supports research implementation and evaluation of Heads of Security (HoS), which learns per-head attention gates to separate trusted instructions from untrusted data. It establishes implementation readiness, maps trust boundaries, and plans controlled training and evaluation for agents consuming external content. Robustness is an empirical property of a tested configuration, not a guarantee against prompt injection.

## When to Use This Skill

- Investigating model-level defenses for agents that read retrieved documents, webpages, or tool outputs.
- Reproducing or adapting HoS with access to model weights and attention internals.
- Auditing token-level trust annotations across chat templates and tool interactions.
- Comparing injection resistance with benign task performance under matched evaluation conditions.
- Determining whether a deployment supports attention interventions before committing training resources.

## Core Capabilities

1. **Verify implementation availability and scope.** Read the versioned paper and follow its implementation links before proposing commands. As checked on 2026-10-28, the paper links `imperial-aisp/heads_of_security`, but fetching that GitHub URL returned 404; repository contents, installation commands, checkpoints, and runnable entry points were not verified. Recheck availability and distinguish inaccessible code from evidence that no implementation exists. Record the source version, repository revision, license, checkpoint provenance, and unresolved dependencies. Do not invent package names, configuration keys, or command-line interfaces.

2. **Establish architecture requirements.** Require editable attention computation and gradient access for gate training; an ordinary hosted text-generation API is insufficient. The paper evaluates Llama-3.3-70B-Instruct and GPT-OSS-20B. Its gates add `log(sigmoid(g[query_channel, key_channel]))` to each head's attention logits before softmax, alongside the existing mask. It uses trusted, untrusted, and template channels; GPT-OSS additionally requires a key-only attention-sink channel. Verify the chosen implementation's handling of grouped-query attention, sinks, causal and padding masks, sliding windows, precision, caching, and fused attention kernels before execution. Treat other architectures as unvalidated adaptations.

3. **Annotate trust boundaries from provenance.** Build channel assignments from application-controlled message structure, not instructions found inside external text. Preserve alignment after template rendering, tokenization, truncation, packing, and cached generation. Explicitly classify quoted data inside user messages, tool-result bodies, template tokens, and generated tokens. The paper typically treats model generation as trusted; record this convention without treating generated claims as evidence that external content has become trustworthy. Ensure attacker-controlled role delimiters cannot change channel assignments. Stop implementation work if token provenance cannot be represented reliably.

4. **Prepare gate training with utility controls.** The paper freezes base weights and trains gate parameters using Direct Preference Optimization with preferred instruction-following responses and rejected injection-following responses. Recover initialization, objective settings, reference-policy handling, optimizer, data processing, and checkpoint selection from verified source material before claiming reproduction. For a new experiment, label choices as adaptations. Split examples by underlying task and attack family to reduce leakage; reserve benign validation examples requiring useful reasoning over untrusted data. Record trainable parameter names, seeds, resource limits, and stopping criteria. Small trainable parameter counts do not eliminate base-model memory and backward-pass costs.

5. **Check implementation invariants before scaling.** Compare outputs with gates bypassed against the unchanged base model. Test that equal gate biases over every effective key, including any explicit sink, preserve attention normalization within a stated numerical tolerance. Confirm masked positions remain masked, gate gradients are finite, frozen weights remain unchanged, and cached decoding retains correct channel labels. Test batches with padding, long tool results, mixed trusted and untrusted spans, and attacker-supplied template-like strings. Compare supported attention backends on a small fixture before enabling an optimized kernel.

6. **Evaluate static and adaptive attacks.** Use isolated tools, synthetic secrets, and simulated destinations within the authorized test environment. Define attacker control, target behavior, success criteria, and query or compute budgets before running. Evaluate both fixed injections and attacks adapted to the defended model; include transfer attacks from other targets because adaptation can fail to discover stronger attacks. RL-Hammer and Claudini are attack methods named in the source finding; verify their implementations and protocols before use. For white-box evaluation, establish that gradients traverse the actual gated model. Keep training, attack development, and final evaluation tasks separate.

7. **Measure utility and interpret results jointly.** Evaluate the same base checkpoint and defended model with matched prompts, decoding, tool schemas, and budgets. Include zero-shot, few-shot, multi-turn, tool-call validity, and benign extraction or reasoning over external data. Define attack success rate as successful attacks divided by the explicitly stated eligible trials; report counts, exclusions, uncertainty, and whether success is per attempt or cumulative within a budget. Separately report benign task success, task success under attack, refusal rate, latency, and resource usage. Predeclare acceptable utility loss and security targets. Preserve benchmark-native metrics when reproducing published experiments and identify any denominator differences.

### Bounded workflow

- **Minutes 0–3:** Inspect available source material, code, model access, and compute. Choose verified reproduction, clearly labeled adaptation, or planning-only mode.
- **Minutes 3–6:** Produce the trust-channel specification and architecture compatibility checklist. Identify missing prerequisites before providing executable instructions.
- **Minutes 6–11:** If a verified implementation and suitable environment already exist, run a bounded smoke check using its documented interface. Otherwise, produce an implementation specification and evaluation matrix; do not substitute guessed commands.
- **Minutes 11–15:** Report evidence, blockers, smoke-check results if available, and the next bounded experiment. Full model training and adaptive attack campaigns are outside the 15-minute completion claim.

Treat the numbered capabilities as a research protocol proposed by this skill except where explicitly attributed to the paper. Use only verified entry points for executable steps. A valid planning-only result identifies what is missing and makes no claim that training, robustness evaluation, or deployment readiness has been demonstrated.

## Inputs / Outputs

**Inputs**

- Versioned paper, accessible implementation or local candidate code, and model/tokenizer identifiers.
- Agent message format, tool schemas, provenance metadata, and the intended trust policy.
- Authorized training and evaluation data, attack surface, success criteria, and held-out split definitions.
- Hardware and software constraints, execution budget, and acceptable security and utility thresholds.

**Outputs**

- Readiness report with source verification status, architecture support, execution mode, and blockers.
- Token-channel specification with representative trusted, untrusted, and template spans and alignment checks.
- Training plan or verified configuration recording frozen versus trainable parameters and data provenance.
- Evaluation matrix covering baseline, gated model, static attacks, adaptive attacks, and benign utility.
- Evidence ledger distinguishing paper-reported claims, locally measured results, proposed experiments, and untested assumptions.

For each measured result, include model and code revisions, dataset split, sample count, seed, budget, metric definition, and observed failures. Mark unexecuted evaluations as `not_run`, never as passing or zero attack success. Consider the workflow complete when these outputs are internally consistent and prerequisites or blockers are explicit; deployment requires separate evidence that the declared acceptance criteria are met. Retain application-level tool authorization and access controls regardless of measured robustness.

## References

- [Heads of Security: Learned Attention Gates Improve Prompt Injection Robustness, arXiv:2610.37657v1](http://arxiv.org/abs/2610.37657v1) — source finding; empirical claims apply to the reported configurations.
- [Versioned full text](https://arxiv.org/html/2610.37657v1) — method in Section 3, sink treatment in Appendix D, and training details in Appendix J.
- [Paper-linked implementation repository](https://github.com/imperial-aisp/heads_of_security) — returned 404 during verification on 2026-10-28; recheck before relying on code or installation instructions.
