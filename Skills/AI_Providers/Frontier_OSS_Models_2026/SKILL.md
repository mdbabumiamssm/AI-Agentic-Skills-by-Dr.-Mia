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

---
name: frontier-oss-models-2026
description: Work with leading open or openly available model ecosystems from Meta, Mistral, DeepSeek, and xAI. Use when evaluating self-hosted or hybrid deployment paths, model tradeoffs, and integration patterns across non-single-vendor stacks.
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Frontier OSS Models (2026)

## Core Capabilities

- Clinical deployment caution for distilled DeepSeek-R1 and open-source reasoning models: evaluate against medical diagnosis tasks before any clinical use, weigh on-premises deployment tradeoffs, require calibration and safety checks, and do not treat model distillation alone as sufficient evidence of clinical readiness.
- Distilled DeepSeek-R1/open-source medical diagnosis models need on-prem validation, calibration against local clinical tasks, local governance approval, and failure-mode review before self-hosted clinical deployment.
- On-premises clinical diagnosis deployments require diagnosis-specific benchmarking, clinical safety gates, privacy controls, and model provenance checks; PubMed 42062641 reports performance challenges for distilled DeepSeek-R1-style open-source LLMs in medical diagnosis workflows.
- Distilled open-source reasoning models for on-prem clinical use require explicit comparison against proprietary diagnostic baselines, documented assessment of diagnostic benchmark degradation, privacy/local-hosting tradeoff review, clinical guardrails, and post-deployment monitoring before use in patient-facing workflows.
- Clinical on-premises deployment caution for distilled DeepSeek-R1 and other open models: validate against local diagnostic tasks, set calibration and abstention thresholds, enforce privacy and IAM constraints, document hardware and latency tradeoffs, and red-team against commercial models before clinical use.
- Clinical workload routing for distilled DeepSeek-R1 and open-source reasoning models: if on-prem performance validation, diagnostic safety benchmarking, or privacy governance is incomplete, route clinical work to governed provider workflows or human-reviewed pathways instead of raw local model outputs.
- Clinical self-hosting caution for distilled DeepSeek-R1 and open-weight medical diagnosis deployments: benchmark against licensed frontier and domain models, require local evaluation sets, calibration and abstention checks, privacy/on-prem controls, and clinical safety gates before use in diagnostic workflows.
- Clinical on-prem deployment caution for distilled DeepSeek-R1/open-source medical diagnosis models: require benchmark validation on local diagnostic cases, document privacy/security tradeoffs of keeping PHI local versus securing local infrastructure, compare model size against measured diagnostic accuracy, and enforce guardrails that prevent distilled reasoning models from being treated as clinically interchangeable with frontier hosted systems.

## Clinical On-Prem Deployment: Distilled Reasoning Models

- Treat local hosting as a privacy and data-residency advantage, not as evidence of diagnostic readiness.
- Before clinical use, validate distilled open-source reasoning models against proprietary model baselines and local diagnostic tasks; document any diagnostic benchmark degradation without inventing unsupported benchmark names or thresholds.
- Define calibration and abstention thresholds before deployment, including conditions that route cases to human review instead of returning a diagnosis.
- Review privacy, identity and access management, hardware capacity, and latency constraints as part of the on-premises deployment decision.
- Route clinical workloads to governed provider workflows or human-reviewed pathways when raw local model outputs have not passed on-prem diagnostic performance validation, safety benchmarking, and privacy governance review.
- Red-team distilled DeepSeek-R1 and other open models against commercial model baselines before using them in clinical workflows.
- Require clinical safety guardrails, escalation paths, human review, and monitoring for drift, failure modes, and inconsistent diagnostic behavior.
- Record local governance sign-off and known diagnostic failure modes before enabling self-hosted use.

## Workflow

1. Identify target deployment mode: hosted API, self-hosted, or hybrid.
2. Review provider-specific docs and repos in `references/sources.md`.
3. Compare model choice by licensing, latency, and hardware profile.
4. Define eval set before changing production model baselines.
5. Roll out through canary traffic with measurable success thresholds.

## Output Requirements

- State selected model family and deployment mode.
- List one licensing/compliance check.
- List one rollback condition tied to quality or cost.

## References

- https://pubmed.ncbi.nlm.nih.gov/42062641/


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
