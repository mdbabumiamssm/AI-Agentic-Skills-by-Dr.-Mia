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
name: deepseek-api-operations-2026
description: Integrate and operate DeepSeek APIs using current docs and compatibility guidance. Use when implementing DeepSeek-powered chat or reasoning workloads, selecting models, or validating OpenAI-compatible client behavior.
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# DeepSeek API Operations (2026)

## Workflow

1. Verify current model, rate-limit, and release status from `references/sources.md`.
2. Choose whether to integrate through the OpenAI-compatible path or a provider-specific feature path.
3. Align model choice with thinking versus non-thinking behavior before coding.
4. Define timeouts, retries, and compatibility checks for the chosen SDK path.
5. Run a focused smoke test with the exact `base_url` and model name used in production.

## Output Requirements

- State the selected DeepSeek model and endpoint path.
- State the SDK compatibility approach being used.
- Include one operational guardrail for retries, fallback, or release risk.
- For on-premises clinical diagnosis use cases, flag the model-size/diagnostic-accuracy tradeoff: distilled DeepSeek-R1 open-source variants underperformed expectations vs. larger closed models in a comparative study (Zhong et al., J Med Syst, May 2026). Recommend rigorous local evaluation before any clinical deployment.

## References

- Zhong W, Fu Y, Peng D, Liu Y, Liu Y. Open-Source Large Language Models Distilled DeepSeek-R1 Pose Challenges for On-Premises Clinical Deployment in Medical Diagnosis: A Comparative Study of Performance. J Med Syst, 2026 May 1. https://pubmed.ncbi.nlm.nih.gov/42062641/


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
