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
name: azure-ai-foundry-operations-2026
description: Implement and operate Azure AI Foundry and Azure OpenAI workloads with enterprise auth, networking, and deployment controls. Use when deploying model endpoints, migrating model versions, or setting up production guardrails on Azure.
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Azure AI Foundry Operations (2026)

## Workflow

1. Confirm the active Azure AI Foundry and Azure OpenAI docs in `references/sources.md`.
2. Choose the deployment path: Foundry model inference, Azure OpenAI, or a mixed Foundry project workflow.
3. Implement Entra ID, service principal, and network isolation requirements before inference code.
4. Add rate-limit handling, tracing, and deployment rollback notes for model updates.
5. Run a focused smoke test against the target endpoint before broader rollout.

## Output Requirements

- State the selected Azure AI service and deployment path.
- State the auth and network boundary approach.
- Include one rollout guardrail for version or endpoint changes.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
