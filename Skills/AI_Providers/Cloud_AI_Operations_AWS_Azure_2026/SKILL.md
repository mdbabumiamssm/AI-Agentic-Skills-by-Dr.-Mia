<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

---
name: cloud-ai-operations-aws-azure-2026
description: Operate AI workloads on AWS Bedrock and Azure AI/Azure OpenAI with production-focused cloud controls. Use when selecting managed model providers, implementing enterprise auth, and designing resilient cloud-native inference pipelines.
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Cloud AI Operations: AWS + Azure (2026)

## Workflow

1. Choose cloud target (AWS, Azure, or both) and capture compliance constraints.
2. Verify current service docs and SDK references in `references/sources.md`.
3. Implement auth first (IAM/STS or Entra/service principal).
4. Add observability hooks before scaling traffic.
5. Validate with low-volume staged inference tests.

## Output Requirements

- Name selected cloud AI service and reason.
- Specify auth pattern and secret-handling approach.
- Include one failover strategy across regions or providers.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->