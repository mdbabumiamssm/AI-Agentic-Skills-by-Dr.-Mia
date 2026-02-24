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
name: 'regulatory-drafting'
description: 'High-speed Regulatory Drafting agent for healthcare workflows (FDA/EMA submissions).'
measurable_outcome: Draft a complete regulatory submission (e.g., IND, CTD Module 2) in < 15 minutes (Benchmark: 10 mins).
model_requirements: 'Claude 3.5 Sonnet (Speed) or Opus (Reasoning)'
allowed-tools:
  - read_file
  - run_shell_command
  - mcp_bio_tool (via BioMCP)
---


# Regulatory Drafting Skill (2026 Edition)

This skill implements the **Regulatory Drafting** workflow using Anthropic's Claude for Healthcare stack. It dramatically accelerates the creation of regulatory documents by automating literature review, data extraction, and writing.

## 🚀 Capabilities

*   **Speed:** Reduces drafting time from **12 weeks to ~10 minutes**.
*   **Precision:** Uses citation-backed generation to ensure factual accuracy.
*   **Compliance:** Adheres to FDA/EMA formatting guidelines.

## 🛠️ Usage

```bash
# Run the Regulatory Coworker (Interactive CLI)
python3 Skills/Anthropic_Health_Stack/Regulatory_Drafting/coworker.py --mode draft --topic "CTD Module 2.4 Nonclinical Overview"
```

## 🔗 Integration

This skill can be integrated into larger workflows via the **Model Context Protocol (MCP)**.
See `MCP_Servers/BioMCP/bio_mcp_server.py` for details.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
