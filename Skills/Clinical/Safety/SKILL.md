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
name: 'safety-monitor'
description: 'Safety Monitor'
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---


# Safety Monitor Agent

Agent responsible for compliance and safety checks. It inspects content for safety violations, toxins, or medical misinformation.

## When to Use This Skill

*   Before returning any medical advice to a user.
*   To audit generated molecules for toxicity.
*   To check for Protected Health Information (PHI) leaks.

## Core Capabilities

1.  **Content Inspection**: Scans text for dangerous content.
2.  **LLM-based Audit**: Uses a separate LLM call to verify safety.
3.  **Status Reporting**: Returns 'approved', 'flagged', or 'rejected'.

## Workflow

1.  **Input**: Text content to check.
2.  **Execute**: Run the safety agent.
3.  **Output**: JSON status report.

## Example Usage

**User**: "Check this drug proposal for safety."

**Agent Action**:
```bash
python3 Skills/Clinical/Safety/safety_agent.py
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->