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
name: 'plan-and-solve-agent'
description: 'Breaks down complex queries into a step-by-step plan before execution, improving performance on multi-hop reasoning tasks.'
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---


# Plan-and-Solve Agent

The **Plan-and-Solve Agent** separates high-level planning from low-level execution. It is ideal for complex scientific inquiries that require multiple distinct steps (e.g., "Find targets for disease X, then design drugs, then check safety").

## When to Use This Skill

*   When a user query is too complex for a single "ReAct" loop.
*   When you need to visualize the reasoning process *before* committing to execution.
*   To orchestrate multiple specialized sub-agents.

## Core Capabilities

1.  **Decomposition**: Splits a goal into linear or parallel sub-tasks.
2.  **Execution**: runs each step sequentially (mocked in this version).
3.  **Reporting**: Summarizes the outputs of all steps.

## Workflow

1.  **Input**: A complex natural language query.
2.  **Plan**: The agent generates a list of `PlanNode` objects.
3.  **Execute**: The agent iterates through nodes, executing them (simulation).

## Example Usage

**User**: "Investigate the impact of variant X on drug response."

**Agent Action**:
```bash
python3 Skills/Agentic_AI/Agent_Architectures/Plan_and_Solve/plan_and_solve.py \
    --query "Investigate the impact of variant X on drug response."
```

```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->