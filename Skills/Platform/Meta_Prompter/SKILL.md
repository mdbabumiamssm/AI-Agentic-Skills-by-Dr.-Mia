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
name: meta-prompter
description: Automatic prompt engineering & optimization
keywords:
  - prompt-engineering
  - optimization
  - meta-prompting
  - llm
  - tuning
measurable_outcome: Improves prompt performance metrics by >15% over baseline.
license: MIT
metadata:
  author: Biomedical OS Team
  version: "1.0.0"
compatibility:
  - system: Python 3.10+
allowed-tools:
  - run_shell_command
  - read_file
---

# Meta-Prompter

The Meta-Prompter is a tool for self-optimizing agent prompts. It analyzes agent performance and iteratively refines system prompts to maximize accuracy and adherence to instructions.

## When to Use This Skill

*   When an agent is consistently failing a specific type of task.
*   When deploying a new agent and needing to tune its persona.
*   When A/B testing different prompting strategies.

## Core Capabilities

1.  **Prompt Optimization**: Rewriting prompts for clarity and effectiveness.
2.  **Performance Evaluation**: Testing prompts against benchmarks.
3.  **Few-Shot Generation**: Creating optimal examples for context.

## Example Usage

**User**: "Optimize the Clinical Reasoning prompt."

**Agent Action**:
```bash
python3 platform/optimizer/meta_prompter.py --target "clinical_reasoning" --iterations 5
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->