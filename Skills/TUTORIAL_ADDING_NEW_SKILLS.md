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

# Tutorial: Adding New Skills (The SKILL.md Standard)

This project adopts the `SKILL.md` standard to make agent capabilities discoverable, self-documenting, and executable by AI systems.

## What is a Skill?

A **Skill** is a discrete capability that an AI agent can "learn" and execute. It wraps code (Python scripts, APIs) or knowledge (Checklists, Rules) into a standardized format.

## The Standard Format

Every skill folder should contain a `SKILL.md` file.

### Directory Structure
```text
Skills/
└── Domain_Name/
    └── Skill_Name/
        ├── SKILL.md          # The Definition (Required)
        ├── impl.py           # The Implementation (Optional)
        └── references/       # Supporting docs (Optional)
            └── rules.md
```

### SKILL.md Template

```markdown
---
name: my-new-skill
description: A one-sentence description of what this skill does.
license: MIT
metadata:
  author: Your Name
  version: "1.0.0"
compatibility:
  - system: Python 3.9+
allowed-tools:
  - run_shell_command
  - read_file
---

# My New Skill

A more detailed description of the skill.

## When to Use This Skill

*   Trigger condition 1
*   Trigger condition 2

## Core Capabilities

1.  **Capability A**: Description.
2.  **Capability B**: Description.

## Workflow

1.  **Step 1**: What the agent should do first.
2.  **Step 2**: The command to run or file to read.

## Example Usage

**User**: "Do X."

**Agent Action**:
```bash
python3 path/to/script.py --arg "value"
```
```

## Step-by-Step Guide

1.  **Identify the Capability**: Is it a Python script? A set of rules? A complex workflow?
2.  **Create the Folder**: Place it in the appropriate category (e.g., `Clinical`, `Computer_Science`, `External_Collections`).
3.  **Draft SKILL.md**: Use the template above. Ensure the `description` is high-quality as it acts as the "prompt" for the agent to select this skill.
4.  **Implement Code**: If the skill requires execution, write the Python/Shell script and ensure it handles CLI arguments.
5.  **Test**: Verify that an agent reading the `SKILL.md` can successfully execute the workflow.

## Best Practices

*   **Atomic**: Skills should do one thing well.
*   **Discoverable**: Use clear names and descriptions.
*   **Safe**: Define `allowed-tools` to limit scope.
*   **Contextual**: Provide `references/` for complex rules or knowledge bases.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->