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

# Tutorial: Managing the "Best of Best" Skills Collection

This repository uses a strict standardization protocol to ensure all 600+ skills are discoverable, measurable, and safe for autonomous agents.

## 1. The Standard: `SKILL.md`

Every skill directory must contain a `SKILL.md` file. This file acts as the "API definition" for the skill.

### Required Frontmatter (YAML)
```yaml
---
name: skill-name-slug
description: A 1-sentence summary of what the skill does.
keywords:
  - tag1
  - tag2
measurable_outcome: A SMART goal (Specific, Measurable, Achievable, Relevant, Time-bound).
allowed-tools:
  - tool_name_1
  - tool_name_2
---
```

### Required Body (Markdown)
The body should explain:
- **When to Use:** Context for the agent.
- **Core Capabilities:** What features it offers.
- **Workflow:** Step-by-step execution plan.
- **Example Usage:** CLI command or code snippet.
- **Guardrails:** Safety and accuracy checks.

## 2. Adding a New Skill

1.  **Create Directory:**
    ```bash
    mkdir -p Skills/Domain/My_New_Skill
    ```
2.  **Create `SKILL.md`:**
    Copy the template above or use an existing skill as a reference.
3.  **Validate:**
    Run the catalog generator to check for errors.
    ```bash
    python3 platform/skills_catalog.py
    ```

## 3. Maintaining the Catalog

The catalog is a JSON file (`skills_catalog.json`) that indexes all skills. It is generated automatically.

**To update the catalog:**
```bash
python3 platform/skills_catalog.py
```

This script will:
- Scan all `SKILL.md` files.
- Validate required fields.
- Report any errors or missing metadata.
- Update `skills_catalog.json`.

## 4. Bulk Updates

If you need to update metadata across many skills (e.g., adding a new allowed tool), modify and run:
```bash
python3 platform/fix_skills.py
```
*(Note: Always backup your Skills directory before running bulk updates)*

## 5. Best Practices

- **Measurable Outcomes:** Avoid vague goals like "Analyze data." Use "Generate a PCA plot from 10k cells within 5 minutes."
- **Allowed Tools:** Least Privilege Principle. Only list tools the agent *absolutely needs*.
- **Keywords:** Use standard biomedical terms (e.g., "scRNA-seq", "CRISPR", "Survival Analysis") for better indexing.

---
**Maintained by MD BABU MIA, PhD**


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->