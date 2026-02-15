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
name: technical-writing-expert
description: Create comprehensive, clear, and structured technical documentation, reports, and whitepapers.
keywords:
  - documentation
  - reporting
  - markdown
  - latex
  - clarity
  - structure
measurable_outcome: Produce a 2000-word technical report with <1% grammatical errors and >90% clarity score (Flesch-Kincaid) within 1 hour.
license: MIT
metadata:
  author: AI Agentic Skills Team
  version: "2.0.0"
compatibility:
  - system: any
allowed-tools:
  - read_file
  - write_file
  - google_web_search
---

# Technical Writing Expert

A specialized skill for generating high-quality technical content, including API documentation, research papers, engineering RFCs, and user manuals.

## When to Use
- **Documentation:** Writing READMEs, API references, or internal wikis.
- **Reporting:** Summarizing complex data analysis or incident reports.
- **Academic:** Drafting sections of research papers (Methods, Results).
- **Architecture:** Creating Request for Comments (RFCs) or design docs.

## Core Capabilities
1.  **Structure Optimization:** Automatically outlines content with logical flow (Introduction -> Context -> Details -> Conclusion).
2.  **Audience Adaptation:** Adjusts tone and complexity for Executive vs. Engineer vs. End-User.
3.  **Format Mastery:** Expert in Markdown, LaTeX, and reStructuredText.
4.  **Visual Integration:** Suggests placement for diagrams/charts (via placeholders).

## Workflow
1.  **Context Gathering:** Analyze input requirements, existing code, or raw data.
2.  **Outlining:** Propose a Table of Contents (ToC) for approval.
3.  **Drafting:** Generate content section-by-section, focusing on clarity and conciseness.
4.  **Review:** Self-critique for passive voice, jargon overuse, and formatting consistency.
5.  **Final Polish:** Add metadata, cross-links, and formatting artifacts.

## Example Usage
```bash
# This skill is typically invoked by an agent with a prompt like:
"Draft a technical README for the new authentication module based on src/auth/README_draft.md"
```

## Guardrails
- **Accuracy:** Always verify technical claims against provided source code/data.
- **Plagiarism:** Do not copy-paste large blocks of external text without attribution.
- **Tone:** Maintain an objective, professional tone (avoid marketing fluff).


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->