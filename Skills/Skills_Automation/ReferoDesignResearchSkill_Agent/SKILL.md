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



<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->

---
name: 'refero-design-research-skill'
description: 'Research-first UI design skill for using Refero MCP to study real app screens and flows before making product design decisions.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Refero Design Research Skill

## Overview

This skill guides AI agents through research-first UI design using Refero, a design corpus and MCP workflow for studying real app screens and flows. Use it to ground product interface decisions in observed patterns before proposing layouts, interactions, information architecture, or visual treatments.

It is intended for evidence gathering and design synthesis, not for copying proprietary screens or treating examples as universal rules.

## When to Use This Skill

- Designing a new product screen, feature flow, onboarding path, dashboard, editor, checkout, settings area, or mobile/web app workflow.
- Redesigning an existing UI where comparable real-world app screens would clarify conventions and tradeoffs.
- Reviewing UI/UX work and needing concrete product-screen evidence instead of generic design advice.
- Comparing how mature apps handle navigation, forms, empty states, permissions, collaboration, search, filters, or dense operational views.
- Preparing a design brief, implementation prompt, or critique that should cite observed interaction and layout patterns.
- Resolving uncertainty about platform conventions, screen density, component hierarchy, or expected flow structure.

## Core Capabilities

1. **Research framing** - Translate the design request into product category, user role, platform, feature type, workflow stage, and pattern questions.
2. **Refero corpus lookup** - Use the available Refero MCP or related access path to locate relevant real app screens and flows from the Refero design corpus.
3. **Pattern extraction** - Identify repeated approaches across examples, including layout structure, navigation model, component hierarchy, state handling, content density, calls to action, and affordances.
4. **Evidence capture** - Record source app names, screen or flow identifiers, links when available, query terms, and concise observations that justify each design conclusion.
5. **Design synthesis** - Convert observed patterns into practical UI recommendations, annotated design requirements, implementation guidance, or critique notes.
6. **Gap handling** - State when Refero access is unavailable, search results are thin, or a conclusion is based on inference rather than direct screen evidence.

## Inputs / Outputs

**Inputs**

- Product or feature description, including target users and primary jobs to be done.
- Platform and context, such as web app, mobile app, desktop app, SaaS dashboard, consumer app, admin tool, or internal workflow.
- Target screen or flow, such as signup, search, project setup, checkout, settings, analytics, editor, collaboration, or approval workflow.
- Constraints, including brand direction, accessibility needs, responsive behavior, technical stack, density requirements, or existing design system.
- Optional comparison set, competitor category, inspiration apps, or specific UI questions to investigate.

**Outputs**

- A concise research brief listing the queries or research angles used.
- Relevant Refero examples with app names, screen or flow references, and links or identifiers when available.
- Observed patterns, tradeoffs, and notable divergences across examples.
- Design implications for layout, navigation, components, states, copy hierarchy, and interaction behavior.
- Concrete recommendations that can be used in a design spec, implementation prompt, or review.
- Open questions and confidence notes where evidence is limited.

## References

- Source finding: https://github.com/referodesign/refero_skill
- Model Context Protocol GitHub organization: https://github.com/modelcontextprotocol
