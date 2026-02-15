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
name: paper-banana
description: Agentic framework for automating the generation of publication-ready academic illustrations and statistical plots.
license: CC-BY-SA-4.0
metadata:
  author: Peking University & Google Cloud AI Research
  version: "1.0.0"
compatibility:
  - system: Python 3.9+
allowed-tools:
  - run_shell_command
  - read_file
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
---

# PaperBanana

PaperBanana is an advanced agentic framework designed to automate the creation of high-quality, publication-ready academic illustrations. It employs a multi-agent architecture to retrieve data, plan visualizations, style figures, and critique the output, ensuring adherence to strict academic standards.

## When to Use This Skill
- You need to generate methodology diagrams from text descriptions.
- You want to create statistical plots (e.g., bar charts, line graphs, scatter plots) that meet academic publication standards.
- You need to refine existing figures for better clarity, aesthetics, or faithfulness to the data.

## Core Capabilities
1.  **Multi-Agent Orchestration**: Coordinates specialized agents (Retriever, Planner, Stylist, Visualizer, Critic) to handle complex illustration tasks.
2.  **Methodology Diagrams**: Generates flowcharts and system architecture diagrams.
3.  **Statistical Plots**: Produces high-quality plots for data visualization.
4.  **Iterative Refinement**: Uses a critic agent to review and improve figures based on academic criteria.

## Workflow
1.  **Input**: Provide a description of the figure or data to be visualized.
2.  **Planning**: The Planner agent breaks down the request into actionable steps.
3.  **Generation**: The Visualizer and Stylist agents create the initial draft.
4.  **Critique & Refine**: The Critic agent reviews the output, and the system iteratively improves it.
5.  **Output**: A high-resolution image file ready for inclusion in a manuscript.

## References
- [Project Website](https://paperbanana.github.io/) (Placeholder based on typical project structure)


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->