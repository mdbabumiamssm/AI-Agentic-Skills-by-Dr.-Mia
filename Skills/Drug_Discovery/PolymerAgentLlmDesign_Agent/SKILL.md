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
name: 'polymer-agent-llm-design'
description: 'LLM-driven agent for polymer design that proposes, evaluates, and iterates candidate polymer structures against target property constraints.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# Polymer-Agent LLM Design Agent

## Overview
Polymer-Agent operationalizes the JCIM 2026 method by Nigam, Chandrasekhar, and Barati Farimani that uses a large language model as the reasoning backbone of an autonomous agent for polymer design. This skill guides the construction of an iterative loop where the LLM proposes monomer/repeat-unit SMILES, invokes property predictors and chemistry tools, and refines candidates against user-specified targets (e.g., glass transition temperature, density, refractive index, mechanical or thermal properties).

## When to Use This Skill
- Designing novel polymers or copolymers conditioned on target physical/chemical properties.
- Exploring monomer libraries with natural-language constraints (e.g., "transparent, high Tg, biodegradable backbone").
- Coupling LLM reasoning with cheminformatics/property prediction tools for materials discovery.
- Generating rationale-rich design traces for inverse polymer design beyond small-molecule agents.
- Benchmarking agentic chemistry workflows on polymer-specific tasks.

## Core Capabilities
1. **Goal parsing** — Translate a natural-language design brief into structured property targets, ranges, and hard constraints (chemistry classes to include/exclude).
2. **Candidate proposal** — Use the LLM to emit polymer repeat units in SMILES (with `[*]` connection points) or BigSMILES, optionally biased by a seed monomer library.
3. **Tool-augmented evaluation** — Call property predictors (e.g., RDKit descriptors, group-contribution methods, or learned models) on each candidate and return numeric scores plus failure modes.
4. **Reflective refinement** — Feed evaluation results back into the LLM so it critiques deviations from targets and proposes structural edits (substitutions, backbone changes, copolymer ratios).
5. **Stop / budget control** — Track iteration count, score history, and convergence; halt when targets are met or a budget is exhausted.
6. **Provenance logging** — Persist each round (prompt, candidate SMILES, tool outputs, LLM critique) for reproducibility and downstream lab handoff.

## Inputs / Outputs
**Inputs**
- Natural-language design brief and/or structured target dict (property → target value or range).
- Optional seed monomer set, allowed/forbidden functional groups, and synthesis-accessibility constraints.
- Tool configuration: predictor endpoints, RDKit availability, LLM model id and credentials.
- Iteration budget and stopping thresholds.

**Outputs**
- Ranked list of candidate polymer repeat units (SMILES/BigSMILES) with predicted properties.
- Per-candidate rationale and edit history showing how the agent reached the design.
- Structured run log (JSON) of prompts, tool calls, and scores for each iteration.
- Optional summary report contrasting top candidates against the original target brief.

## References
- Source paper: Nigam V, Chandrasekhar A, Barati Farimani A. *Polymer-Agent: Large Language Model Agent for Polymer Design.* J Chem Inf Model, 2026 Apr 28. https://pubmed.ncbi.nlm.nih.gov/42048526/
- BigSMILES specification for stochastic polymer representation: Lin et al., *ACS Cent. Sci.* 2019. https://pubs.acs.org/doi/10.1021/acscentsci.9b00476
- RDKit cheminformatics toolkit: https://www.rdkit.org/
- Related agentic chemistry frameworks: ChemCrow (Bran et al., 2023, https://arxiv.org/abs/2304.05376) and Coscientist (Boiko et al., *Nature* 2023).
