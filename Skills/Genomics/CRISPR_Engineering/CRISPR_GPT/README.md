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

# CRISPR-GPT Agent

## Overview
**CRISPR-GPT Agent** is an autonomous system for designing gene-editing experiments. Inspired by the 2025 Stanford tool, it orchestrates the entire workflow from target selection to validation primer design, acting as a "digital lab assistant" for genome engineering.

## Modules
1.  **Guide Design**: Selects optimal sgRNAs using off-target prediction (e.g., DeepCRISPR, CCLMoff).
2.  **Donor Design**: Constructs HDR templates for knock-ins or base editing windows.
3.  **Primer Design**: Generates PCR primers for genotyping and Sanger sequencing validation.
4.  **Protocol Generation**: Writes a step-by-step wet-lab protocol for the experiment.

## Interaction
- **Chat Interface**: "Design a knockout for TP53 in HEK293T cells."
- **Output**: A comprehensive report including oligonucleotide sequences and experimental conditions.

## Reference
- *CRISPR-GPT: An LLM Agent for Automated Gene Editing Design (Nature BME 2025)*


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->