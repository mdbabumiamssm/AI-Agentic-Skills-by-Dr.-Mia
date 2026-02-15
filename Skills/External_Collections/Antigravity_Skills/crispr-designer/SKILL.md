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

---
name: 'crispr-designer'
description: 'Designs guide RNA (gRNA) sequences for CRISPR-Cas9 editing, including off-target analysis. Use when a user needs to edit a gene or asks for gRNA sequences.'
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# CRISPR gRNA Designer

This skill designs high-efficiency guide RNAs for gene editing experiments.

## When to use this skill
- When a user provides a gene name (e.g., "TP53"), Ensembl ID, or DNA sequence.
- When the user wants to perform "knockout", "activation" (CRISPRa), or "interference" (CRISPRi).
- When specificity and off-target minimization are requested.

## How to use it
1.  **Identify Target Locus:**
    -   Resolve the gene name to the current reference genome (GRCh38 for humans).
    -   Identify functional domains (exons) that are constitutively expressed.
2.  **Design gRNAs:**
    -   Select 20nt targets adjacent to NGG PAM sites.
    -   Prioritize 5' constitutive exons for knockouts to ensure early truncation.
3.  **Score Candidates:**
    -   Calculate **On-Target Efficiency** (e.g., Rule Set 2 score).
    -   Calculate **Off-Target Specificity** (CFD score) by searching the whole genome for mismatches.
4.  **Output Table:**
    -   Return a markdown table with: Sequence, PAM, On-Target Score, Off-Target Score, and Genomic Location.
    -   Recommend the top 3 guides.


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->