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
name: 'dmmr-crc-histopathology-agent'
description: 'Predict dMMR risk in colorectal cancer histopathology workflows using tumor, non-tumor, and low-magnification WSI regions with validation handoff.'
measurable_outcome: 'Execute skill workflow successfully with valid output within 15 minutes.'
allowed-tools:
  - read_file
  - run_shell_command
  - web_fetch
---

# dMMR CRC Histopathology Agent

## Overview

This skill guides computational pathology workflows for predicting mismatch repair deficiency (dMMR) from colorectal cancer histopathology whole-slide images. It emphasizes region-aware review, including tumor, non-tumor, and low-magnification regions, because the source finding reports that these areas can contribute to dMMR prediction. Use it to structure analysis, validation, reporting, and clinical-pathology handoff without treating model output as a standalone diagnosis.

## When to Use This Skill

- Assessing colorectal cancer histopathology slides for dMMR or MSI-related risk signals.
- Designing or reviewing a WSI pipeline that includes tumor, non-tumor, stromal, mucosal, or low-magnification context.
- Comparing region-selection strategies for biomarker prediction in colorectal cancer.
- Preparing a validation checklist before deploying a dMMR histopathology model.
- Translating computational pathology outputs into pathology-facing reports that can be reconciled with IHC, PCR, or NGS results.

## Core Capabilities

1. **Case and slide intake**
   Confirm diagnosis context, specimen type, slide identifiers, stain type, scanner metadata, resolution, and available ground-truth dMMR/MSI labels.

2. **Quality control**
   Screen for missing tissue, blur, folds, pen marks, staining artifacts, scanner artifacts, tissue fragmentation, and label mismatches before inference or model review.

3. **Region-aware tissue parsing**
   Separate tumor-rich regions from non-tumor tissue, invasive front, stroma, lymphoid aggregates, normal mucosa, necrosis, and low-magnification context when annotations or segmentation are available.

4. **Low-magnification feature review**
   Preserve slide-level architecture and broad tissue context during analysis instead of restricting review to only high-power tumor tiles.

5. **Multi-region sampling and attribution audit**
   Sample tumor, non-tumor, and low-magnification regions explicitly; aggregate their predictions at slide level; run region-ablation comparisons and attention audits to assess each region's contribution; validate findings on external cohorts; and check scanner, staining, fixation, sectioning, and other tissue-processing variables to guard against shortcut learning.

6. **Inference workflow support**
   Run or audit a dMMR prediction pipeline by documenting preprocessing, magnification levels, tile sampling, aggregation logic, calibration method, and model version.

7. **Validation and error analysis**
   Compare predictions against reference testing, stratify errors by tissue region and slide quality, and flag cases where non-tumor or low-magnification signals may drive model behavior.

8. **Clinical-pathology handoff**
   Produce a concise report with case identifiers, model inputs, region basis, limitations, confidence or risk category if supplied by the model, and recommended correlation with standard molecular or immunohistochemical testing.

## Inputs / Outputs

### Inputs

- Colorectal cancer whole-slide images or tile exports, preferably H&E-stained.
- Slide metadata including scanner, objective magnification, microns-per-pixel, stain, block, and case identifiers.
- Optional tumor annotations, tissue masks, region labels, or segmentation outputs.
- Optional reference labels from mismatch repair immunohistochemistry, MSI PCR, NGS, or curated clinical records.
- Optional trained model artifacts, inference scripts, calibration tables, and validation cohort definitions.

### Outputs

- A slide- or case-level dMMR risk summary suitable for computational pathology review.
- Region-level notes identifying whether tumor, non-tumor, or low-magnification areas contributed to interpretation.
- Quality-control flags and exclusions that may affect reliability.
- Validation notes comparing model outputs with available reference testing.
- A pathology-facing handoff statement that clearly separates computational prediction from diagnostic confirmation.

## References

- PubMed source finding: https://pubmed.ncbi.nlm.nih.gov/41875848/
