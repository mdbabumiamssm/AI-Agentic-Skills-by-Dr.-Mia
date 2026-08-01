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

9. **Comparative evidence screening workflow**
   Explicitly sample non-tumor and low-magnification regions alongside tumor regions; compare tumor-only, non-tumor-ablated, low-magnification-ablated, and full-region strategies; aggregate tile evidence with documented logic; inspect attention or attribution outputs for interpretability; validate performance on external cohorts; and report predictions only as screening signals requiring confirmatory diagnostic testing.

10. **Evidence-aligned dMMR risk prioritization**
    Use multi-region sampling across tumor, non-tumor, and low-magnification tissue; enforce leakage-safe patient-level train, validation, and test splits; compare magnification-level and region ablations; review attention maps and attribution outputs for plausible regional evidence; perform external validation; and frame outputs as dMMR risk prioritization for confirmatory testing rather than a definitive diagnosis.

11. **Leakage-controlled multi-magnification tiling plan**
    Include tumor, non-tumor, and low-magnification WSI tiles when data permit; keep region and magnification strata explicit through sampling, aggregation, error analysis, and reporting; split cohorts at patient or case level before tiling to prevent leakage; validate tumor-only versus multi-region strategies against reference dMMR/MSI labels; and hand off results as screening evidence that must be reconciled with standard clinical testing.

12. **Non-tumor low-magnification validation pattern**
    Use the reported colorectal cancer dMMR histopathology finding as a validation pattern: deliberately sample non-tumor and low-magnification WSI regions alongside tumor regions, require region-level attribution or ablation to show how regions influence slide-level risk, document slide sampling and aggregation logic, route positive or uncertain outputs to clinical confirmatory testing, and flag deployment caveats such as site-specific staining, scanner, preprocessing, and shortcut-learning risks.

13. **Biomarker model interpretation beyond tumor-only tiles**
    Interpret dMMR biomarker models with explicit attention to non-tumor and low-magnification WSI regions reported by the 2026 colorectal cancer histopathology finding; avoid assuming tumor-only tiles are the sole valid evidence source; document region-selection rationale, attribution or ablation support, and validation handoff for confirmatory clinical testing.

14. **Region-level dMMR evidence handoff guardrails**
    Preserve the mapping from non-tumor, tumor, and low-magnification WSI evidence to slide-level risk scores; validate region contributions before relying on the biomarker signal; hand off positive, uncertain, or discordant predictions for standard dMMR/MSI confirmation; and check site, scanner, stain, preprocessing, and tissue-background effects for shortcut learning.

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
