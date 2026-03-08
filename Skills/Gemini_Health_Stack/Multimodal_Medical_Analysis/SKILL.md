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
name: gemini-multimodal-medical-analysis
description: Analyzes multimodal medical inputs (X-ray, MRI, Pathology slides + Clinical notes) using Gemini.
keywords:
  - multimodal
  - radiology
  - pathology
  - vision
  - gemini
  - health-stack
measurable_outcome: Generates multimodal diagnostic insights in under 10 seconds.
license: MIT
metadata:
  author: AI Agentic Skills Team
  version: "1.0.0"
compatibility:
  - system: Python 3.10+
allowed-tools:
  - read_file
  - gemini_vision_api
---

# Gemini Multimodal Medical Analysis

A specialized skill leveraging the **Gemini Health Stack** (e.g., Gemini 1.5 Pro) to analyze complex multimodal medical data, such as combining radiological imaging with structured EHR data.

## When to Use This Skill

*   When the user provides both text (clinical history) and an image (X-ray, slide).
*   For generating preliminary differential diagnoses from multimodal context.
*   To extract specific biomarkers from stained pathology slides.

## Example Usage

**User**: "Analyze this chest X-ray (`patient_cxr.png`) along with the patient's EHR notes (`ehr.txt`) and flag any abnormalities."

**Agent Action**:
```bash
python3 Skills/Gemini_Health_Stack/Multimodal_Medical_Analysis/multimodal_agent.py --image "patient_cxr.png" --notes "ehr.txt"
```
