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

# Radiology Agent

**ID:** `biomedical.clinical.radiology`
**Version:** 0.5.0
**Status:** Experimental
**Category:** Clinical / Imaging

---

## Overview

The **Radiology Agent** utilizes multimodal LLMs (Vision-Language Models) to assist in the interpretation of medical images (X-ray, CT, MRI). It is designed to act as a "second reader," generating preliminary reports and answering visual questions.

## Key Capabilities

- **Report Generation:** Drafts structured radiology reports describing findings, impressions, and recommendations.
- **Visual QA:** Answers specific questions like "Is there a pleural effusion?" or "Locate the mass in the left lung."
- **Modality Support:**
  - Chest X-rays (CXR)
  - Brain MRI (Tumor segmentation/identification)
  - CT Scans (Abdominal/Thoracic)

## Models & Tools

- **RadGPT:** Specializes in conversational radiology report generation.
- **VILA-M3:** Advanced multimodal model for medical imaging.
- **MONAI:** Underlying framework for image transformations and deep learning tasks.

## References
- *RadGPT* (Stanford, 2025)
- *MIMIC-CXR Benchmark*


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->