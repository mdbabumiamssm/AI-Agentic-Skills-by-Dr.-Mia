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
name: 'multimodal-medical-imaging'
description: 'Analyzes medical images (X-ray, MRI, CT) using multimodal LLMs to identify anomalies and generate reports.'
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---


# Multimodal Medical Imaging Analysis

The **Multimodal Medical Imaging Analysis Skill** leverages state-of-the-art Vision-Language Models (VLMs) like Gemini 1.5 Pro and GPT-4o to interpret medical imagery alongside clinical text.

## When to Use This Skill

*   When you need a preliminary screening of medical images.
*   When correlating visual findings with textual clinical notes.
*   To generate structured reports (DICOM-SR-like) from raw images.

## Core Capabilities

1.  **Anomaly Detection**: Identify potential pathologies in X-rays, CTs, etc.
2.  **Report Generation**: Draft radiology reports in standard formats.
3.  **VQA (Visual Question Answering)**: Answer specific questions about an image (e.g., "Is there a fracture in the left femur?").
4.  **Dermatology and Dermoscopy Evaluation Profile**: For suspected basal cell carcinoma and predefined common mimicker sets, evaluate paired clinical and dermoscopic images when available; stratify image quality for focus, lighting, scale, lesion completeness, and artifacts or occlusion before interpretation; record available lesion metadata; frame outputs as triage support or a ranked lesion-level differential diagnosis with calibrated confidence, not as a standalone definitive diagnosis; abstain when image quality or confidence is insufficient; report clinician-confirmed evaluation metrics such as sensitivity and specificity for evaluation datasets; analyze lesion-level errors, including confusion among mimickers; and require dermatologist review before clinical use.

## Workflow

1.  **Input**: Provide an image file path (JPG, PNG) and a specific clinical question or "generate report" instruction.
2.  **Analyze**: The agent sends the image and prompt to the VLM.
3.  **Output**: Returns a JSON object with findings, confidence scores, and reasoning.

## Example Usage

**User**: "Analyze this chest X-ray for pneumonia."

**Agent Action**:
```bash
python3 Skills/Clinical/Medical_Imaging/Multimodal_Analysis/multimodal_agent.py \
    --image "/path/to/cxr.jpg" \
    --prompt "Check for signs of pneumonia and consolidation."
```



## References

- https://pubmed.ncbi.nlm.nih.gov/41952838/

<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
