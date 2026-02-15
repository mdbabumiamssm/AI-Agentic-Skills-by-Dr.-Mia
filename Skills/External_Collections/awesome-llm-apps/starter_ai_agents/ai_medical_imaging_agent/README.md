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

# 🩻 Medical Imaging Diagnosis Agent

A Medical Imaging Diagnosis Agent build on agno powered by Gemini 2.0 Flash that provides AI-assisted analysis of medical images of various scans. The agent acts as a medical imaging diagnosis expert to analyze various types of medical images and videos, providing detailed diagnostic insights and explanations.

## Features

- **Comprehensive Image Analysis**
  - Image Type Identification (X-ray, MRI, CT scan, ultrasound)
  - Anatomical Region Detection
  - Key Findings and Observations
  - Potential Abnormalities Detection
  - Image Quality Assessment
  - Research and Reference

## How to Run

1. **Setup Environment**
   ```bash
   # Clone the repository
   git clone https://github.com/Shubhamsaboo/awesome-llm-apps.git
   cd ai_agent_tutorials/ai_medical_imaging_agent

   # Install dependencies
   pip install -r requirements.txt
   ```

2. **Configure API Keys**
   - Get Google API key from [Google AI Studio](https://aistudio.google.com)

3. **Run the Application**
   ```bash
   streamlit run ai_medical_imaging.py
   ```

## Analysis Components

- **Image Type and Region**
  - Identifies imaging modality
  - Specifies anatomical region

- **Key Findings**
  - Systematic listing of observations
  - Detailed appearance descriptions
  - Abnormality highlighting

- **Diagnostic Assessment**
  - Potential diagnoses ranking
  - Differential diagnoses
  - Severity assessment

- **Patient-Friendly Explanations**
  - Simplified terminology
  - Detailed first-principles explanations
  - Visual reference points

## Notes

- Uses Gemini 2.0 Flash for analysis
- Requires stable internet connection
- Free API usage costs -  1,500 free requests per day by google!
- For educational and development purposes only
- Not a replacement for professional medical diagnosis

## Disclaimer

This tool is for educational and informational purposes only. All analyses should be reviewed by qualified healthcare professionals. Do not make medical decisions based solely on this analysis.

<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->