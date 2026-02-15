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

# Digital Biomarker Agent (SKDH)

**ID:** `biomedical.clinical.digital_biomarkers`
**Version:** 1.0.0
**Status:** Beta
**Category:** Clinical / Wearables

---

## Overview

The **Digital Biomarker Agent** utilizes **SciKit Digital Health (SKDH)** to process high-frequency inertial sensor data (accelerometer, gyroscope) from wearables like Apple Watch or Fitbit. It extracts clinically relevant endpoints for movement disorders, sleep, and physical activity.

## Key Capabilities

- **Gait Analysis:** Stride length, gait speed, asymmetry (e.g., for Parkinson's monitoring).
- **Sleep Detection:** Sleep/Wake classification from actigraphy.
- **Tremor Quantification:** Detecting and quantifying resting vs. action tremors.

## Pipeline

1.  **Ingest:** Raw CSV/Parquet data from wearable.
2.  **Preprocess:** Calibration and gravity filtering.
3.  **detect:** Event detection (e.g., "Step", "Sit-to-Stand").
4.  **Extract:** Calculate summary metrics.

## References
- *Pfizer SciKit Digital Health*
- *Digital Biomarker Discovery Pipeline (DBDP)*


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->