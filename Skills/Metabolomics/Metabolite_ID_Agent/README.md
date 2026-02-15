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

# Metabolite ID Agent

## Overview
The **Metabolite ID Agent** automates the identification of unknown metabolites from Mass Spectrometry (LC-MS/MS) data. It integrates spectral similarity scoring, fragmentation tree analysis, and machine learning to annotate features with high confidence.

## Core Capabilities
- **Spectral Matching**: Queries local and remote databases (GNPS, MassBank, HMDB).
- **In Silico Fragmentation**: Uses CSI:FingerID and SIRIUS-like logic to predict structures for unknowns.
- **Spec2Vec Integration**: Utilizes spectral embeddings for finding structurally related compounds.
- **Retention Time Prediction**: Validates candidates using deep learning-based RT prediction.

## Workflow
1.  **Input**: .mzML raw data or feature lists (mass/charge, retention time).
2.  **Preprocessing**: Peak picking and alignment (via mzmine or pyopenms).
3.  **Annotation**: AI-driven spectral matching and substructure prediction.
4.  **Output**: Annotated metabolite list with confidence scores.

## Dependencies
- `pyopenms`
- `matchms`
- `spec2vec`
- `pubchempy`


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->