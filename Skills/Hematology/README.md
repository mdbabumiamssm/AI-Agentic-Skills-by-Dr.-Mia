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

# Hematology Skills

**Category:** Hematology / Blood Disorders / AI Diagnostics
**Version:** 1.0.0

---

## Overview

The **Hematology Skills** collection provides AI-powered diagnostic and treatment planning capabilities for blood disorders including leukemias, lymphomas, myelodysplastic syndromes, plasma cell neoplasms, and clonal hematopoiesis. These skills integrate cutting-edge machine learning with clinical hematopathology to deliver faster, more accurate diagnoses and personalized treatment recommendations.

---

## Available Skills

### 1. Flow Cytometry AI Agent
- **Path:** `Flow_Cytometry_AI/`
- **Purpose:** Automated flow cytometry analysis using deep learning
- **Key Features:**
  - DeepFlow integration (95% accuracy, 100x faster)
  - Acute leukemia classification
  - MRD detection
  - MDS flow scoring (Ogata score)

### 2. Clonal Hematopoiesis Agent
- **Path:** `Clonal_Hematopoiesis/`
- **Purpose:** CHIP/CCUS identification and risk assessment
- **Key Features:**
  - CHIP vs CCUS classification
  - MDS/AML progression risk
  - Cardiovascular risk assessment
  - Surveillance planning

### 3. AML Classification Agent
- **Path:** `AML_Classification/`
- **Purpose:** WHO 2022 AML classification and ELN risk stratification
- **Key Features:**
  - WHO 2022 subtyping
  - ELN 2022 risk stratification
  - Treatment recommendations
  - MRD monitoring design

### 4. MDS Diagnosis Agent
- **Path:** `MDS_Diagnosis/`
- **Purpose:** MDS diagnosis and IPSS-M risk scoring
- **Key Features:**
  - WHO 2022 MDS classification
  - IPSS-M calculation
  - Treatment selection
  - Progression monitoring

### 5. Multiple Myeloma AI Agent
- **Path:** `Multiple_Myeloma_AI/`
- **Purpose:** Myeloma diagnosis, staging, and treatment planning
- **Key Features:**
  - MGUS/SMM/MM differentiation
  - R-ISS staging
  - Treatment sequencing
  - MRD monitoring

### 6. Blood Smear Analysis Agent
- **Path:** `Blood_Smear_Analysis/`
- **Purpose:** AI-powered peripheral blood smear analysis
- **Key Features:**
  - Automated cell classification
  - Morphologic abnormality detection
  - FDA-approved system integration

### 7. Lymphoma AI Agent
- **Path:** `Lymphoma_AI/`
- **Purpose:** Lymphoma classification and treatment
- **Key Features:**
  - WHO classification support
  - Risk stratification
  - CAR-T eligibility assessment

---

## Integrated AI Systems

| System | FDA Status | Application | Performance |
|--------|------------|-------------|-------------|
| **DeepFlow** | Research | Flow cytometry | 95% accuracy |
| **CellaVision** | FDA Approved | Blood smear | Automated counting |
| **Scopio X100** | FDA Approved | Digital morphology | High resolution |
| **ALMA** | Research | Leukemia methylome | 27 WHO subtypes |

---

## Clinical Workflows

### Acute Leukemia Workup
```
Presentation (cytopenias, blasts)
    ↓
Peripheral Blood Smear (Blood_Smear_Analysis)
    ↓
Flow Cytometry (Flow_Cytometry_AI)
    ↓
Lineage Determination
    ├── Myeloid → AML_Classification
    └── Lymphoid → Lymphoma_AI
    ↓
Cytogenetics + Molecular
    ↓
WHO Classification + Risk Stratification
    ↓
Treatment + MRD Monitoring
```

### Cytopenia Evaluation
```
Unexplained Cytopenias
    ↓
Peripheral Blood + Bone Marrow
    ↓
Flow Cytometry (MDS-flow score)
    ↓
NGS Panel
    ├── Mutations detected → Clonal_Hematopoiesis
    │   ├── CHIP (no cytopenias)
    │   ├── CCUS (cytopenias, no dysplasia)
    │   └── MDS (dysplasia) → MDS_Diagnosis
    └── No mutations → Other causes
```

---

## Dependencies

### Required Data Sources
- Complete blood count (CBC)
- Peripheral blood smear
- Bone marrow aspirate/biopsy
- Flow cytometry (FCS files)
- Cytogenetics (karyotype, FISH)
- Molecular panel (NGS)

### Required Tools
- FlowJo or equivalent FCS analysis
- MiXCR/IgBLAST for receptor sequencing
- Standard hematopathology databases

---

## References

- WHO Classification of Tumours of Haematopoietic and Lymphoid Tissues, 5th Edition (2022)
- ELN 2022 AML Guidelines
- IPSS-M Molecular Scoring System
- NCCN Guidelines for Hematologic Malignancies

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->