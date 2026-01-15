# Multiple Myeloma AI Agent

**ID:** `biomedical.hematology.multiple_myeloma_ai`
**Version:** 1.0.0
**Status:** Production
**Category:** Hematology / Multiple Myeloma / Diagnosis & Treatment

---

## Overview

The **Multiple Myeloma AI Agent** provides comprehensive diagnostic evaluation, risk stratification, and treatment planning for multiple myeloma (MM) and related plasma cell neoplasms. Multiple myeloma is characterized by clonal plasma cell proliferation, monoclonal protein production, and end-organ damage (CRAB criteria).

This agent integrates deep learning models for bone marrow assessment, survival prediction algorithms, and MRD detection strategies to deliver precision medicine approaches for myeloma management.

---

## Key Capabilities

### 1. Diagnostic Classification

| Entity | Criteria | Risk of Progression |
|--------|----------|---------------------|
| **MGUS** | M-protein <3g/dL, <10% PCs, no CRAB | 1%/year to MM |
| **Smoldering MM** | M-protein ≥3g/dL or ≥10% PCs, no CRAB | Variable (20-2-20) |
| **Multiple Myeloma** | ≥10% PCs + CRAB or SLiM criteria | Active disease |
| **Plasma cell leukemia** | >2x10^9/L circulating PCs | Aggressive |

### 2. CRAB/SLiM Criteria

| Criterion | Definition |
|-----------|------------|
| **C** - Calcium | >11 mg/dL or >1 mg/dL above ULN |
| **R** - Renal | Creatinine >2 mg/dL or CrCl <40 |
| **A** - Anemia | Hgb <10 g/dL or >2 below LLN |
| **B** - Bone | ≥1 lytic lesion on imaging |
| **S** - Sixty | ≥60% clonal plasma cells |
| **Li** - Light chains | sFLC ratio ≥100 |
| **M** - MRI | >1 focal lesion ≥5mm |

### 3. Risk Stratification (R-ISS)

| Stage | Criteria | Median OS |
|-------|----------|-----------|
| **I** | ISS I + standard-risk cytogenetics + LDH normal | ~80 months |
| **II** | Not I or III | ~50 months |
| **III** | ISS III + high-risk cytogenetics or LDH elevated | ~25 months |

### 4. AI-Powered Analysis

| Application | Method | Purpose |
|-------------|--------|---------|
| **PET/CT analysis** | Deep learning | Bone lesion detection |
| **Survival prediction** | ML models | Personalized prognosis |
| **MRD assessment** | Flow cytometry AI | Response depth |
| **Treatment selection** | Decision support | Precision therapy |

---

## Usage

### Example Prompt

```text
Evaluate this newly diagnosed multiple myeloma patient and provide treatment recommendations.

Patient: 65-year-old male

Findings:
- M-protein: IgG kappa 4.2 g/dL
- Serum FLC: Kappa 856 mg/L, Lambda 12 mg/L (ratio 71)
- Bone marrow: 45% clonal plasma cells
- Cytogenetics: t(4;14) detected by FISH
- Hemoglobin: 9.2 g/dL
- Creatinine: 1.8 mg/dL
- Calcium: 10.8 mg/dL
- LDH: Normal
- Beta-2 microglobulin: 5.8 mg/L
- Albumin: 3.2 g/dL
- PET/CT: Multiple lytic lesions

Tasks:
1. Confirm diagnosis and stage
2. Risk stratification
3. Treatment recommendation
4. MRD monitoring strategy
5. Prognosis assessment
```

### Expected Output

```
## Multiple Myeloma Evaluation Report

### Patient Information
- **Age:** 65 years
- **Sex:** Male
- **Performance status:** [Not provided]

### Diagnostic Evaluation

#### Monoclonal Protein Assessment

| Parameter | Value | Interpretation |
|-----------|-------|----------------|
| M-protein type | IgG kappa | Heavy chain type |
| M-protein | 4.2 g/dL | Elevated |
| Serum kappa FLC | 856 mg/L | Markedly elevated |
| Serum lambda FLC | 12 mg/L | Normal |
| FLC ratio | 71.3 | Abnormal (>100 = SLiM) |

#### Bone Marrow Findings

| Finding | Result |
|---------|--------|
| Plasma cells | 45% |
| Clonality | Kappa-restricted |
| Cytogenetics | t(4;14) by FISH |

#### End-Organ Damage (CRAB)

| Criterion | Finding | Status |
|-----------|---------|--------|
| **C**alcium | 10.8 mg/dL | Not elevated |
| **R**enal | Cr 1.8 mg/dL | Borderline |
| **A**nemia | Hgb 9.2 g/dL | Present (<10) |
| **B**one | Multiple lytic lesions | Present |

**CRAB criteria met:** Anemia + Bone lesions

#### SLiM Criteria

| Criterion | Finding | Status |
|-----------|---------|--------|
| **S** ≥60% PCs | 45% | Not met |
| **Li** FLC ratio ≥100 | 71 | Not met |
| **M** MRI lesions | PET positive | Met |

### Diagnosis

```
┌─────────────────────────────────────────────────────────────┐
│                                                             │
│   DIAGNOSIS: MULTIPLE MYELOMA                               │
│                                                             │
│   IgG Kappa Subtype                                         │
│                                                             │
│   Diagnostic criteria met:                                  │
│   ✓ Clonal plasma cells ≥10% (45%)                         │
│   ✓ CRAB criteria (Anemia, Bone lesions)                   │
│                                                             │
│   Cytogenetic abnormality: t(4;14) = HIGH RISK             │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

### Risk Stratification

#### ISS Stage

| Parameter | Value | Points |
|-----------|-------|--------|
| Beta-2 microglobulin | 5.8 mg/L | ≥5.5 = Stage III |
| Albumin | 3.2 g/dL | <3.5 |

**ISS Stage: III** (B2M ≥5.5)

#### Cytogenetic Risk

| Abnormality | Status | Risk |
|-------------|--------|------|
| **t(4;14)** | Present | **High risk** |
| t(14;16) | Not detected | - |
| t(14;20) | Not detected | - |
| del(17p) | Not detected | - |
| 1q gain | [Recommend testing] | - |

**Cytogenetic risk category:** HIGH RISK (t(4;14))

#### R-ISS Stage

| Criterion | Status |
|-----------|--------|
| ISS Stage | III |
| High-risk cytogenetics | Yes (t(4;14)) |
| LDH | Normal |

```
┌─────────────────────────────────────────────────────────────┐
│                                                             │
│   R-ISS STAGE: III                                          │
│                                                             │
│   ISS III + High-risk cytogenetics [t(4;14)]               │
│                                                             │
│   Median OS: ~25-30 months (historical)                     │
│   Note: Newer therapies significantly improve outcomes      │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

#### Additional Risk Factors

| Factor | Status | Impact |
|--------|--------|--------|
| Age >75 | No | Favorable |
| ISS III | Yes | Adverse |
| High-risk cytogenetics | Yes | Adverse |
| Extramedullary disease | No | Favorable |
| High LDH | No | Favorable |
| del(17p) | No | Favorable (absent) |

### Treatment Recommendations

#### Transplant Eligibility Assessment

| Factor | Status | Transplant eligible |
|--------|--------|---------------------|
| Age | 65 years | Yes |
| Performance status | [Assess] | - |
| Comorbidities | [Assess] | - |
| Renal function | Borderline | Yes (if stable) |

**Assessment:** LIKELY TRANSPLANT ELIGIBLE

#### Recommended Induction Regimen

**For High-Risk Transplant-Eligible MM:**

```
VRd (Bortezomib, Lenalidomide, Dexamethasone)
+ Daratumumab (Dara-VRd)

Regimen: Dara-VRd
- Daratumumab 16 mg/kg weekly x 8, then q2wks, then monthly
- Bortezomib 1.3 mg/m² SC days 1, 4, 8, 11
- Lenalidomide 25 mg days 1-14
- Dexamethasone 40 mg weekly

Cycles: 4-6 induction cycles

Rationale:
- GRIFFIN trial: Dara-VRd superior to VRd
- High-risk cytogenetics benefit from quadruplet
- t(4;14) may benefit from bortezomib
```

#### Treatment Sequence

```
1. INDUCTION (4-6 cycles)
   Dara-VRd

2. STEM CELL COLLECTION
   - Collect after 4+ cycles
   - Plead + G-CSF mobilization
   - Target: ≥4 x 10^6 CD34+/kg (2 transplants)

3. AUTOLOGOUS SCT
   - High-dose melphalan 200 mg/m²
   - Consider tandem SCT for high-risk

4. CONSOLIDATION (2-4 cycles post-SCT)
   - Continue Dara-VRd or VRd

5. MAINTENANCE
   - Lenalidomide + bortezomib (high-risk)
   - Continue until progression
```

### MRD Monitoring Strategy

#### MRD Assessment Methods

| Method | Sensitivity | Timing |
|--------|-------------|--------|
| **Next-gen flow** | 10^-5 | Post-induction, post-SCT |
| **NGS (clonoSEQ)** | 10^-6 | Post-SCT, during maintenance |
| **PET/CT** | Anatomic | Post-induction, post-SCT |

#### Monitoring Schedule

| Timepoint | Assessment |
|-----------|------------|
| Post-induction | MRD (flow/NGS), M-protein, PET/CT |
| Day +100 post-SCT | MRD, M-protein |
| 6 months post-SCT | MRD, M-protein |
| Annually on maintenance | MRD if in CR |

**MRD Goal:** MRD negativity (10^-5 or better)
- MRD negative high-risk: Better outcomes
- MRD positive: Consider intensification

### Prognosis

#### Expected Outcomes (with modern therapy)

| Outcome | Standard Risk | High Risk (t(4;14)) |
|---------|---------------|---------------------|
| ORR | >95% | 90% |
| CR/sCR rate | 60-70% | 50-60% |
| MRD negativity | 50-60% | 30-40% |
| 5-year PFS | 60-70% | 40-50% |
| 5-year OS | 70-80% | 55-65% |

**Note:** Outcomes improving with quadruplet induction and extended maintenance.

### Additional Recommendations

| Recommendation | Rationale |
|----------------|-----------|
| Bone health | Zoledronic acid monthly x 2 years |
| VTE prophylaxis | Aspirin or anticoagulation |
| Infection prophylaxis | Acyclovir, consider IVIG if recurrent |
| Vaccinations | Pneumococcal, influenza, COVID-19 |

### Summary

```
DIAGNOSIS: Multiple Myeloma, IgG Kappa
STAGE: R-ISS III (High-risk)
CYTOGENETICS: t(4;14) - High risk

TREATMENT PLAN:
1. Dara-VRd induction x 4-6 cycles
2. Autologous SCT (consider tandem)
3. Dara-VRd or VRd consolidation
4. Lenalidomide + bortezomib maintenance

MRD MONITORING:
- Target MRD negativity
- NGS-based monitoring post-SCT

PROGNOSIS:
- 5-year PFS: ~40-50% (high-risk with modern therapy)
- 5-year OS: ~55-65%
```
```

### LLM Agent Integration

```python
@tool
def diagnose_plasma_cell_neoplasm(
    m_protein: dict,
    bone_marrow: dict,
    imaging: dict,
    labs: dict
) -> str:
    """
    Diagnoses plasma cell neoplasm (MGUS/SMM/MM).

    Args:
        m_protein: M-protein and FLC values
        bone_marrow: Plasma cell percentage, clonality
        imaging: PET/CT, MRI findings
        labs: CRAB criteria labs

    Returns:
        Diagnosis with staging
    """
    pass


@tool
def stratify_myeloma_risk(
    iss_stage: int,
    cytogenetics: list[str],
    ldh: str,
    other_factors: dict = None
) -> str:
    """
    Calculates R-ISS and other risk scores.

    Args:
        iss_stage: ISS stage (1, 2, 3)
        cytogenetics: Cytogenetic abnormalities
        ldh: LDH status (normal/elevated)
        other_factors: Additional risk factors

    Returns:
        R-ISS stage with prognosis
    """
    pass


@tool
def recommend_myeloma_treatment(
    risk_category: str,
    transplant_eligible: bool,
    patient_age: int,
    comorbidities: list[str] = None
) -> str:
    """
    Recommends myeloma treatment regimen.

    Args:
        risk_category: Standard vs high risk
        transplant_eligible: Transplant eligibility
        patient_age: Patient age
        comorbidities: Relevant comorbidities

    Returns:
        Treatment recommendation with sequence
    """
    pass


@tool
def assess_myeloma_response(
    baseline_values: dict,
    current_values: dict,
    mrd_status: str = None
) -> str:
    """
    Assesses myeloma treatment response (IMWG criteria).

    Args:
        baseline_values: Baseline M-protein, FLC, PCs
        current_values: Current values
        mrd_status: MRD result if available

    Returns:
        Response category (sCR, CR, VGPR, PR, etc.)
    """
    pass
```

---

## Prerequisites

### Required Data

| Data Type | Purpose |
|-----------|---------|
| M-protein quantification | Diagnosis, monitoring |
| Serum free light chains | Diagnosis, prognosis |
| Bone marrow biopsy | Plasma cell %, cytogenetics |
| FISH panel | Risk stratification |
| Imaging (PET/CT, MRI) | Bone disease |

### Dependencies

```
pandas>=2.0
numpy>=1.24
```

---

## Clinical Applications

### Diagnosis
- MGUS vs SMM vs MM differentiation
- CRAB/SLiM criteria assessment
- Subtype classification

### Risk Stratification
- R-ISS staging
- Cytogenetic risk
- Survival prediction

### Treatment
- Transplant eligibility
- Regimen selection
- Maintenance planning

### Monitoring
- Response assessment (IMWG)
- MRD monitoring
- Relapse detection

---

## Related Skills

- **Flow Cytometry AI Agent:** Plasma cell detection
- **MRD Detection Agent:** MRD assessment
- **Bone Imaging AI Agent:** Lytic lesion detection
- **Cytogenetics Agent:** FISH interpretation

---

## References

- **Rajkumar et al. (2014):** "International Myeloma Working Group updated criteria for the diagnosis of multiple myeloma." *Lancet Oncology*
- **Palumbo et al. (2015):** "Revised International Staging System for Multiple Myeloma." *JCO*
- **Voorhees et al. (2020):** "Daratumumab, lenalidomide, bortezomib, and dexamethasone for transplant-eligible newly diagnosed multiple myeloma: GRIFFIN." *Blood*

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu
