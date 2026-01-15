# Clonal Hematopoiesis Analysis Agent

**ID:** `biomedical.hematology.clonal_hematopoiesis`
**Version:** 1.0.0
**Status:** Production
**Category:** Hematology / Clonal Hematopoiesis / Risk Assessment

---

## Overview

The **Clonal Hematopoiesis Analysis Agent** identifies and characterizes clonal hematopoiesis of indeterminate potential (CHIP) and clonal cytopenia of undetermined significance (CCUS). Clonal hematopoiesis refers to the age-related acquisition of somatic mutations in hematopoietic stem cells, which is associated with increased risk of hematologic malignancy, cardiovascular disease, and overall mortality.

This agent integrates genomic analysis, clinical risk assessment, and AI-powered progression modeling to provide actionable insights for clinical management of patients with clonal hematopoiesis.

---

## Key Capabilities

### 1. CHIP/CCUS Classification

| Entity | VAF Threshold | Cytopenias | Morphologic Dysplasia |
|--------|---------------|------------|----------------------|
| **CHIP** | ≥2% | None | None |
| **CCUS** | ≥2% | Present | None/minimal |
| **Low-risk MDS** | ≥2% | Present | Present |
| **High-risk MDS** | ≥2% | Present | Present + excess blasts |

### 2. Common CHIP Mutations

| Gene | Function | Frequency | Clinical Significance |
|------|----------|-----------|----------------------|
| **DNMT3A** | DNA methylation | 40-60% | Lowest malignancy risk |
| **TET2** | DNA demethylation | 20-30% | Cardiovascular risk |
| **ASXL1** | Chromatin modification | 10-15% | Higher malignancy risk |
| **JAK2** | Signaling (V617F) | 5-10% | MPN association |
| **TP53** | Tumor suppressor | 5-8% | Highest malignancy risk |
| **SF3B1** | Splicing | 3-5% | MDS-RS association |
| **SRSF2** | Splicing | 3-5% | MDS/CMML association |
| **PPM1D** | DNA damage response | 2-4% | Therapy-related |

### 3. Risk Stratification Models

| Model | Variables | Outcome |
|-------|-----------|---------|
| **CHIP-score** | Gene, VAF, age | 10-year AML risk |
| **CCUS progression** | Mutations, cytopenias | 5-year MDS risk |
| **CHRS** | Cardiovascular factors | CV event risk |
| **MDS-flow AI** | Flow + genetics | CHIP → MDS progression |

### 4. Clinical Associations

- **Hematologic malignancy:** AML, MDS, MPN
- **Cardiovascular disease:** Atherosclerosis, heart failure
- **Inflammation:** Cytokine activation, inflammasome
- **Therapy-related:** Post-chemo/radiation clones
- **All-cause mortality:** Independent risk factor

---

## Usage

### Example Prompt

```text
Analyze this NGS panel result for a 68-year-old patient with unexplained mild anemia.

Patient information:
- Age: 68, Female
- Hemoglobin: 10.8 g/dL (mild anemia)
- MCV: 102 fL (macrocytic)
- Platelets: Normal
- WBC: Normal
- No prior chemotherapy
- History of CAD, statin therapy

NGS findings:
- DNMT3A R882H VAF 18%
- TET2 Q1034* VAF 12%

Tasks:
1. Classify as CHIP vs CCUS
2. Assess progression risk to MDS/AML
3. Evaluate cardiovascular risk implications
4. Recommend surveillance strategy
5. Provide clinical management guidance
```

### Expected Output

```
## Clonal Hematopoiesis Analysis Report

### Patient Profile
- **Age:** 68 years
- **Sex:** Female
- **Chief complaint:** Unexplained mild anemia
- **Comorbidities:** Coronary artery disease

### Laboratory Findings

| Parameter | Value | Reference | Interpretation |
|-----------|-------|-----------|----------------|
| Hemoglobin | 10.8 g/dL | 12-16 | Mild anemia |
| MCV | 102 fL | 80-100 | Macrocytic |
| Platelets | Normal | 150-400 | Normal |
| WBC | Normal | 4-11 | Normal |
| Reticulocytes | [Recommend] | - | Pending |
| B12/Folate | [Recommend] | - | Pending |

### Genomic Findings

#### Detected Mutations

| Gene | Mutation | VAF | Pathogenicity | CHIP Gene |
|------|----------|-----|---------------|-----------|
| **DNMT3A** | R882H | 18% | Pathogenic | Yes |
| **TET2** | Q1034* | 12% | Pathogenic | Yes |

#### Mutation Analysis

**DNMT3A R882H (VAF 18%):**
- Most common DNMT3A hotspot
- Loss of DNA methyltransferase activity
- Associated with clonal expansion
- Lower malignancy risk among CHIP genes

**TET2 Q1034* (VAF 12%):**
- Truncating mutation (stop gain)
- Loss of TET2 function
- Associated with myeloid bias
- Cardiovascular risk implications

### Classification

```
┌─────────────────────────────────────────────────────────┐
│                                                         │
│   CLASSIFICATION: CCUS                                  │
│   (Clonal Cytopenia of Undetermined Significance)       │
│                                                         │
│   Criteria met:                                         │
│   ✓ Somatic mutation VAF ≥2%                           │
│   ✓ Cytopenia present (Hgb 10.8)                       │
│   ✗ Morphologic dysplasia (bone marrow not done)       │
│                                                         │
│   Note: Bone marrow recommended to exclude MDS          │
│                                                         │
└─────────────────────────────────────────────────────────┘
```

### Progression Risk Assessment

#### CHIP to MDS/AML Risk

**Risk factors present:**
| Factor | Status | Impact |
|--------|--------|--------|
| Age >65 | Yes | Increased risk |
| VAF >10% | Yes (18% DNMT3A) | Increased risk |
| Multiple mutations | Yes (2) | Increased risk |
| High-risk genes | No (DNMT3A/TET2 = moderate) | Moderate |
| TP53/SRSF2/ASXL1 | No | Favorable |
| Cytopenia | Yes | Increased risk (CCUS) |

**Estimated 10-year progression risk:**

| Outcome | Risk | Confidence |
|---------|------|------------|
| MDS | 25-35% | Moderate |
| AML | 5-10% | Moderate |
| Stable CCUS | 55-70% | Moderate |

**Risk category:** INTERMEDIATE

#### Comparison to Population Risks

| Population | 10-year MDS/AML Risk |
|------------|---------------------|
| No mutations (age 68) | 0.5% |
| Single DNMT3A VAF <10% | 2-3% |
| **This patient (CCUS, 2 mutations)** | **25-35%** |
| CCUS with ASXL1/SRSF2 | 40-50% |
| CCUS with TP53 | 60-70% |

### Cardiovascular Risk Assessment

#### CHIP and Cardiovascular Disease

| Finding | Clinical Implication |
|---------|---------------------|
| **TET2 mutation** | 1.9x increased CVD risk |
| **Pre-existing CAD** | Compounding risk factor |
| **Inflammatory pathway** | TET2 loss → IL-1β/IL-6 activation |

**CHRS (Clonal Hematopoiesis Risk Score) Assessment:**

| Factor | Points |
|--------|--------|
| Age >60 | +1 |
| TET2 mutation | +2 |
| VAF >10% | +1 |
| Existing CAD | +2 |
| **Total** | **6 points** |

**Interpretation:** High cardiovascular risk category
- 2.4x increased risk of MACE vs no CHIP
- Consider aggressive cardiovascular risk modification

### Recommended Workup

#### Immediate (within 2 weeks)

| Test | Rationale | Priority |
|------|-----------|----------|
| **Bone marrow biopsy** | Exclude MDS | High |
| **Reticulocyte count** | Anemia workup | High |
| **B12, folate, iron studies** | Exclude other causes | High |
| **LDH, haptoglobin** | Hemolysis screen | Moderate |

#### If Bone Marrow Normal (confirms CCUS)

| Test | Frequency | Purpose |
|------|-----------|---------|
| CBC with differential | Q3 months x1 year | Monitor cytopenias |
| Peripheral smear | Q6 months | Dysplasia screening |
| Repeat NGS panel | Q12 months | Clone evolution |

### Management Recommendations

#### Hematologic Management

```
SURVEILLANCE PROTOCOL (CCUS)

Year 1:
- CBC every 3 months
- Clinical assessment for symptoms
- Peripheral smear review if cytopenias worsen

Year 2+:
- CBC every 4-6 months
- Annual NGS panel (track clone evolution)
- Annual consideration for bone marrow

INTERVENTION TRIGGERS:
- Hgb <10 g/dL (progressive anemia)
- Platelet <100 K
- ANC <1.5 K
- New mutations (especially TP53, SF3B1, SRSF2)
- VAF increase >50%
- Morphologic dysplasia
```

#### Cardiovascular Management

| Recommendation | Rationale | Evidence |
|----------------|-----------|----------|
| **Statin optimization** | Continue/intensify | CHIP + CVD risk |
| **Aspirin** | Continue | Secondary prevention |
| **BP control** | Target <130/80 | Risk modification |
| **Inflammatory markers** | Check hsCRP | CHIP inflammation |
| **Consider SGLT2i/GLP-1** | CV benefit if diabetic | Emerging data |

**Emerging therapies (research context):**
- IL-1β inhibition (canakinumab): Reduced CV events in TET2 CHIP
- IL-6 inhibition: Under investigation

### Clinical Summary

```
FINAL ASSESSMENT:

DIAGNOSIS: CCUS (Clonal Cytopenia of Undetermined Significance)
- Two CHIP mutations: DNMT3A R882H (18%), TET2 Q1034* (12%)
- Mild macrocytic anemia
- Bone marrow pending to exclude MDS

RISK STRATIFICATION:
- MDS/AML progression: Intermediate (25-35% 10-year)
- Cardiovascular: High (TET2 + pre-existing CAD)

KEY RECOMMENDATIONS:
1. Bone marrow biopsy to exclude MDS
2. Q3-month CBC surveillance (year 1)
3. Annual NGS to monitor clone evolution
4. Aggressive cardiovascular risk modification
5. Consider hematology-cardiology collaboration
6. Patient education regarding signs of progression
```
```

### LLM Agent Integration

```python
@tool
def classify_clonal_hematopoiesis(
    mutations: list[dict],
    cbc_values: dict,
    bone_marrow_done: bool = False
) -> str:
    """
    Classifies as CHIP, CCUS, or MDS.

    Args:
        mutations: List of detected mutations with VAF
        cbc_values: Complete blood count values
        bone_marrow_done: Whether bone marrow was evaluated

    Returns:
        Classification with criteria evaluation
    """
    pass


@tool
def assess_chip_progression_risk(
    mutations: list[dict],
    patient_age: int,
    cytopenias: bool = False
) -> str:
    """
    Assesses risk of progression to MDS/AML.

    Args:
        mutations: Detected CHIP mutations
        patient_age: Patient age in years
        cytopenias: Presence of cytopenias

    Returns:
        Progression risk with confidence intervals
    """
    pass


@tool
def assess_chip_cardiovascular_risk(
    mutations: list[dict],
    cv_history: dict
) -> str:
    """
    Assesses cardiovascular risk from CHIP.

    Args:
        mutations: Detected CHIP mutations
        cv_history: Cardiovascular history and risk factors

    Returns:
        Cardiovascular risk assessment with recommendations
    """
    pass


@tool
def generate_chip_surveillance_plan(
    classification: str,
    risk_level: str,
    mutations: list[dict]
) -> str:
    """
    Generates surveillance protocol for CHIP/CCUS.

    Args:
        classification: CHIP, CCUS, or low-risk MDS
        risk_level: Low, intermediate, or high
        mutations: Specific mutations detected

    Returns:
        Surveillance protocol with monitoring schedule
    """
    pass
```

---

## Prerequisites

### Required Data

| Data Type | Purpose |
|-----------|---------|
| NGS panel results | Mutation detection |
| Complete blood count | Cytopenia assessment |
| Clinical history | Risk stratification |
| Bone marrow (if available) | Dysplasia evaluation |

### Dependencies

```
pandas>=2.0
numpy>=1.24
scipy>=1.11
```

---

## Methodology

### CHIP Analysis Pipeline

```
NGS Panel Results
    ↓
Mutation Filtering
├── VAF ≥2%
├── CHIP-associated genes
└── Pathogenicity assessment
    ↓
Clinical Correlation
├── CBC analysis
├── Cytopenia identification
└── Bone marrow status
    ↓
Classification
├── CHIP (no cytopenias)
├── CCUS (cytopenias, no dysplasia)
└── MDS (dysplasia present)
    ↓
Risk Stratification
├── Hematologic progression
├── Cardiovascular risk
└── Overall prognosis
    ↓
Management Plan
├── Surveillance protocol
├── Intervention triggers
└── Specialist referrals
```

---

## Clinical Applications

### Hematology
- MDS screening in cytopenic patients
- Post-chemotherapy surveillance
- Bone marrow donor screening
- Unexplained cytopenias workup

### Cardiology
- Novel CV risk factor assessment
- Inflammatory pathway targeting
- Risk-adapted management

### Oncology
- Therapy-related clone monitoring
- Cancer survivorship care
- Pre-stem cell transplant evaluation

---

## Related Skills

- **MDS Diagnosis Agent:** MDS classification
- **Flow Cytometry AI Agent:** MDS flow scoring
- **AML Classification Agent:** CHIP → AML progression
- **Cardiovascular Risk Agent:** CHIP-CV integration

---

## References

- **Jaiswal et al. (2014):** "Age-related clonal hematopoiesis associated with adverse outcomes." *NEJM*
- **Steensma et al. (2015):** "Clonal hematopoiesis of indeterminate potential and its distinction from myelodysplastic syndromes." *Blood*
- **Fuster et al. (2017):** "Clonal hematopoiesis associated with TET2 deficiency accelerates atherosclerosis development in mice." *Science*
- [CHIP Clinical Guidelines](https://www.hematology.org/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu
