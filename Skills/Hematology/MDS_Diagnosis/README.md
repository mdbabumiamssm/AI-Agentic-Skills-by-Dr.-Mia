# MDS Diagnosis Agent

**ID:** `biomedical.hematology.mds_diagnosis`
**Version:** 1.0.0
**Status:** Production
**Category:** Hematology / Myelodysplastic Syndromes / Diagnosis

---

## Overview

The **MDS Diagnosis Agent** provides comprehensive diagnostic classification and risk stratification for myelodysplastic syndromes (MDS) according to WHO 2022 and ICC 2022 criteria. MDS encompasses a heterogeneous group of clonal hematopoietic stem cell disorders characterized by ineffective hematopoiesis, cytopenias, and increased risk of transformation to acute myeloid leukemia.

This agent integrates morphologic evaluation, flow cytometry AI scoring, cytogenetics, and molecular analysis to deliver accurate MDS diagnosis with IPSS-M risk stratification for treatment planning.

---

## Key Capabilities

### 1. WHO 2022 MDS Classification

| Category | Defining Features | Blast Threshold |
|----------|-------------------|-----------------|
| **MDS with low blasts (MDS-LB)** | Cytopenias, dysplasia, <5% blasts | <5% BM, <2% PB |
| **MDS with increased blasts (MDS-IB)** | 5-19% blasts or Auer rods | 5-19% BM |
| **MDS with defining genetic abnormality** | SF3B1, del(5q), TP53 | Variable |
| **MDS, hypoplastic** | <25% cellularity | <5% BM |
| **MDS with fibrosis** | Grade 2-3 fibrosis | Variable |

### 2. Key Genetic Abnormalities

| Abnormality | Frequency | Prognostic Impact | WHO 2022 Defining |
|-------------|-----------|-------------------|-------------------|
| **SF3B1** | 20-30% | Favorable | Yes (MDS-SF3B1) |
| **del(5q)** | 10-15% | Favorable (isolated) | Yes (MDS with del(5q)) |
| **TP53** | 5-10% | Very adverse | Yes (MDS with biallelic TP53) |
| **Complex karyotype** | 10-15% | Adverse | Risk modifier |
| **ASXL1** | 15-20% | Adverse | Risk modifier |
| **RUNX1** | 10% | Adverse | Risk modifier |
| **EZH2** | 5% | Adverse | Risk modifier |

### 3. Risk Stratification Models

| Model | Variables | Categories | Use |
|-------|-----------|------------|-----|
| **IPSS** | Blasts, karyotype, cytopenias | 4 | Historical |
| **IPSS-R** | Blasts, karyotype, cytopenias (refined) | 5 | Standard |
| **IPSS-M** | IPSS-R + molecular | 6 | Current standard |
| **MDS-flow AI** | Flow cytometry features | Continuous | Adjunct |

### 4. AI-Powered Diagnostics

| Tool | Application | Performance |
|------|-------------|-------------|
| **MDS-flow AI** | Flow cytometry scoring | AUC 0.935 |
| **Ogata score** | Immunophenotype dysplasia | Sensitivity 91.8% |
| **Morphology AI** | Dysplasia quantification | Under development |

---

## Usage

### Example Prompt

```text
Evaluate this patient for MDS and provide risk stratification.

Patient: 72-year-old female with refractory anemia

Findings:
- CBC: Hgb 8.2 g/dL, MCV 108 fL, Plt 98K, WBC 3.2K
- Bone marrow: 45% cellularity, 3% blasts
- Morphology: Dyserythropoiesis (megaloblastoid changes, nuclear budding), dysgranulopoiesis (hypogranular neutrophils)
- Cytogenetics: 46,XX,del(5)(q13q33)[18]/46,XX[2]
- Molecular: SF3B1 K700E (VAF 38%), TET2 Q1034* (VAF 22%)

Tasks:
1. Confirm MDS diagnosis
2. WHO 2022 classification
3. IPSS-M risk stratification
4. Treatment recommendations
5. Monitoring plan
```

### Expected Output

```
## MDS Diagnosis and Risk Stratification Report

### Patient Information
- **Age:** 72 years
- **Sex:** Female
- **Presentation:** Refractory anemia

### Diagnostic Findings

#### Peripheral Blood

| Parameter | Value | Reference | Interpretation |
|-----------|-------|-----------|----------------|
| Hemoglobin | 8.2 g/dL | 12-16 | Moderate anemia |
| MCV | 108 fL | 80-100 | Macrocytic |
| Platelets | 98 K/μL | 150-400 | Mild thrombocytopenia |
| WBC | 3.2 K/μL | 4-11 | Leukopenia |
| ANC | 1.4 K/μL | 1.5-8 | Neutropenia |

**Cytopenias present:** 3/3 (anemia, thrombocytopenia, neutropenia)

#### Bone Marrow

| Finding | Result | Significance |
|---------|--------|--------------|
| Cellularity | 45% | Normal for age (expected 30-40%) |
| Blasts (BM) | 3% | <5% = low blasts |
| Blasts (PB) | 1% | <2% |
| Ring sideroblasts | 32% | ≥15% = significant |

#### Dysplasia Assessment

| Lineage | Features | % Dysplastic | Threshold |
|---------|----------|--------------|-----------|
| **Erythroid** | Megaloblastoid, nuclear budding | 45% | ≥10% |
| **Granulocytic** | Hypogranular neutrophils | 22% | ≥10% |
| **Megakaryocytic** | Not specified | - | ≥10% |

**Multilineage dysplasia:** 2 lineages confirmed

#### Cytogenetics

```
46,XX,del(5)(q13q33)[18]/46,XX[2]

Interpretation:
- del(5q) isolated abnormality (90% of metaphases)
- MDS-defining cytogenetic abnormality
- WHO 2022: MDS with isolated del(5q)
- IPSS-R cytogenetic risk: GOOD
```

#### Molecular Profile

| Gene | Mutation | VAF | MDS-Associated | Impact |
|------|----------|-----|----------------|--------|
| **SF3B1** | K700E | 38% | Yes | Favorable |
| **TET2** | Q1034* | 22% | Yes | Neutral |

### WHO 2022 Classification

```
┌─────────────────────────────────────────────────────────────┐
│                                                             │
│   WHO 2022 DIAGNOSIS:                                       │
│                                                             │
│   MDS WITH LOW BLASTS AND ISOLATED del(5q)                 │
│   (MDS-5q)                                                  │
│                                                             │
│   Additional finding: SF3B1 mutation                        │
│   - Compatible with MDS-5q                                  │
│   - Associated with ring sideroblasts (32%)                │
│                                                             │
│   Note: del(5q) is the primary defining feature            │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

#### Classification Criteria Met

| Criterion | Finding | Met |
|-----------|---------|-----|
| Cytopenia(s) | Pancytopenia | ✓ |
| Dysplasia ≥10% in ≥1 lineage | Erythroid 45%, Myeloid 22% | ✓ |
| Blasts <5% BM | 3% | ✓ |
| Isolated del(5q) | Present (90%) | ✓ |
| No AML-defining mutations | Not detected | ✓ |

### Risk Stratification

#### IPSS-R Score

| Variable | Finding | Points |
|----------|---------|--------|
| Cytogenetic risk | Good (del(5q) isolated) | 0 |
| BM blast % | 0-2% (use 3%) | 0.5 |
| Hemoglobin | <8 g/dL | 1.5 |
| Platelets | 50-100 K | 0.5 |
| ANC | <0.8 K | 0.5 |

**IPSS-R Score: 3.0**
**IPSS-R Category: INTERMEDIATE**

#### IPSS-M Score (Molecular-Adjusted)

| Factor | Value | Adjustment |
|--------|-------|------------|
| IPSS-R baseline | 3.0 | - |
| SF3B1 mutation | Present | -0.5 (favorable) |
| TET2 mutation | Present | 0 (neutral) |
| TP53 | Absent | 0 |
| ASXL1 | Absent | 0 |

**IPSS-M Score: 2.5**
**IPSS-M Category: LOW**

#### Risk Summary

| Model | Category | Median Survival | AML Risk (5y) |
|-------|----------|-----------------|---------------|
| IPSS-R | Intermediate | 3.0 years | 25% |
| **IPSS-M** | **Low** | **5.5 years** | **12%** |

**Note:** IPSS-M is more accurate with molecular data available.

### Treatment Recommendations

#### First-Line Therapy

| Option | Indication | Evidence | Recommendation |
|--------|------------|----------|----------------|
| **Lenalidomide** | del(5q) MDS | MDS-003, MDS-004 | **First choice** |
| ESA | Anemia, low EPO | ESAM | If EPO <500 |
| Luspatercept | Ring sideroblasts | MEDALIST | Second-line |

**Recommended treatment:**

```
LENALIDOMIDE 10 mg daily x 21 days, 28-day cycle

Rationale:
- del(5q) isolated = excellent lenalidomide response
- Expected transfusion independence: 67%
- Expected cytogenetic response: 45%

Monitoring:
- CBC weekly x 8 weeks, then q2-4 weeks
- Watch for: Neutropenia, thrombocytopenia
- DVT prophylaxis consideration
```

#### Alternative/Subsequent Therapy

| Scenario | Option |
|----------|--------|
| Lenalidomide failure | Luspatercept (SF3B1+, RS+) |
| Progression to higher-risk | Azacitidine |
| Young, fit, high-risk | Allogeneic SCT |

### Monitoring Plan

#### Response Assessment

| Timepoint | Assessment | Criteria |
|-----------|------------|----------|
| 8 weeks | CBC, transfusion requirement | Hematologic response |
| 16 weeks | BM biopsy | Cytogenetic response |
| 6 months | Full restaging | Complete response |
| Ongoing | CBC q4-8 weeks | Durability |

#### Surveillance for Progression

| Sign | Action |
|------|--------|
| Rising blasts (>5%) | Restage, reconsider treatment |
| New cytogenetic abnormalities | Higher-risk therapy |
| TP53 acquisition | Poor prognosis, SCT if eligible |
| Transfusion dependence | Consider luspatercept |

### Prognosis

#### MDS with del(5q) + SF3B1

| Outcome | Expectation |
|---------|-------------|
| Median survival | 5-7 years |
| Lenalidomide response | 67% transfusion independence |
| Cytogenetic remission | 45% |
| AML transformation | 15% at 5 years |
| Favorable modifiers | SF3B1 mutation, isolated del(5q) |

### Summary

```
FINAL DIAGNOSIS:
MDS with low blasts and isolated del(5q) (MDS-5q)
- SF3B1 mutated
- IPSS-M: Low risk

TREATMENT:
Lenalidomide 10 mg daily (21/28 cycle)

PROGNOSIS:
- Favorable for del(5q) MDS
- Expected response to lenalidomide
- Median survival ~5-7 years

MONITORING:
- CBC weekly x8 weeks, then q2-4 weeks
- BM restaging at 16 weeks
- Watch for lenalidomide-related cytopenias
```
```

### LLM Agent Integration

```python
@tool
def diagnose_mds(
    morphology: dict,
    flow_cytometry: dict,
    cytogenetics: str,
    molecular: dict,
    blast_percentage: float
) -> str:
    """
    Diagnoses MDS according to WHO 2022.

    Args:
        morphology: Dysplasia findings
        flow_cytometry: Immunophenotype (for Ogata score)
        cytogenetics: Karyotype result
        molecular: Molecular panel results
        blast_percentage: BM blast percentage

    Returns:
        WHO 2022 MDS classification
    """
    pass


@tool
def calculate_ipss_m(
    cytopenias: dict,
    cytogenetics: str,
    blasts: float,
    molecular: dict
) -> str:
    """
    Calculates IPSS-M risk score.

    Args:
        cytopenias: CBC values
        cytogenetics: Cytogenetic risk group
        blasts: Bone marrow blast percentage
        molecular: Molecular findings

    Returns:
        IPSS-M score and risk category
    """
    pass


@tool
def calculate_mds_flow_score(
    flow_data: str
) -> str:
    """
    Calculates AI-assisted MDS flow score.

    Args:
        flow_data: Path to flow cytometry data

    Returns:
        MDS flow score with interpretation
    """
    pass


@tool
def recommend_mds_treatment(
    who_classification: str,
    ipss_m_risk: str,
    patient_age: int,
    molecular_features: dict
) -> str:
    """
    Recommends MDS treatment.

    Args:
        who_classification: WHO 2022 diagnosis
        ipss_m_risk: IPSS-M risk category
        patient_age: Patient age
        molecular_features: Actionable mutations

    Returns:
        Treatment recommendations
    """
    pass
```

---

## Prerequisites

### Required Data

| Data Type | Purpose |
|-----------|---------|
| CBC | Cytopenia assessment |
| Bone marrow | Blasts, dysplasia |
| Cytogenetics | Risk stratification |
| Molecular panel | IPSS-M, classification |

### Dependencies

```
pandas>=2.0
numpy>=1.24
```

---

## Clinical Applications

### Diagnosis
- WHO 2022 classification
- Differential from CHIP/CCUS
- AML exclusion

### Prognosis
- IPSS-M risk stratification
- Survival estimation
- AML transformation risk

### Treatment
- Risk-adapted therapy
- Targeted therapy selection
- SCT decision support

---

## Related Skills

- **CHIP Analysis Agent:** CHIP vs MDS
- **AML Classification Agent:** MDS vs AML
- **Flow Cytometry AI Agent:** MDS-flow scoring
- **Cytogenetics Agent:** Karyotype interpretation

---

## References

- **Khoury et al. (2022):** "WHO Classification of Tumours of Haematopoietic and Lymphoid Tissues." 5th edition
- **Bernard et al. (2022):** "Molecular International Prognostic Scoring System for Myelodysplastic Syndromes." *NEJM Evidence*
- [IPSS-M Calculator](https://www.mds-risk-model.com/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu
