# AML Classification Agent

**ID:** `biomedical.hematology.aml_classification`
**Version:** 1.0.0
**Status:** Production
**Category:** Hematology / Acute Myeloid Leukemia / Diagnosis

---

## Overview

The **AML Classification Agent** provides comprehensive diagnostic classification and risk stratification for acute myeloid leukemia (AML) according to the WHO 2022 and ICC 2022 classifications. AML is a highly heterogeneous malignancy requiring integration of morphology, immunophenotype, cytogenetics, and molecular genetics (MICM) for accurate diagnosis and prognosis.

This agent integrates AI-powered diagnostic tools including the Acute Leukemia Methylome Atlas (ALMA) and deep learning classifiers to accelerate diagnosis, enable precise WHO subtyping, and provide ELN 2022 risk stratification for treatment planning.

---

## Key Capabilities

### 1. WHO 2022 AML Classification

| Category | Defining Features | Examples |
|----------|-------------------|----------|
| **AML with defining genetic abnormalities** | Recurrent translocations/mutations | t(8;21), inv(16), t(15;17), NPM1 |
| **AML, myelodysplasia-related** | MDS-related changes | Prior MDS, MDS mutations |
| **AML with other defined genetic alterations** | Other recurrent changes | CEBPA, TP53, KMT2A |
| **AML, defined by differentiation** | Morphology-based | M0-M7 equivalents |

### 2. Molecular Markers

| Gene | Frequency | Prognostic Impact | Targeted Therapy |
|------|-----------|-------------------|------------------|
| **NPM1** | 30% | Favorable (if FLT3-) | - |
| **FLT3-ITD** | 25% | Adverse | Midostaurin, gilteritinib |
| **FLT3-TKD** | 7% | Neutral | Midostaurin, gilteritinib |
| **CEBPA** (biallelic) | 5% | Favorable | - |
| **IDH1** | 8% | Neutral | Ivosidenib |
| **IDH2** | 12% | Neutral | Enasidenib |
| **RUNX1** | 10% | Adverse | - |
| **ASXL1** | 10% | Adverse | - |
| **TP53** | 10% | Very adverse | - |

### 3. ELN 2022 Risk Stratification

| Risk | Genetic Features | 5-year OS |
|------|------------------|-----------|
| **Favorable** | t(8;21), inv(16), NPM1 mutated (no FLT3-ITD), bZIP CEBPA | 60-70% |
| **Intermediate** | NPM1 wild-type without adverse, t(9;11), other | 40-50% |
| **Adverse** | t(6;9), t(v;11q23.3) except t(9;11), TP53, complex karyotype, monosomal karyotype, ASXL1, RUNX1, high FLT3-ITD | 10-20% |

### 4. AI-Powered Diagnostics

| Tool | Function | Performance |
|------|----------|-------------|
| **ALMA** | Methylome-based subtyping | 27 WHO subtypes |
| **DeepFlow** | Flow cytometry classification | 95% accuracy |
| **AI karyotyping** | Cytogenetic analysis | Cryptic abnormalities |
| **ML prognostic** | Survival prediction | AUC 0.78-0.85 |

---

## Usage

### Example Prompt

```text
Classify this AML case and provide risk stratification for treatment planning.

Patient: 52-year-old male with newly diagnosed acute leukemia

Findings:
- Bone marrow: 78% blasts
- Morphology: Myeloid differentiation with Auer rods
- Flow cytometry: CD34+, CD117+, CD13+, CD33+, MPO+, HLA-DR+
- Cytogenetics: 46,XY,t(8;21)(q22;q22)[20]
- Molecular: RUNX1-RUNX1T1 fusion detected, KIT D816V mutation

Tasks:
1. WHO 2022 classification
2. ELN 2022 risk stratification
3. Treatment recommendations
4. MRD monitoring strategy
5. Prognosis assessment
```

### Expected Output

```
## AML Classification and Risk Stratification Report

### Patient Information
- **Age:** 52 years
- **Sex:** Male
- **Presentation:** De novo AML

### Diagnostic Findings Summary

| Category | Finding | Significance |
|----------|---------|--------------|
| **Morphology** | 78% blasts, Auer rods | AML confirmed |
| **Immunophenotype** | Myeloid (CD34+, CD117+, CD13+, CD33+, MPO+) | Myeloid lineage |
| **Cytogenetics** | t(8;21)(q22;q22) | WHO-defining |
| **Molecular** | RUNX1-RUNX1T1, KIT D816V | Core binding factor AML |

### WHO 2022 Classification

```
┌─────────────────────────────────────────────────────────────┐
│                                                             │
│   WHO 2022 DIAGNOSIS:                                       │
│                                                             │
│   AML WITH t(8;21)(q22;q22); RUNX1::RUNX1T1                │
│                                                             │
│   Category: AML with defining genetic abnormalities         │
│   Subtype: Core binding factor AML (CBF-AML)               │
│                                                             │
│   Note: Diagnosis established regardless of blast %         │
│   (WHO 2022 allows <20% blasts if defining abnormality)    │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

#### Classification Evidence

| Criterion | Finding | Met |
|-----------|---------|-----|
| t(8;21)(q22;q22) | Present | ✓ |
| RUNX1-RUNX1T1 fusion | Detected | ✓ |
| Myeloid blasts | 78% | ✓ |
| Auer rods | Present | ✓ |

### Molecular Profile Analysis

#### Detected Mutations

| Gene | Mutation | VAF | Clinical Impact |
|------|----------|-----|-----------------|
| **RUNX1-RUNX1T1** | Fusion | N/A | WHO-defining |
| **KIT** | D816V | 28% | Adverse modifier |

#### Not Detected (Relevant Negatives)

| Gene | Implication |
|------|-------------|
| FLT3-ITD | Favorable (absent) |
| NPM1 | Expected negative in t(8;21) |
| TP53 | Favorable (absent) |
| ASXL1 | Favorable (absent) |

### ELN 2022 Risk Stratification

```
┌─────────────────────────────────────────────────────────────┐
│                                                             │
│   ELN 2022 RISK CATEGORY:                                   │
│                                                             │
│   FAVORABLE (with adverse modifier)                         │
│                                                             │
│   Primary classification: FAVORABLE                         │
│   • t(8;21)(q22;q22); RUNX1::RUNX1T1                       │
│                                                             │
│   Modifying factor: KIT D816V mutation                      │
│   • Associated with increased relapse risk                  │
│   • Some guidelines consider INTERMEDIATE                   │
│                                                             │
│   Recommendation: Treat as favorable, close MRD monitoring  │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

#### Risk Factor Analysis

| Factor | Status | Impact |
|--------|--------|--------|
| t(8;21) | Present | Favorable |
| KIT mutation | D816V (28%) | Adverse modifier |
| FLT3-ITD | Absent | Favorable |
| TP53 | Wild-type | Favorable |
| Complex karyotype | No | Favorable |

### Prognosis Assessment

#### CBF-AML with KIT Mutation

| Metric | Without KIT | With KIT D816V | This Patient |
|--------|-------------|----------------|--------------|
| CR rate | 90-95% | 85-90% | ~88% |
| 5-year OS | 60-70% | 45-55% | ~50% |
| 5-year RFS | 55-65% | 35-45% | ~40% |
| Relapse risk | 30-40% | 45-55% | ~50% |

**Note:** KIT mutations (especially D816V) in CBF-AML are associated with higher relapse rates.

### Treatment Recommendations

#### Induction Therapy

| Regimen | Recommendation | Rationale |
|---------|----------------|-----------|
| **7+3** (standard) | First choice | Standard for favorable risk |
| **Gemtuzumab ozogamicin** | Add to induction | CD33+ AML benefit |
| **Midostaurin** | Not indicated | FLT3 negative |

**Recommended induction:**
```
Cytarabine 200 mg/m²/day CI x 7 days
+ Daunorubicin 60 mg/m²/day x 3 days
+ Gemtuzumab ozogamicin 3 mg/m² day 1, 4, 7
```

#### Consolidation Options

| Option | Indication | This Patient |
|--------|------------|--------------|
| **HiDAC x 3-4 cycles** | Standard for favorable | Recommended |
| **Allogeneic SCT in CR1** | Adverse/intermediate risk | Consider given KIT |
| **Autologous SCT** | Favorable, MRD negative | Option if MRD- |

**Recommendation:**
Given KIT D816V mutation (adverse modifier), consider:
1. HiDAC consolidation with close MRD monitoring
2. If MRD positive after consolidation → allogeneic SCT
3. If MRD negative → complete HiDAC, continue monitoring

### MRD Monitoring Strategy

#### MRD Markers for t(8;21) AML

| Method | Target | Sensitivity | Timing |
|--------|--------|-------------|--------|
| **RT-qPCR** | RUNX1-RUNX1T1 | 10^-4 to 10^-5 | Standard |
| **Flow cytometry** | Aberrant phenotype | 10^-3 to 10^-4 | Complementary |
| **NGS** | KIT D816V | 10^-3 | If applicable |

#### Monitoring Schedule

| Timepoint | Test | Threshold for Action |
|-----------|------|---------------------|
| Post-induction | BM PCR + flow | MRD ≥10^-3 = poor |
| Post-consolidation 1 | BM PCR | MRD ≥10^-3 = consider SCT |
| Post-consolidation 2 | BM PCR | Rising MRD = intervention |
| Maintenance | Blood PCR q3mo | Rising from nadir |
| Year 1-2 | BM PCR q3-4mo | Any increase |
| Year 3+ | q6mo | Relapse monitoring |

**MRD interpretation for CBF-AML:**
- >3 log reduction after induction: Favorable
- <3 log reduction: Consider SCT
- Rising MRD (≥1 log increase): Pre-emptive intervention

### Additional Recommendations

| Recommendation | Rationale |
|----------------|-----------|
| HLA typing | In case SCT needed |
| Cardiac evaluation | Anthracycline use |
| Fertility preservation | Age-appropriate |
| Tumor lysis prophylaxis | High blast count |
| CNS evaluation | Low risk, not routine |

### Summary

```
DIAGNOSTIC SUMMARY:

WHO 2022: AML with t(8;21)(q22;q22); RUNX1::RUNX1T1
ELN 2022: Favorable (with adverse KIT modifier)

TREATMENT PLAN:
1. Induction: 7+3 + gemtuzumab ozogamicin
2. MRD assessment post-induction (PCR + flow)
3. Consolidation: HiDAC cycles
4. MRD-guided decision for SCT vs continued chemo

PROGNOSIS:
- CR probability: ~88%
- 5-year OS: ~50% (reduced by KIT mutation)
- 5-year RFS: ~40%

KEY MONITORING:
- RUNX1-RUNX1T1 PCR is gold standard MRD marker
- Rising MRD = early intervention trigger
```
```

### LLM Agent Integration

```python
@tool
def classify_aml_who2022(
    cytogenetics: str,
    molecular: dict,
    morphology: dict,
    flow_cytometry: dict
) -> str:
    """
    Classifies AML according to WHO 2022.

    Args:
        cytogenetics: Karyotype result
        molecular: Molecular findings (fusions, mutations)
        morphology: Morphologic findings
        flow_cytometry: Immunophenotype

    Returns:
        WHO 2022 classification with supporting evidence
    """
    pass


@tool
def stratify_aml_eln2022(
    who_classification: str,
    molecular_profile: dict,
    cytogenetics: str
) -> str:
    """
    Performs ELN 2022 risk stratification.

    Args:
        who_classification: WHO 2022 diagnosis
        molecular_profile: All molecular findings
        cytogenetics: Karyotype

    Returns:
        ELN 2022 risk category with prognosis
    """
    pass


@tool
def recommend_aml_treatment(
    eln_risk: str,
    patient_age: int,
    comorbidities: list[str],
    molecular_targets: list[str]
) -> str:
    """
    Recommends AML treatment based on risk and targets.

    Args:
        eln_risk: ELN 2022 risk category
        patient_age: Patient age
        comorbidities: Relevant comorbidities
        molecular_targets: Targetable mutations

    Returns:
        Treatment recommendations with rationale
    """
    pass


@tool
def design_mrd_monitoring(
    aml_subtype: str,
    molecular_markers: list[str],
    treatment_phase: str
) -> str:
    """
    Designs MRD monitoring strategy for AML.

    Args:
        aml_subtype: WHO classification
        molecular_markers: Available MRD markers
        treatment_phase: Current treatment phase

    Returns:
        MRD monitoring protocol
    """
    pass
```

---

## Prerequisites

### Required Data

| Data Type | Purpose |
|-----------|---------|
| Cytogenetics | WHO classification |
| Molecular panel | Risk stratification |
| Flow cytometry | Immunophenotype |
| Morphology | Lineage confirmation |

### Dependencies

```
pandas>=2.0
numpy>=1.24
```

---

## Methodology

### AML Diagnostic Pipeline

```
Clinical Presentation (cytopenias, blasts)
    ↓
Bone Marrow Evaluation
├── Morphology + cytochemistry
├── Flow cytometry
├── Cytogenetics
└── Molecular panel
    ↓
WHO 2022 Classification
├── Defining genetic abnormalities
├── Myelodysplasia-related changes
├── Other genetic alterations
└── Differentiation-based
    ↓
ELN 2022 Risk Stratification
├── Favorable
├── Intermediate
└── Adverse
    ↓
Treatment Planning
├── Induction selection
├── Consolidation strategy
├── SCT consideration
└── Targeted therapy
    ↓
MRD Monitoring
├── Marker selection
├── Schedule design
└── Intervention triggers
```

---

## Clinical Applications

### Diagnosis
- WHO 2022 classification
- Lineage determination
- Subtype identification

### Prognosis
- ELN 2022 risk stratification
- Survival prediction
- Relapse risk assessment

### Treatment
- Induction selection
- Consolidation planning
- Targeted therapy matching
- SCT decision support

### Monitoring
- MRD marker selection
- Response assessment
- Relapse prediction

---

## Related Skills

- **Flow Cytometry AI Agent:** Immunophenotyping
- **MDS Diagnosis Agent:** AML vs MDS
- **MRD Detection Agent:** Molecular MRD
- **CHIP Analysis Agent:** Antecedent clones

---

## References

- **Döhner et al. (2022):** "Diagnosis and management of AML in adults: 2022 recommendations." *Blood*
- **Khoury et al. (2022):** "The 5th edition of the WHO Classification of Haematolymphoid Tumours: Myeloid and Histiocytic/Dendritic Neoplasms." *Leukemia*
- **ALMA (2025):** AI tool for acute leukemia methylome-based diagnosis, UF Health
- [ELN Guidelines](https://www.leukemia-net.org/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu
