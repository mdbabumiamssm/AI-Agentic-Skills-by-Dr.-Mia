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

# Multi-Cancer Early Detection Agent

**ID:** `biomedical.clinical.multi_cancer_early_detection`
**Version:** 1.0.0
**Status:** Production
**Category:** Clinical / Oncology / Cancer Screening

---

## Overview

The **Multi-Cancer Early Detection (MCED) Agent** orchestrates the analysis and interpretation of blood-based tests designed to detect multiple cancer types from a single blood draw. MCED tests represent a paradigm shift in cancer screening, enabling detection of cancers that lack standard screening tests and identifying cancer signal origin (CSO) for tissue-of-origin prediction.

This agent integrates methylation-based cancer detection, machine learning classification, and clinical decision support to provide comprehensive MCED analysis with appropriate clinical context.

---

## Key Capabilities

### 1. MCED Platforms

| Platform | Company | Cancers | Technology | FDA Status |
|----------|---------|---------|------------|------------|
| **Galleri** | Grail | 50+ types | cfDNA methylation | FDA Breakthrough |
| **CancerSEEK** | Exact Sciences | 8 types | cfDNA + protein | Research |
| **PanSeer** | Singlera | 5 types | cfDNA methylation | Research |
| **cfMeDIP-seq** | Academic | Multiple | Methylated cfDNA | Research |
| **DELFI** | Delfi | Multiple | Fragmentomics | Research |

### 2. Cancer Types Detectable

| Cancer Type | Stage I Sensitivity | Stage II+ Sensitivity |
|-------------|--------------------|-----------------------|
| Pancreatic | 63% | 83% |
| Ovarian | 67% | 90% |
| Liver (HCC) | 68% | 87% |
| Esophageal | 52% | 78% |
| Lung | 24% | 67% |
| Colorectal | 17% | 65% |
| Breast | 21% | 64% |
| Prostate | 14% | 52% |

### 3. Detection Methodology

| Method | Principle | Sensitivity | Specificity |
|--------|-----------|-------------|-------------|
| **Methylation patterns** | Cancer-specific CpG signatures | Moderate-High | Very High (99.5%) |
| **Fragmentomics** | cfDNA fragment length/coverage | Moderate | High |
| **Protein biomarkers** | Circulating tumor proteins | Low-Moderate | Moderate |
| **Multi-analyte** | Combined approaches | Highest | High |

### 4. Clinical Workflow Integration

- **Cancer signal detected:** Proceed to diagnostic workup
- **Cancer signal origin:** Guide diagnostic imaging/procedures
- **No signal detected:** Standard screening continues
- **Indeterminate:** Repeat testing or enhanced surveillance

---

## Usage

### Example Prompt

```text
Interpret this Galleri multi-cancer early detection test result for a 58-year-old patient.

Patient information:
- Age: 58, Male
- Risk factors: Former smoker (20 pack-years), family history of pancreatic cancer
- No current symptoms
- Standard screenings: Colonoscopy (3 years ago), PSA normal

Test result:
- Cancer signal: DETECTED
- Cancer signal origin: Pancreas/Gallbladder

Tasks:
1. Interpret the MCED result
2. Explain the cancer signal origin prediction
3. Recommend diagnostic workup
4. Provide guidance on false positive considerations
5. Suggest follow-up protocol
```

### Expected Output

```
## Multi-Cancer Early Detection Result Interpretation

### Patient Profile
- **Age/Sex:** 58-year-old male
- **Risk factors:** Former smoker (20 pack-years), family history pancreatic cancer
- **Symptoms:** None (screening context)
- **Prior screening:** Colonoscopy current, PSA normal

### MCED Test Result

```
┌─────────────────────────────────────────────────┐
│                                                 │
│           CANCER SIGNAL DETECTED                │
│                                                 │
│   Cancer Signal Origin:                         │
│   PRIMARY: Pancreas/Gallbladder                 │
│   (Confidence: 89%)                             │
│                                                 │
│   SECONDARY: Upper GI tract                     │
│   (Probability: 7%)                             │
│                                                 │
└─────────────────────────────────────────────────┘
```

### Result Interpretation

#### Cancer Signal Detection
- **Interpretation:** The test detected a methylation pattern consistent with cancer
- **False positive rate:** ~0.5% (specificity 99.5%)
- **PPV (given prevalence):** ~43% for asymptomatic screening
- **Clinical significance:** Requires diagnostic confirmation

#### Cancer Signal Origin (CSO) Analysis

| Predicted Origin | Probability | Anatomic Sites |
|-----------------|-------------|----------------|
| **Pancreas/Gallbladder** | 89% | Pancreas, gallbladder, bile ducts |
| Upper GI | 7% | Esophagus, stomach |
| Other | 4% | Lower probability sites |

**CSO Accuracy:** For pancreas/gallbladder prediction, accuracy is ~93% when cancer is confirmed.

#### Risk Assessment Integration

| Factor | Status | Impact on Interpretation |
|--------|--------|-------------------------|
| Age (58) | Appropriate for MCED | Baseline PPV applies |
| Former smoker | Elevated pancreatic risk | Increases clinical suspicion |
| Family history (pancreatic) | High-risk | Significantly increases PPV |
| Asymptomatic | Screen-detected | Earlier stage more likely |

**Adjusted PPV with risk factors:** ~58%

### Diagnostic Workup Recommendations

#### Immediate Actions (Within 2 Weeks)

```
PRIORITY: HIGH

1. Pancreatic Protocol CT (with contrast)
   - Triple-phase CT of abdomen/pelvis
   - Focus on pancreas, bile ducts, gallbladder
   - Alternative: MRI/MRCP if contrast contraindicated

2. Tumor Markers
   - CA 19-9 (pancreatic/biliary)
   - CEA (GI malignancies)
   - AFP (if liver involvement suspected)

3. Clinical Assessment
   - Complete symptom review
   - Weight history
   - Stool changes, jaundice
   - Back pain
```

#### If CT Negative or Indeterminate

| Finding | Next Step | Rationale |
|---------|-----------|-----------|
| Normal CT | EUS | Higher sensitivity for small lesions |
| Suspicious lesion | EUS-FNA | Tissue diagnosis |
| Indeterminate | MRI/MRCP + EUS | Complementary imaging |
| Dilated bile duct | ERCP vs EUS | Evaluate obstruction |

### Imaging Interpretation Guidelines

#### CT Findings to Report

| Finding | Significance | Action |
|---------|-------------|--------|
| Pancreatic mass | High suspicion | EUS-FNA for tissue |
| Pancreatic duct dilation | Concerning | Further evaluation |
| Bile duct dilation | Obstruction possible | MRCP or EUS |
| Liver lesions | Metastatic disease | Biopsy if isolated |
| Lymphadenopathy | Staging | Include in workup |
| Normal pancreas | False positive possible | Consider EUS |

### False Positive Considerations

#### Conditions That May Cause False Positives

| Condition | Mechanism | Frequency |
|-----------|-----------|-----------|
| Chronic pancreatitis | Methylation changes | Occasional |
| Cirrhosis | Liver methylation | Occasional |
| Recent infection/inflammation | Transient changes | Rare |
| Benign tumors | Abnormal methylation | Rare |
| Unknown | Non-cancer methylation | 0.3-0.5% |

#### This Patient's FP Risk Assessment
- **Risk factors for FP:** Smoking history (mild inflammation)
- **Mitigating factors:** No known chronic pancreatitis, no cirrhosis
- **Overall assessment:** False positive less likely given high CSO confidence

### Follow-Up Protocol

#### If Cancer Confirmed

```
1. Staging workup
   - CT chest (lung metastases)
   - Consider PET-CT
   - Staging laparoscopy if surgical candidate

2. Multidisciplinary tumor board

3. Treatment planning based on stage

4. Genetic counseling (family history)
```

#### If Cancer NOT Confirmed After Workup

```
1. Document all negative findings

2. Repeat MCED in 6 months
   - May have very early cancer not yet imageable
   - Recommended given persistent signal possibility

3. Enhanced surveillance protocol:
   - CA 19-9 every 3 months (year 1)
   - CT or MRI every 6 months (year 1)
   - Annual MCED thereafter

4. Re-evaluate if symptoms develop
```

### Communication Guidance

#### Key Points for Patient Discussion

1. **Cancer signal detected does NOT equal cancer diagnosis**
   - Confirmation requires imaging/biopsy
   - False positive rate is low but exists

2. **High CSO confidence is favorable for workup**
   - 89% confidence narrows diagnostic focus
   - Faster path to diagnosis if cancer present

3. **Risk factors support pursuing workup**
   - Family history increases prior probability
   - Early detection in pancreatic cancer is valuable

4. **Timeline expectations**
   - CT should occur within 2 weeks
   - Most cases resolved within 4-6 weeks
   - Anxiety management is important

### Clinical Decision Summary

| Decision Point | Recommendation | Urgency |
|----------------|---------------|---------|
| **Immediate workup?** | YES | High |
| **Primary imaging?** | Pancreatic protocol CT | Within 2 weeks |
| **Tumor markers?** | CA 19-9, CEA | With imaging |
| **If CT negative?** | EUS | Within 4 weeks |
| **If all negative?** | Repeat MCED 6 months | Required |
```

### LLM Agent Integration

```python
@tool
def interpret_mced_result(
    test_result: dict,
    patient_info: dict,
    test_platform: str = "galleri"
) -> str:
    """
    Interprets MCED test results with clinical context.

    Args:
        test_result: Dict with signal_detected, cso, confidence
        patient_info: Patient demographics and risk factors
        test_platform: MCED platform used

    Returns:
        Clinical interpretation with PPV adjustment
    """
    pass


@tool
def recommend_diagnostic_workup(
    mced_result: dict,
    cancer_signal_origin: str,
    patient_factors: dict
) -> str:
    """
    Recommends diagnostic workup based on MCED result.

    Args:
        mced_result: MCED test result
        cancer_signal_origin: Predicted tissue of origin
        patient_factors: Patient-specific factors

    Returns:
        Prioritized diagnostic workup protocol
    """
    pass


@tool
def calculate_adjusted_ppv(
    baseline_ppv: float,
    risk_factors: list[str],
    cancer_type: str
) -> str:
    """
    Calculates risk-adjusted PPV for MCED result.

    Args:
        baseline_ppv: Test baseline PPV
        risk_factors: Patient risk factors
        cancer_type: Predicted cancer type

    Returns:
        Adjusted PPV with confidence interval
    """
    pass


@tool
def generate_followup_protocol(
    workup_result: str,
    mced_result: dict,
    cancer_confirmed: bool
) -> str:
    """
    Generates follow-up protocol based on workup outcome.

    Args:
        workup_result: Outcome of diagnostic workup
        mced_result: Original MCED result
        cancer_confirmed: Whether cancer was confirmed

    Returns:
        Follow-up surveillance protocol
    """
    pass
```

---

## Prerequisites

### Clinical Requirements

| Requirement | Purpose |
|-------------|---------|
| MCED test result | Signal/CSO data |
| Patient demographics | Risk stratification |
| Medical history | PPV adjustment |
| Imaging access | Diagnostic workup |

### Dependencies

```
pandas>=2.0
numpy>=1.24
scipy>=1.11
```

---

## Methodology

### MCED Clinical Workflow

```
Blood Draw (asymptomatic patient)
    ↓
cfDNA Extraction & Methylation Analysis
    ↓
Machine Learning Classification
├── Cancer vs No Cancer
└── Cancer Signal Origin
    ↓
Result Reporting
├── Signal Detected / Not Detected
└── CSO Prediction (if positive)
    ↓
Clinical Integration
├── Risk factor adjustment
├── PPV calculation
└── Workup prioritization
    ↓
Diagnostic Workup
├── Imaging (CT/MRI/EUS)
├── Biomarkers
└── Tissue diagnosis
    ↓
Resolution
├── Cancer confirmed → Treatment
└── Cancer ruled out → Surveillance
```

### Performance Characteristics (Galleri)

| Metric | Value | Notes |
|--------|-------|-------|
| Overall sensitivity | 51.5% | All stages |
| Stage I sensitivity | 16.8% | Early stage |
| Stage II sensitivity | 40.4% | Localized |
| Stage III sensitivity | 77.0% | Regional |
| Stage IV sensitivity | 90.1% | Metastatic |
| Specificity | 99.5% | Low false positives |
| CSO accuracy | 88.7% | When cancer confirmed |

---

## Clinical Applications

### Cancer Screening
- Complement standard screening (colonoscopy, mammography, low-dose CT)
- Screen for cancers without available tests (pancreatic, ovarian)
- High-risk population monitoring

### Early Detection
- Stage shift opportunity
- Improved survival for early-detected cancers
- Reduced cancer mortality (long-term goal)

### Clinical Trials
- PATHFINDER study outcomes
- NHS-Galleri trial (UK)
- MCED screening implementation studies

---

## Related Skills

- **Liquid Biopsy Analysis Agent:** Comprehensive cfDNA analysis
- **Cancer Risk Assessment Agent:** Risk stratification
- **Diagnostic Workup Agent:** Imaging and biopsy guidance
- **Cancer Screening Agent:** Standard screening integration

---

## References

- **Klein et al. (2021):** "Clinical validation of a targeted methylation-based multi-cancer early detection test." *Annals of Oncology*
- **Schrag et al. (2023):** "Blood-based tests for multicancer early detection (PATHFINDER)." *Lancet*
- [Grail Galleri](https://www.galleri.com/)
- [MCED Screening Guidelines](https://www.nccn.org/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->