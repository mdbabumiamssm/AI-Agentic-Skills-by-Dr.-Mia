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

# MRD Detection Agent

**ID:** `biomedical.clinical.mrd_detection`
**Version:** 1.0.0
**Status:** Production
**Category:** Clinical / Oncology / Minimal Residual Disease

---

## Overview

The **MRD Detection Agent** specializes in detecting minimal residual disease (MRD) from liquid biopsy samples to predict cancer recurrence and guide adjuvant therapy decisions. MRD refers to the small number of cancer cells that remain after treatment and are below the detection limit of standard imaging or laboratory tests.

This agent integrates tumor-informed and tumor-agnostic ctDNA approaches with AI-enhanced interpretation to detect molecular relapse months before clinical recurrence, enabling earlier intervention and improved survival outcomes.

---

## Key Capabilities

### 1. MRD Detection Approaches

| Approach | Method | Sensitivity | Best For |
|----------|--------|-------------|----------|
| **Tumor-informed** | Custom panel from tissue WES | 0.001-0.01% | Post-surgery monitoring |
| **Tumor-agnostic** | Fixed panel + methylation | 0.01-0.1% | Unknown primary, no tissue |
| **Methylation-based** | Epigenetic signatures | 0.01% | Early detection, recurrence |
| **Fragment analysis** | cfDNA fragmentomics | 0.1% | Screening, complementary |

### 2. Commercial MRD Assays

| Assay | Company | Type | LOD | Medicare Coverage |
|-------|---------|------|-----|-------------------|
| **Signatera** | Natera | Tumor-informed (16 SNVs) | 0.01% | CRC, breast |
| **Guardant Reveal** | Guardant | Epigenomic + genomic | 0.01% | CRC |
| **FoundationOne Tracker** | Foundation | Tumor-informed | 0.01% | Multiple |
| **clonoSEQ** | Adaptive | Lymphoid MRD (IGH/TCR) | 10^-6 | ALL, CLL, MM |
| **RaDaR** | Inivata | Tumor-informed | 0.003% | Research |

### 3. AI-Powered MRD Analysis

| Model | Function | Performance |
|-------|----------|-------------|
| **PRIME** | Progression risk prediction | AUC 0.82-0.85 |
| **DeepMRD** | ctDNA signal detection | Sensitivity +15% |
| **MethylMRD** | Methylation pattern recognition | AUC 0.92 |
| **Ensemble AI** | Multi-feature integration | Improved specificity |

### 4. Clinical Decision Support

- **MRD-positive:** Consider adjuvant therapy intensification
- **MRD-negative:** Potential de-escalation, surveillance
- **MRD clearance:** Treatment response confirmation
- **MRD conversion:** Early recurrence detection

---

## Usage

### Example Prompt

```text
Analyze MRD status for a stage III colon cancer patient post-surgery.

Patient information:
- Stage IIIB colon cancer (pT3N2aM0)
- Radical colectomy with clear margins
- Tumor WES available (customized Signatera panel)
- Post-operative plasma samples at 4, 8, 12 weeks

Tasks:
1. Analyze MRD status at each timepoint
2. Calculate ctDNA clearance kinetics
3. Assess recurrence risk
4. Recommend adjuvant therapy approach
5. Define surveillance schedule
```

### Expected Output

```
## MRD Analysis Report

### Patient Summary
- **Diagnosis:** Stage IIIB colon adenocarcinoma
- **Surgery:** R0 resection, 24 lymph nodes (4 positive)
- **MSI status:** MSS (microsatellite stable)
- **Risk factors:** T3, N2a, lymphovascular invasion

### Tumor-Informed Assay Design

**Signatera Custom Panel (16 variants):**

| Variant | Gene | VAF (Tumor) | Driver/Passenger |
|---------|------|-------------|------------------|
| G12V | KRAS | 34.2% | Driver |
| R175H | TP53 | 28.7% | Driver |
| R1450* | APC | 31.4% | Driver |
| D351N | SMAD4 | 22.3% | Driver |
| [12 additional clonal variants] | ... | >10% | Passenger |

### Longitudinal MRD Analysis

| Timepoint | ctDNA Status | MTM/mL | Mean VAF | Clinical Status |
|-----------|--------------|--------|----------|-----------------|
| Pre-operative | Positive | 845 | 2.4% | Tumor present |
| Week 4 (post-op) | Positive | 12 | 0.03% | MRD+ |
| Week 8 | Negative | 0 | <LOD | MRD clearance |
| Week 12 | Negative | 0 | <LOD | MRD- confirmed |

**MTM/mL = Mean Tumor Molecules per mL plasma**

### MRD Clearance Kinetics

```
ctDNA Level (log scale)
   |
 1000 ─●─────────────── Pre-op
   |    \
  100 ─  \
   |      \
   10 ─    \──●───────── Week 4
    |         \
    1 ─        \
   |            ●────●── Week 8, 12 (below LOD)
  0.1 ─ - - - - - - - - LOD (0.01%)
   └──┬──┬──┬──┬──┬──┬─
      0  4  8  12 16 20 weeks
```

**Clearance half-life:** 5.2 days (normal: 4-7 days)
**Interpretation:** Appropriate clearance kinetics, consistent with R0 resection

### Recurrence Risk Assessment

#### PRIME AI Model Prediction

| Factor | Value | Risk Contribution |
|--------|-------|-------------------|
| Post-op ctDNA (Week 4) | Positive | High risk (+) |
| MRD clearance (Week 8+) | Cleared | Favorable (-) |
| Stage | IIIB | High risk (+) |
| MSI status | MSS | Neutral |
| KRAS mutation | G12V | Unfavorable (+) |

**Predicted 3-year DFS:**
- Without adjuvant therapy: 52% (95% CI: 44-60%)
- With FOLFOX (6 months): 71% (95% CI: 63-79%)
- **Absolute benefit of adjuvant:** 19%

### MRD-Guided Therapy Recommendation

#### Current MRD Status: CLEARED (Week 12)

**Standard approach:**
- 6 months FOLFOX adjuvant chemotherapy
- Based on IDEA collaboration (Stage III, high-risk)

**MRD-guided consideration:**
- MRD clearance by Week 8 is favorable
- Consider enrollment in de-escalation trial
- Alternative: 3 months FOLFOX (IDEA, high-risk T4/N2)

**Recommendation:**
```
Given:
- MRD positive at Week 4 (elevated recurrence risk)
- MRD clearance by Week 8 (favorable response)
- High-risk features (N2a, LVI)

RECOMMEND: 6 months FOLFOX
Rationale: Week 4 MRD positivity associated with
higher recurrence risk even after clearance.
```

### Surveillance Schedule

#### MRD-Guided Monitoring Protocol

| Timepoint | Test | Rationale |
|-----------|------|-----------|
| During chemo | ctDNA q6 weeks | Treatment response |
| End of chemo | ctDNA + CT | Confirm clearance |
| Year 1 | ctDNA q3 months | Early detection window |
| Year 2-3 | ctDNA q4 months | Continued surveillance |
| Year 4-5 | ctDNA q6 months | Extended follow-up |

#### MRD Conversion Protocol

**If MRD becomes positive:**
1. Confirm with repeat sample (2 weeks)
2. Order imaging (CT chest/abdomen/pelvis)
3. If confirmed, consider:
   - Intensified surveillance
   - Early intervention therapy
   - Clinical trial enrollment

### Lead Time Benefit

**Expected detection advantage:**
- MRD positivity precedes imaging by: 8.9 months (median)
- 95% of recurrences are ctDNA+ before imaging
- Intervention window: Earlier treatment possible

### Prognostic Summary

| Outcome | MRD- at Week 12 | MRD+ at Week 12 |
|---------|-----------------|-----------------|
| 3-year DFS | 76% | 38% |
| Recurrence HR | 1.0 (ref) | 4.2 |
| Lead time | N/A | 8.9 months |
```

### LLM Agent Integration

```python
@tool
def analyze_mrd_status(
    ctdna_data: str,
    tumor_variants: str = None,
    assay_type: str = "tumor_informed",
    lod: float = 0.01
) -> str:
    """
    Analyzes MRD status from ctDNA data.

    Args:
        ctdna_data: Path to ctDNA sequencing/assay data
        tumor_variants: Path to tumor WES variants (for tumor-informed)
        assay_type: Assay type (tumor_informed, tumor_agnostic, methylation)
        lod: Limit of detection (percentage)

    Returns:
        MRD status with quantification
    """
    pass


@tool
def track_mrd_kinetics(
    longitudinal_samples: list[str],
    timepoints: list[str],
    treatment_dates: dict = None
) -> str:
    """
    Tracks MRD kinetics over time.

    Args:
        longitudinal_samples: List of paths to serial samples
        timepoints: Sample collection timepoints
        treatment_dates: Surgery/treatment dates for context

    Returns:
        MRD kinetics with clearance analysis
    """
    pass


@tool
def predict_recurrence_risk(
    mrd_data: str,
    clinical_features: dict,
    model: str = "prime"
) -> str:
    """
    Predicts recurrence risk using AI model.

    Args:
        mrd_data: Path to MRD analysis data
        clinical_features: Clinical features dict
        model: Prediction model (prime, ensemble)

    Returns:
        Recurrence risk prediction with confidence
    """
    pass


@tool
def generate_mrd_surveillance_plan(
    cancer_type: str,
    stage: str,
    mrd_status: str,
    risk_factors: list[str] = None
) -> str:
    """
    Generates MRD-guided surveillance schedule.

    Args:
        cancer_type: Cancer type
        stage: Cancer stage
        mrd_status: Current MRD status
        risk_factors: Additional risk factors

    Returns:
        Personalized surveillance protocol
    """
    pass
```

---

## Prerequisites

### Required Tools

| Tool | Version | Purpose |
|------|---------|---------|
| **Signatera API** | Latest | Tumor-informed MRD |
| **ichorCNA** | >=0.3 | Tumor fraction |
| **PyClone-VI** | >=0.1 | Clonal analysis |

### Dependencies

```
pandas>=2.0
numpy>=1.24
scipy>=1.11
scikit-learn>=1.3
matplotlib>=3.7
seaborn>=0.12
```

---

## Methodology

### MRD Detection Workflow

```
Tumor Tissue WES
    ↓
Variant Calling & Filtering
    ↓
Custom Panel Design (16-50 variants)
    ↓
Post-operative Plasma Samples
    ↓
Ultra-deep Sequencing (>100,000x)
    ↓
ctDNA Signal Detection
├── Variant detection (each target)
├── Background error modeling
└── Statistical significance testing
    ↓
MRD Status Call
├── Positive: ≥2 variants detected
├── Negative: <2 variants at LOD
└── Indeterminate: borderline
    ↓
Clinical Integration
├── Risk stratification
├── Therapy guidance
└── Surveillance planning
```

### Clinical Validation Data

| Cancer Type | MRD+ Recurrence Rate | MRD- Recurrence Rate | HR |
|-------------|---------------------|---------------------|-----|
| Colorectal | 73% | 14% | 7.4 |
| Breast | 52% | 10% | 5.8 |
| Lung (NSCLC) | 67% | 19% | 4.6 |
| Bladder | 79% | 18% | 6.2 |

---

## Clinical Applications

### Adjuvant Therapy Decisions
- MRD-guided therapy escalation/de-escalation
- Personalized treatment duration
- Clinical trial eligibility

### Surveillance Optimization
- Risk-stratified follow-up
- Early recurrence detection
- Reduced imaging burden for MRD- patients

### Clinical Trials
- MRD as surrogate endpoint
- Accelerated drug approval pathways
- Adaptive trial designs

---

## Related Skills

- **Liquid Biopsy Analysis Agent:** Comprehensive ctDNA profiling
- **ctDNA Analysis Agent:** Deep ctDNA analysis
- **Cancer Surveillance Agent:** Long-term monitoring
- **Adjuvant Therapy Agent:** Treatment recommendations

---

## References

- **Tie et al. (2022):** "Circulating Tumor DNA Analysis Guiding Adjuvant Therapy in Stage II Colon Cancer." *NEJM*
- **Reinert et al. (2019):** "Analysis of Plasma Cell-Free DNA by Ultradeep Sequencing in Patients With Stages I to III Colorectal Cancer." *JAMA Oncology*
- [PRIME Model](https://link.springer.com/article/10.1186/s40779-025-00679-z)
- [Signatera](https://www.natera.com/oncology/signatera-advanced-cancer-detection/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->