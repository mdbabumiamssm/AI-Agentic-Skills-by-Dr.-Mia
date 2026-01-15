# Liquid Biopsy Analysis Agent

**ID:** `biomedical.clinical.liquid_biopsy_analysis`
**Version:** 1.0.0
**Status:** Production
**Category:** Clinical / Oncology / Liquid Biopsy

---

## Overview

The **Liquid Biopsy Analysis Agent** orchestrates the analysis of circulating biomarkers from blood samples for cancer detection, monitoring, and treatment selection. Liquid biopsy has emerged as a minimally invasive alternative to tissue biopsy, enabling real-time monitoring of tumor dynamics through analysis of circulating tumor DNA (ctDNA), circulating tumor cells (CTCs), extracellular vesicles (EVs), and cell-free RNA (cfRNA).

This agent integrates AI-enhanced analysis pipelines for multi-analyte liquid biopsy, combining genomic, epigenomic, and proteomic data for comprehensive cancer profiling with up to 95% accuracy in detecting mutations.

---

## Key Capabilities

### 1. Liquid Biopsy Analytes

| Analyte | Source | Information Content | Clinical Use |
|---------|--------|---------------------|--------------|
| **ctDNA** | Cell-free DNA | Mutations, CNVs, methylation | Mutation profiling, MRD |
| **CTCs** | Whole cells | Morphology, expression, WGS | Enumeration, profiling |
| **cfRNA** | Cell-free RNA | Gene expression, fusions | Pathway activity |
| **EVs/Exosomes** | Vesicles | miRNA, proteins, metabolites | TME insights |
| **TEPs** | Platelets | RNA from tumor education | Cancer detection |
| **Proteins** | Plasma | Tumor markers | Early detection |

### 2. Analysis Platforms

| Platform | Analyte | Application | Sensitivity |
|----------|---------|-------------|-------------|
| **Guardant360** | ctDNA | Comprehensive genomic profiling | 0.1% VAF |
| **FoundationOne Liquid** | ctDNA | Pan-cancer mutations | 0.5% VAF |
| **Guardant Reveal** | ctDNA + methylation | MRD detection | 0.01% |
| **Grail Galleri** | Methylation | Multi-cancer early detection | Stage I: 17% |
| **CellSearch** | CTCs | CTC enumeration | 1 CTC/7.5mL |
| **Signatera** | ctDNA (tumor-informed) | MRD, recurrence | 0.01% |

### 3. AI-Powered Analysis

| Method | Application | Performance |
|--------|-------------|-------------|
| **Deep learning** | Variant calling | 90% sensitivity |
| **Machine learning** | Cancer classification | 95% accuracy (7 cancers) |
| **Neural networks** | Methylation patterns | AUC 0.85-0.92 |
| **XAI** | Clinical interpretation | Transparent predictions |

### 4. Clinical Applications

- **Cancer screening:** Multi-cancer early detection (MCED)
- **Diagnosis:** Tissue-of-origin identification
- **Treatment selection:** Actionable mutation profiling
- **Monitoring:** Treatment response assessment
- **MRD detection:** Minimal residual disease
- **Resistance:** Acquired resistance mutations

---

## Usage

### Example Prompt

```text
Analyze this liquid biopsy sample from a metastatic colorectal cancer patient.

Sample information:
- Plasma cfDNA extracted
- Sequencing: Guardant360 CDx panel (74 genes)
- Prior treatment: FOLFOX + bevacizumab (progressed)
- Current: Considering second-line therapy

Tasks:
1. Analyze mutation profile from ctDNA
2. Calculate tumor fraction
3. Identify actionable alterations
4. Recommend targeted therapy options
5. Provide longitudinal monitoring strategy
```

### Expected Output

```
## Liquid Biopsy Analysis Report

### Patient Information
- **Sample type:** Plasma (10mL blood draw)
- **Cancer type:** Metastatic colorectal cancer (mCRC)
- **Prior therapy:** FOLFOX + bevacizumab
- **Clinical status:** Disease progression

### Sample Quality Metrics

| Metric | Value | Status |
|--------|-------|--------|
| cfDNA yield | 28.4 ng/mL | Adequate |
| cfDNA fragment size | 167 bp (peak) | Normal |
| Sequencing depth | 15,000x mean | Excellent |
| Unique molecules | 8.2M | Sufficient |

### Mutation Profile

#### Somatic Mutations Detected

| Gene | Mutation | VAF | Classification | Actionability |
|------|----------|-----|----------------|---------------|
| **KRAS** | G12D | 12.4% | Pathogenic | Predictive (anti-EGFR resistance) |
| **APC** | R1450* | 8.7% | Pathogenic | Diagnostic |
| **TP53** | R175H | 10.2% | Pathogenic | Prognostic |
| **PIK3CA** | E545K | 3.2% | Pathogenic | Targetable (PI3K inhibitors) |
| **SMAD4** | Q311* | 5.6% | Pathogenic | Prognostic |
| **BRAF** | Not detected | - | - | - |
| **NRAS** | Not detected | - | - | - |

#### Estimated Tumor Fraction
- **Maximum VAF:** 12.4% (KRAS G12D)
- **Estimated tumor fraction:** ~25% (ichorCNA)
- **Interpretation:** High tumor burden

#### Copy Number Alterations

| Gene | Status | Log2 Ratio | Clinical Significance |
|------|--------|------------|----------------------|
| ERBB2 (HER2) | Amplified | +2.8 | Targetable |
| MYC | Amplified | +1.9 | Prognostic |
| CDKN2A | Loss | -1.2 | Prognostic |

### Actionable Alterations & Therapy Options

#### Primary Findings

| Alteration | Drug/Class | Evidence | FDA Status |
|------------|------------|----------|------------|
| **HER2 amplification** | Trastuzumab + pertuzumab | DESTINY-CRC01 | Off-label |
| **PIK3CA E545K** | Alpelisib | Emerging | Not approved (mCRC) |
| **KRAS G12D** | Sotorasib (G12C only) | Not applicable | - |

#### Treatment Recommendations

**First-line recommendation:**
- **FOLFIRI + aflibercept** (VELOUR trial)
- Rationale: KRAS mutant, progressed on oxaliplatin

**Alternative option:**
- **HER2-directed therapy** (trastuzumab + pertuzumab + chemotherapy)
- Rationale: HER2 amplification detected
- Evidence: HERACLES, MyPathway trials
- Response rate: ~30% in HER2+ mCRC

**Clinical trial consideration:**
- KRAS G12D targeted therapies (emerging)
- Current trials: NCT04699188 (MRTX1133)

### Contraindicated Therapies

| Therapy | Reason | Evidence |
|---------|--------|----------|
| Cetuximab | KRAS G12D mutation | Level 1A |
| Panitumumab | KRAS G12D mutation | Level 1A |

### Longitudinal Monitoring Strategy

#### Recommended Schedule
| Timepoint | Analysis | Purpose |
|-----------|----------|---------|
| Baseline | Full panel (current) | Comprehensive profiling |
| 4 weeks | ctDNA quantification | Early response |
| 8 weeks | Limited panel | Response assessment |
| 12 weeks | Full panel | Resistance screening |
| Progression | Full panel + methylation | Mechanism identification |

#### Response Assessment Criteria
- **Molecular response:** >50% reduction in maximum VAF
- **Stable disease:** <25% change in VAF
- **Molecular progression:** >25% increase or new mutations

### Resistance Monitoring

**High-risk acquired resistance mutations to monitor:**
- KRAS amplification
- NRAS mutations (Q61)
- BRAF V600E
- EGFR S492R (if anti-EGFR used)
- MET amplification
- PIK3CA additional mutations
```

### LLM Agent Integration

```python
@tool
def analyze_liquid_biopsy(
    sample_data: str,
    panel: str = "comprehensive",
    include_cnv: bool = True,
    include_methylation: bool = False
) -> str:
    """
    Analyzes liquid biopsy data for mutation profiling.

    Args:
        sample_data: Path to sequencing data or VCF
        panel: Analysis panel (comprehensive, targeted, mrd)
        include_cnv: Include copy number analysis
        include_methylation: Include methylation analysis

    Returns:
        Mutation profile with clinical annotations
    """
    pass


@tool
def estimate_tumor_fraction(
    sample_data: str,
    method: str = "ichorcna"
) -> str:
    """
    Estimates tumor fraction from cfDNA data.

    Args:
        sample_data: Path to cfDNA sequencing data
        method: Estimation method (ichorcna, maximum_vaf, combined)

    Returns:
        Tumor fraction estimate with confidence interval
    """
    pass


@tool
def identify_actionable_alterations(
    mutation_profile: str,
    cancer_type: str,
    therapy_databases: list[str] = ["oncokb", "civic"]
) -> str:
    """
    Identifies actionable alterations for therapy selection.

    Args:
        mutation_profile: Path to mutation profile
        cancer_type: Cancer type for context
        therapy_databases: Databases for therapy matching

    Returns:
        Actionable alterations with therapy recommendations
    """
    pass


@tool
def monitor_treatment_response(
    baseline_sample: str,
    followup_sample: str,
    tracked_mutations: list[str] = None
) -> str:
    """
    Monitors treatment response via ctDNA dynamics.

    Args:
        baseline_sample: Baseline liquid biopsy data
        followup_sample: Follow-up sample data
        tracked_mutations: Specific mutations to track

    Returns:
        Response assessment with molecular dynamics
    """
    pass
```

---

## Prerequisites

### Required Tools

| Tool | Version | Purpose |
|------|---------|---------|
| **ichorCNA** | >=0.3 | Tumor fraction estimation |
| **GATK** | >=4.4 | Variant calling |
| **VarDict** | >=1.8 | Low-VAF variant detection |
| **OncoKB** | API | Therapy annotation |

### Dependencies

```
ichor>=0.3.0
gatk>=4.4.0
pandas>=2.0
numpy>=1.24
pysam>=0.21
cyvcf2>=0.30
requests>=2.28
```

---

## Methodology

### Liquid Biopsy Pipeline

```
Blood Sample (10-20 mL)
    ↓
Plasma Separation (within 4 hours)
    ↓
cfDNA Extraction (Streck tubes)
    ↓
Library Preparation
├── Hybrid capture (targeted)
└── WGS (shallow/deep)
    ↓
Sequencing (>10,000x for ctDNA)
    ↓
Bioinformatics Analysis
├── Variant calling (VarDict, GATK Mutect2)
├── CNV analysis (ichorCNA)
├── Methylation (if applicable)
└── Tumor fraction estimation
    ↓
Clinical Annotation
├── OncoKB
├── CIViC
└── Therapy matching
    ↓
Clinical Report
```

### Performance Metrics (2025)

| Application | Sensitivity | Specificity | PPV |
|-------------|-------------|-------------|-----|
| Mutation detection (>0.5% VAF) | 95% | 99% | 97% |
| Mutation detection (>0.1% VAF) | 85% | 98% | 92% |
| CNV detection | 90% | 95% | 88% |
| MRD detection | 88% | 95% | 90% |

---

## Clinical Applications

### Cancer Detection
- Multi-cancer early detection (MCED)
- Screening high-risk populations
- Tissue-of-origin identification

### Treatment Selection
- Comprehensive genomic profiling
- Actionable mutation identification
- Resistance mechanism detection

### Monitoring
- Treatment response assessment
- MRD detection
- Recurrence surveillance

---

## Related Skills

- **ctDNA Analysis Agent:** Deep ctDNA analysis
- **MRD Detection Agent:** Minimal residual disease
- **Tumor Profiling Agent:** Comprehensive tumor analysis
- **Therapy Matching Agent:** Treatment recommendations

---

## References

- **Corcoran & Chabner (2018):** "Application of Cell-free DNA Analysis to Cancer Treatment." *NEJM*
- **Keller et al. (2021):** "Clinical relevance of blood-based ctDNA analysis." *Clinical Cancer Research*
- [OncoKB](https://www.oncokb.org/)
- [Guardant Health](https://www.guardanthealth.com/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu
