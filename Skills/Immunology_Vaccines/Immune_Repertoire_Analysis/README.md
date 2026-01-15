# Immune Repertoire Analysis Agent

**ID:** `biomedical.immunology.immune_repertoire_analysis`
**Version:** 1.0.0
**Status:** Production
**Category:** Immunology / Adaptive Immunity / Repertoire Analysis

---

## Overview

The **Immune Repertoire Analysis Agent** performs comprehensive analysis of T cell receptor (TCR) and B cell receptor (BCR) repertoires from high-throughput sequencing data. Immune repertoire sequencing (Rep-seq) provides deep insights into adaptive immune system diversity, clonal dynamics, and antigen-specific responses critical for understanding immunity, autoimmunity, cancer, and vaccine responses.

This agent integrates cutting-edge bioinformatics tools and machine learning approaches including MAL-ID (Machine Learning for Immunological Diagnosis) to extract clinically actionable information from immune repertoire data.

---

## Key Capabilities

### 1. Repertoire Sequencing Methods

| Method | Target | Resolution | Best For |
|--------|--------|------------|----------|
| **Bulk TCR-seq** | TCRα/β chains | Clonotype-level | Population diversity |
| **Bulk BCR-seq** | BCR heavy/light | Clonotype-level | Antibody repertoire |
| **Single-cell TCR** | Paired αβ chains | Single-cell | TCR clonality + phenotype |
| **Single-cell BCR** | Paired H/L chains | Single-cell | Antibody discovery |
| **CITE-seq + VDJ** | TCR/BCR + protein | Single-cell | Multi-modal profiling |

### 2. Analysis Capabilities

| Analysis | Description | Tools |
|----------|-------------|-------|
| **V(D)J annotation** | Gene segment assignment | IgBLAST, MiXCR, IMGT |
| **Clonotype identification** | CDR3 clustering | MiXCR, Immunarch |
| **Diversity metrics** | Shannon, Simpson, Chao1 | Immunarch, tcR |
| **Clonal expansion** | Clone size distribution | Immunarch, VDJtools |
| **Public clones** | Shared clonotypes | GLIPH2, tcrdist3 |
| **Antigen specificity** | Epitope prediction | ERGO, DeepTCR |
| **Disease diagnosis** | Multi-disease classification | MAL-ID |

### 3. MAL-ID Framework

The **Machine Learning for Immunological Diagnosis (MAL-ID)** framework enables:
- Simultaneous screening for multiple illnesses
- Detection of specific infections, autoimmune disorders
- Vaccine response assessment
- Disease severity classification

### 4. Clinical Applications

- **Infectious disease:** Pathogen exposure history
- **Autoimmunity:** Autoreactive clone detection
- **Cancer:** TIL repertoire analysis
- **Vaccine development:** Response monitoring
- **Transplant:** Rejection risk assessment

---

## Usage

### Example Prompt

```text
Analyze this bulk TCR-seq dataset from COVID-19 patients to characterize the SARS-CoV-2 specific immune response.

Dataset:
- 50 COVID-19 patients (25 mild, 25 severe)
- 25 healthy controls
- TCRβ chain sequencing
- Timepoint: 4 weeks post-infection

Tasks:
1. Perform quality control and preprocessing
2. Calculate diversity metrics
3. Identify expanded clones
4. Find SARS-CoV-2 specific TCRs
5. Compare repertoires between severity groups
6. Identify prognostic features
```

### Expected Output

```
## Immune Repertoire Analysis Report

### Dataset Summary
- **Samples:** 75 (50 COVID-19, 25 controls)
- **Sequencing:** Bulk TCRβ chain
- **Mean reads/sample:** 2.4M
- **Mean clonotypes/sample:** 124,567

### Quality Control

| Metric | Median | Range | Status |
|--------|--------|-------|--------|
| Total reads | 2.4M | 1.2-4.1M | Pass |
| In-frame sequences | 92.3% | 89-95% | Pass |
| Unique clonotypes | 124,567 | 78K-189K | Pass |
| CDR3 length distribution | Normal | - | Pass |

**QC failures:** 2 samples removed (low read count)
**Final dataset:** 73 samples

### Diversity Analysis

#### Alpha Diversity (Within-Sample)

| Group | Shannon Index | Simpson Index | Chao1 |
|-------|---------------|---------------|-------|
| Healthy | 11.2 ± 0.8 | 0.998 ± 0.001 | 156,789 |
| COVID Mild | 10.4 ± 1.1 | 0.996 ± 0.002 | 134,234 |
| COVID Severe | 8.7 ± 1.4 | 0.989 ± 0.008 | 98,456 |

**Statistical comparison:**
- Severe vs Healthy: p < 0.001 (Shannon)
- Severe vs Mild: p = 0.003 (Shannon)
- Mild vs Healthy: p = 0.024 (Shannon)

**Interpretation:** Severe COVID-19 associated with repertoire contraction

#### Clonal Expansion Analysis

| Group | Top Clone % | Gini Index | Expanded Clones (>1%) |
|-------|-------------|------------|----------------------|
| Healthy | 0.8 ± 0.4% | 0.42 | 2.1 |
| COVID Mild | 2.3 ± 1.2% | 0.58 | 8.4 |
| COVID Severe | 8.7 ± 4.2% | 0.74 | 24.6 |

**Key finding:** Severe COVID shows oligoclonal expansion

### SARS-CoV-2 Specific TCR Identification

#### Public Clonotype Analysis

**Method:** Comparison to VDJdb + COVID-19 specific databases

| Epitope | Protein | HLA | Public Clones Found | Patients |
|---------|---------|-----|---------------------|----------|
| YLQPRTFLL | Spike | A*02:01 | 234 | 42/50 |
| TTDPSFLGRY | Spike | A*03:01 | 156 | 28/50 |
| KLPDDFTGCV | Spike | A*02:01 | 189 | 38/50 |
| NYNYLYRLF | Nucleocapsid | A*24:02 | 112 | 21/50 |
| LLLDRLNQL | ORF1ab | A*02:01 | 78 | 18/50 |

**Total unique SARS-CoV-2 TCRs identified:** 1,847

#### Top Expanded COVID-Specific Clones

| Rank | CDR3β Sequence | TRBV | TRBJ | Epitope | Expansion |
|------|----------------|------|------|---------|-----------|
| 1 | CASSIRSSYEQYF | 19 | 2-7 | YLQPRTFLL | 12.4% |
| 2 | CASSLGQGYEQYF | 6-1 | 2-7 | YLQPRTFLL | 8.7% |
| 3 | CASSLAPGATNEKLFF | 6-5 | 1-4 | KLPDDFTGCV | 6.2% |
| 4 | CASSFQGTSGYTF | 12-3 | 1-2 | TTDPSFLGRY | 5.8% |
| 5 | CASSYRGTEAFF | 6-1 | 1-1 | NYNYLYRLF | 4.3% |

### Severity Group Comparison

#### Differential Features

| Feature | Mild | Severe | P-value | Direction |
|---------|------|--------|---------|-----------|
| Shannon diversity | 10.4 | 8.7 | <0.001 | Severe ↓ |
| Top clone size | 2.3% | 8.7% | <0.001 | Severe ↑ |
| COVID-specific TCRs | 3.2% | 8.9% | 0.002 | Severe ↑ |
| TRBV28 usage | 4.1% | 8.7% | 0.008 | Severe ↑ |
| Public clone overlap | 12.3% | 18.7% | 0.034 | Severe ↑ |

#### GLIPH2 Analysis: Specificity Groups

| Cluster | Size | Motif | Predicted Specificity | Association |
|---------|------|-------|----------------------|-------------|
| C1 | 456 | %SIRS% | Spike (S269) | Severe (OR=3.2) |
| C2 | 312 | %LAPG% | Spike (S1) | Both |
| C3 | 234 | %QGTS% | Nucleocapsid | Mild (OR=2.1) |
| C4 | 189 | %RGTE% | ORF1ab | Severe (OR=2.8) |

### Prognostic Features

#### Machine Learning Model (MAL-ID Approach)

**Features selected:**
1. Repertoire diversity (Shannon)
2. Top clone frequency
3. SARS-CoV-2 TCR fraction
4. V-gene usage entropy
5. CDR3 length distribution

**Model Performance:**

| Metric | Value | 95% CI |
|--------|-------|--------|
| AUC (severity prediction) | 0.87 | 0.79-0.94 |
| Sensitivity | 84% | 72-92% |
| Specificity | 81% | 68-90% |
| PPV | 79% | - |
| NPV | 86% | - |

#### Risk Stratification

| Risk Group | Criteria | Progression to Severe |
|------------|----------|----------------------|
| Low | Shannon > 10, Top clone < 3% | 8% |
| Intermediate | Shannon 9-10 OR Top clone 3-6% | 34% |
| High | Shannon < 9 AND Top clone > 6% | 72% |

### Clinical Insights

**Key Findings:**

1. **Severity biomarkers:**
   - Repertoire contraction strongly associated with severe disease
   - Oligoclonal expansion indicates immune dysregulation
   - TRBV28 bias in severe patients (ORF1ab response)

2. **Protective signatures:**
   - Preserved diversity associated with mild disease
   - Balanced response to multiple epitopes
   - Higher TRBV6-1 usage (Spike-specific)

3. **Prognostic utility:**
   - Early repertoire features predict disease trajectory
   - Potential for therapeutic intervention guidance

### Visualization Outputs

Generated files:
1. `diversity_comparison.png` - Diversity metrics by group
2. `clonal_expansion.png` - Clone size distributions
3. `vj_usage_heatmap.png` - V/J gene usage
4. `covid_tcr_network.html` - Interactive TCR similarity network
5. `severity_roc.png` - Prognostic model performance
```

### LLM Agent Integration

```python
@tool
def analyze_tcr_repertoire(
    sample_files: list[str],
    chain: str = "TRB",
    reference: str = "imgt",
    min_reads: int = 10
) -> str:
    """
    Analyzes TCR repertoire from sequencing data.

    Args:
        sample_files: Paths to TCR sequencing files
        chain: TCR chain (TRA, TRB, paired)
        reference: V(D)J reference database
        min_reads: Minimum reads for clonotype inclusion

    Returns:
        Comprehensive repertoire analysis
    """
    pass


@tool
def analyze_bcr_repertoire(
    sample_files: list[str],
    chain: str = "IGH",
    include_shm: bool = True,
    include_isotype: bool = True
) -> str:
    """
    Analyzes BCR repertoire including somatic hypermutation.

    Args:
        sample_files: Paths to BCR sequencing files
        chain: BCR chain (IGH, IGK, IGL)
        include_shm: Analyze somatic hypermutation
        include_isotype: Include isotype analysis

    Returns:
        BCR repertoire analysis with SHM and isotype
    """
    pass


@tool
def calculate_diversity_metrics(
    repertoire_data: str,
    metrics: list[str] = ["shannon", "simpson", "chao1"]
) -> str:
    """
    Calculates repertoire diversity metrics.

    Args:
        repertoire_data: Path to repertoire data
        metrics: Diversity metrics to calculate

    Returns:
        Diversity metrics with confidence intervals
    """
    pass


@tool
def find_antigen_specific_tcrs(
    repertoire_data: str,
    antigen: str = None,
    database: str = "vdjdb",
    method: str = "sequence_matching"
) -> str:
    """
    Identifies antigen-specific TCRs.

    Args:
        repertoire_data: Path to repertoire data
        antigen: Specific antigen to search (None = all)
        database: Reference database (vdjdb, mcpas, iedb)
        method: Matching method (sequence_matching, gliph2, tcrdist)

    Returns:
        Antigen-specific TCRs with annotations
    """
    pass


@tool
def compare_repertoires(
    group1_files: list[str],
    group2_files: list[str],
    comparison_type: str = "differential"
) -> str:
    """
    Compares repertoires between groups.

    Args:
        group1_files: Files for group 1
        group2_files: Files for group 2
        comparison_type: Type (differential, overlap, clustering)

    Returns:
        Differential features and statistics
    """
    pass


@tool
def diagnose_from_repertoire(
    repertoire_data: str,
    conditions: list[str] = None,
    model: str = "malid"
) -> str:
    """
    Performs disease diagnosis from repertoire using ML.

    Args:
        repertoire_data: Path to repertoire data
        conditions: Conditions to test (None = all supported)
        model: Diagnostic model (malid, custom)

    Returns:
        Disease predictions with probabilities
    """
    pass
```

---

## Prerequisites

### Required Tools

| Tool | Version | Purpose |
|------|---------|---------|
| **MiXCR** | >=4.0 | V(D)J annotation |
| **Immunarch** | >=0.9 | Repertoire analysis (R) |
| **GLIPH2** | Latest | TCR clustering |
| **tcrdist3** | >=0.2 | TCR distance calculation |

### Dependencies

```
mixcr>=4.0.0
pandas>=2.0
numpy>=1.24
scipy>=1.11
scikit-learn>=1.3
networkx>=3.0
```

---

## Methodology

### Repertoire Analysis Pipeline

```
Raw Sequencing Data (FASTQ)
    ↓
V(D)J Annotation (MiXCR/IgBLAST)
├── V gene assignment
├── D gene assignment (TCRβ/IGH)
├── J gene assignment
└── CDR3 extraction
    ↓
Clonotype Identification
├── CDR3 nucleotide clustering
├── Error correction
└── Clone size quantification
    ↓
Quality Control
├── Read count filtering
├── In-frame filtering
└── Productive chain selection
    ↓
Diversity Analysis
├── Alpha diversity (within sample)
├── Beta diversity (between samples)
└── Rarefaction analysis
    ↓
Clonal Analysis
├── Expansion detection
├── Public clone identification
└── Specificity prediction
    ↓
Clinical Integration
├── Disease association
├── Prognostic modeling
└── Biomarker discovery
```

---

## Clinical Applications

### Infectious Disease
- Pathogen exposure history from TCR/BCR
- Vaccine response monitoring
- Correlates of protection identification

### Cancer Immunology
- TIL repertoire analysis
- Immunotherapy response prediction
- Neoantigen-specific T cell detection

### Autoimmunity
- Autoreactive clone detection
- Disease-specific signatures
- Treatment response monitoring

### Transplantation
- Alloreactive clone monitoring
- Rejection risk assessment
- Immune reconstitution tracking

---

## Related Skills

- **TCR-Epitope Agent:** TCR-pMHC binding prediction
- **CAR-T Design Agent:** Engineered receptor design
- **Neoantigen Vaccine Agent:** T cell epitope selection
- **Immunarch Analysis Agent:** R-based repertoire analysis

---

## References

- **Emerson et al. (2017):** "Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA-mediated effects on the T cell repertoire." *Nature Genetics*
- **Minervina et al. (2022):** "Convergent antigen-specific T cell responses in SARS-CoV-2 infection." *Cell*
- **MAL-ID Paper (2024):** "Disease diagnostics using machine learning of B cell and T cell receptor sequences." *Science*
- [Immunarch](https://immunarch.com/)
- [VDJdb](https://vdjdb.cdr3.net/)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu
