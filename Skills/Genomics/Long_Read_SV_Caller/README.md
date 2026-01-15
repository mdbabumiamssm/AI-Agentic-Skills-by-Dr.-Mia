# Long-Read Structural Variant Caller Agent

**ID:** `biomedical.genomics.long_read_sv_caller`
**Version:** 1.0.0
**Status:** Production
**Category:** Genomics / Structural Variants

---

## Overview

The **Long-Read Structural Variant Caller Agent** orchestrates the detection and analysis of structural variants (SVs) from long-read sequencing data (PacBio HiFi, Oxford Nanopore). Structural variants, including deletions, insertions, inversions, duplications, and translocations, are critical in understanding human disease but historically difficult to detect with short-read sequencing.

This agent integrates state-of-the-art SV calling tools (Sniffles2, SVIM, cuteSV, pbsv) with intelligent ensemble methods and AI-enhanced filtering to deliver high-confidence SV calls for clinical and research applications.

---

## Key Capabilities

### 1. Multi-Tool SV Calling Pipeline

| Tool | Best For | Technology | Key Feature |
|------|----------|------------|-------------|
| **Sniffles2** | All SV types | ONT/HiFi | Repeat-aware clustering, 11.8x faster |
| **SVIM** | High sensitivity | ONT | Probabilistic model, excels in repeats |
| **cuteSV** | Balanced performance | ONT/HiFi | Signature-based clustering |
| **pbsv** | HiFi data | PacBio | Optimized for HiFi accuracy |
| **Severus** | Somatic SVs | ONT/HiFi | Cancer-focused detection |
| **DeBreak** | Complex SVs | Long-read | De novo assembly approach |

### 2. Structural Variant Types Detected

| SV Type | Size Range | Clinical Relevance |
|---------|------------|-------------------|
| **Deletions (DEL)** | 50bp - 5Mb | Gene disruption, CNV syndromes |
| **Insertions (INS)** | 50bp - 10kb | Mobile elements, disease genes |
| **Inversions (INV)** | 1kb - 100Mb | Gene disruption, chromosomal disorders |
| **Duplications (DUP)** | 1kb - 50Mb | Copy number changes, drug response |
| **Translocations (BND)** | N/A | Cancer fusions, balanced rearrangements |

### 3. Ensemble SV Merging

- **combiSV:** Combines outputs from multiple callers into superior call set
- **NextSV:** Runs cuteSV + Sniffles with multiple aligners (minimap2, NGMLR)
- **SURVIVOR:** Merges and compares SV callsets with customizable distance thresholds
- **Jasmine:** Graph-based SV merging with population-level support

### 4. AI-Enhanced Filtering

- Deep learning-based false positive filtering
- MAMnet integration for enhanced accuracy
- Coverage and read support adaptive thresholds
- Quality score recalibration

---

## Usage

### Example Prompt

```text
Call structural variants from this Oxford Nanopore whole-genome sequencing dataset.

Input:
- Aligned BAM file: patient_001.sorted.bam
- Reference genome: GRCh38
- Sequencing depth: ~30x ONT R10.4.1

Please:
1. Run Sniffles2 and SVIM in parallel
2. Merge calls using SURVIVOR (max distance 500bp)
3. Annotate SVs with gene overlap and clinical significance
4. Generate a filtered high-confidence callset (PASS only)
5. Identify any SVs overlapping known disease genes
```

### Expected Output

```
## Structural Variant Analysis Report

### Sample Information
- **Sample:** patient_001
- **Technology:** Oxford Nanopore R10.4.1
- **Coverage:** 32.4x mean
- **Read N50:** 18.2 kb

### SV Calling Summary

| Caller | DEL | INS | INV | DUP | BND | Total |
|--------|-----|-----|-----|-----|-----|-------|
| Sniffles2 | 4,521 | 3,892 | 156 | 312 | 89 | 8,970 |
| SVIM | 5,234 | 4,156 | 178 | 398 | 112 | 10,078 |
| **Merged (SURVIVOR)** | 4,892 | 4,012 | 145 | 334 | 78 | 9,461 |
| **High-confidence** | 4,156 | 3,567 | 128 | 289 | 62 | 8,202 |

### Size Distribution

| Size Range | Deletions | Insertions |
|------------|-----------|------------|
| 50-500 bp | 2,341 | 1,892 |
| 500bp-5kb | 1,456 | 1,423 |
| 5-50 kb | 312 | 234 |
| >50 kb | 47 | 18 |

### Clinical Annotations

#### High-Priority SVs (Disease Genes)

| Rank | Type | Position | Size | Gene | Clinical |
|------|------|----------|------|------|----------|
| 1 | DEL | chr7:117,120,017-117,308,719 | 188.7kb | CFTR | Pathogenic CF allele |
| 2 | DUP | chr17:41,196,312-41,277,500 | 81.2kb | BRCA1 | VUS, exon duplication |
| 3 | INV | chr9:133,589,268-133,763,062 | 173.8kb | ABL1 | BCR-ABL breakpoint region |

### Mobile Element Insertions

| Family | Count | Median Size |
|--------|-------|-------------|
| Alu | 423 | 312 bp |
| LINE-1 | 156 | 1,234 bp |
| SVA | 34 | 1,890 bp |
| HERV | 12 | 2,456 bp |

### Quality Metrics

| Metric | Value |
|--------|-------|
| Precision (benchmarked) | 94.2% |
| Recall (benchmarked) | 91.8% |
| F1 Score | 93.0% |
| Ti/Tv ratio | 0.98 |
```

### LLM Agent Integration

```python
@tool
def call_structural_variants(
    bam_file: str,
    reference: str = "GRCh38",
    callers: list[str] = ["sniffles2", "svim"],
    min_support: int = 3,
    min_sv_length: int = 50,
    ensemble_method: str = "survivor"
) -> str:
    """
    Calls structural variants from long-read sequencing data.

    Args:
        bam_file: Path to aligned BAM file
        reference: Reference genome (GRCh38, GRCh37, mm10)
        callers: SV callers to run (sniffles2, svim, cutesv, pbsv)
        min_support: Minimum read support for SV calls
        min_sv_length: Minimum SV length in bp
        ensemble_method: Merging method (survivor, jasmine, combisv)

    Returns:
        SV callset with annotations and quality metrics
    """
    pass


@tool
def annotate_structural_variants(
    vcf_file: str,
    gene_annotations: str = "gencode_v44",
    include_clinical: bool = True,
    include_population: bool = True
) -> str:
    """
    Annotates structural variants with gene and clinical information.

    Args:
        vcf_file: Path to SV VCF file
        gene_annotations: Gene annotation source
        include_clinical: Add ClinVar/ClinGen annotations
        include_population: Add gnomAD-SV frequencies

    Returns:
        Annotated VCF with clinical interpretations
    """
    pass


@tool
def benchmark_sv_calls(
    query_vcf: str,
    truth_vcf: str,
    reference: str,
    sv_types: list[str] = ["DEL", "INS"]
) -> str:
    """
    Benchmarks SV calls against a truth set.

    Args:
        query_vcf: VCF with SV calls to evaluate
        truth_vcf: Truth set VCF (e.g., GIAB)
        reference: Reference genome
        sv_types: SV types to evaluate

    Returns:
        Precision, recall, F1 metrics by SV type
    """
    pass
```

---

## Prerequisites

### Required Tools

| Tool | Version | Purpose | Installation |
|------|---------|---------|--------------|
| **Sniffles2** | >=2.2 | Primary SV caller | `conda install sniffles` |
| **SVIM** | >=2.0 | High-sensitivity calling | `pip install svim` |
| **cuteSV** | >=2.1 | Balanced performance | `pip install cuteSV` |
| **SURVIVOR** | >=1.0.7 | VCF merging | `conda install survivor` |
| **minimap2** | >=2.26 | Read alignment | `conda install minimap2` |
| **samtools** | >=1.17 | BAM processing | `conda install samtools` |

### Dependencies

```
sniffles>=2.2.0
svim>=2.0.0
cutesv>=2.1.0
pysam>=0.21
pandas>=2.0
numpy>=1.24
cyvcf2>=0.30
```

---

## Methodology

### SV Calling Pipeline

```
Long-read FASTQ
    ↓
Alignment (minimap2/NGMLR)
    ↓
BAM Processing (samtools sort/index)
    ↓
┌─────────────┬─────────────┬─────────────┐
│  Sniffles2  │    SVIM     │   cuteSV    │
└─────────────┴─────────────┴─────────────┘
    ↓               ↓               ↓
    └───────────────┴───────────────┘
                    ↓
            Ensemble Merging (SURVIVOR)
                    ↓
            Quality Filtering
                    ↓
            Annotation (AnnotSV, VEP)
                    ↓
            Clinical Interpretation
```

### Benchmarking Performance (2025 Studies)

| Caller | F1 (DEL) | F1 (INS) | Speed | Memory |
|--------|----------|----------|-------|--------|
| Sniffles2 | 90% | 88% | 11.8x faster | Low |
| SVIM | 48% | 52% | Moderate | Low |
| cuteSV | 87% | 85% | Fast | Low |
| pbsv | 89% | 86% | Moderate | Medium |

### Recommended Strategy

1. **Clinical/High-precision:** Sniffles2 + cuteSV ensemble
2. **High-sensitivity:** Add SVIM for maximum recall
3. **Somatic/Cancer:** Include Severus for somatic detection
4. **Complex SVs:** Add DeBreak for assembly-based calling

---

## Clinical Applications

### Rare Disease Diagnosis
- Detection of pathogenic SVs missed by short-read sequencing
- Resolution of complex rearrangements in repeat regions
- Identification of disease-causing mobile element insertions

### Cancer Genomics
- Detection of somatic structural rearrangements
- Identification of gene fusions (BCR-ABL, EML4-ALK)
- Copy number analysis from long-read data

### Pharmacogenomics
- Resolution of CYP2D6 structural variants
- SMN1/SMN2 copy number determination
- HLA typing from structural information

---

## Related Skills

- **Variant Annotation Agent:** Functional annotation of SVs
- **Copy Number Analysis Agent:** CNV detection and analysis
- **Cancer Fusion Detection:** Gene fusion identification
- **Clinical Variant Interpretation:** ACMG classification

---

## References

- **Smolka et al. (2024):** "Detection of mosaic and population-level structural variants with Sniffles2." *Nature Biotechnology*
- **Sedlazeck et al. (2018):** "Accurate detection of complex structural variations using single-molecule sequencing." *Nature Methods*
- [Sniffles2 GitHub](https://github.com/fritzsedlazeck/Sniffles)
- [GIAB SV Benchmark](https://www.nist.gov/programs-projects/genome-bottle)

---

## Author

**MD BABU MIA**
*Artificial Intelligence Group*
*Icahn School of Medicine at Mount Sinai*
md.babu.mia@mssm.edu
