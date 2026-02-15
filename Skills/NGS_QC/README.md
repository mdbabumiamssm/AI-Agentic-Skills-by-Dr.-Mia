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

# NGS Quality Control Skills

Skills for next-generation sequencing data quality assessment and preprocessing.

## Sub-categories

| Category | Description | Skills |
|----------|-------------|--------|
| **read-qc** | Quality reports and filtering | FastQC, MultiQC, adapter trimming |
| **read-alignment** | Short-read alignment | BWA, Bowtie2, STAR, HISAT2 |

## Key Tools

- **FastQC** - Quality assessment
- **MultiQC** - Report aggregation
- **Trimmomatic/fastp** - Read trimming
- **BWA/Bowtie2** - Short-read alignment
- **STAR/HISAT2** - RNA-seq alignment

## Typical Workflow

```bash
# Quality check
fastqc -o qc_reports/ *.fastq.gz

# Adapter trimming
fastp -i input_R1.fq.gz -I input_R2.fq.gz \
      -o clean_R1.fq.gz -O clean_R2.fq.gz

# Alignment
bwa mem reference.fa clean_R1.fq.gz clean_R2.fq.gz | \
    samtools sort -o aligned.bam
```

---
*Source: bioSkills collection - integrated into BioMedical Skills Library*


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->