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

---
name: bio-tcr-bcr-analysis-mixcr-analysis
description: Perform V(D)J alignment and clonotype assembly from TCR-seq or BCR-seq data using MiXCR. Use when processing raw immune repertoire sequencing data to identify clonotypes and their frequencies.
tool_type: cli
primary_tool: MiXCR
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# MiXCR Analysis

## Complete Workflow (Recommended)

```bash
mixcr analyze generic-tcr-amplicon \
    --species human \
    --rna \
    --rigid-left-alignment-boundary \
    --floating-right-alignment-boundary C \
    input_R1.fastq.gz input_R2.fastq.gz \
    output_prefix

mixcr analyze 10x-vdj-tcr \
    input_R1.fastq.gz input_R2.fastq.gz \
    output_prefix
```

## Step-by-Step Workflow

### Step 1: Align Reads

```bash
mixcr align \
    --species human \
    --preset generic-tcr-amplicon-umi \
    input_R1.fastq.gz input_R2.fastq.gz \
    alignments.vdjca

mixcr align \
    --species human \
    --rna \
    -OallowPartialAlignments=true \
    input_R1.fastq.gz input_R2.fastq.gz \
    alignments.vdjca
```

### Step 2: Refine and Assemble

```bash
mixcr refineTagsAndSort alignments.vdjca alignments_refined.vdjca

mixcr assemble alignments_refined.vdjca clones.clns
```

### Step 3: Export Results

```bash
mixcr exportClones \
    --chains TRB \
    --preset full \
    clones.clns \
    clones.tsv

mixcr exportClones \
    --chains TRB \
    -cloneId -readCount -readFraction \
    -nFeature CDR3 -aaFeature CDR3 \
    -vGene -dGene -jGene \
    clones.clns \
    clones_custom.tsv
```

## Preset Protocols

| Protocol | Use Case |
|----------|----------|
| `generic-tcr-amplicon` | TCR amplicon sequencing |
| `generic-bcr-amplicon` | BCR amplicon sequencing |
| `generic-tcr-amplicon-umi` | TCR amplicon with UMIs |
| `rnaseq-tcr` | TCR extraction from bulk RNA-seq |
| `rnaseq-bcr` | BCR extraction from bulk RNA-seq |
| `10x-vdj-tcr` | 10x Genomics TCR enrichment |
| `10x-vdj-bcr` | 10x Genomics BCR enrichment |
| `takara-human-tcr-v2` | Takara SMARTer kit |

## Species Support

```bash
mixcr align --species human ...
mixcr align --species mmu ...

# Available: human, mmu, rat, rhesus, dog, pig, rabbit, chicken
```

## Output Format

| Column | Description |
|--------|-------------|
| cloneId | Unique clone identifier |
| readCount | Number of reads |
| cloneFraction | Proportion of repertoire |
| nSeqCDR3 | Nucleotide CDR3 sequence |
| aaSeqCDR3 | Amino acid CDR3 sequence |
| allVHitsWithScore | V gene assignments |
| allDHitsWithScore | D gene assignments |
| allJHitsWithScore | J gene assignments |

## Quality Metrics

```bash
mixcr exportReports alignments.vdjca

# Key metrics:
# - Successfully aligned reads (>80% is good)
# - CDR3 found (>70% of aligned)
# - Clonotype count (varies by sample type)
```

## Parse MiXCR Output in Python

```python
import pandas as pd

def load_mixcr_clones(filepath):
    df = pd.read_csv(filepath, sep='\t')
    df = df.rename(columns={
        'readCount': 'count',
        'cloneFraction': 'frequency',
        'aaSeqCDR3': 'cdr3_aa',
        'nSeqCDR3': 'cdr3_nt'
    })
    return df
```

## Related Skills

- vdjtools-analysis - Downstream diversity analysis
- scirpy-analysis - Single-cell VDJ integration
- repertoire-visualization - Visualize MiXCR output


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->