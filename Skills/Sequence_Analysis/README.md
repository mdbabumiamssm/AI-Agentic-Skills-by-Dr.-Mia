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

# Sequence Analysis Skills

Core skills for biological sequence manipulation, alignment, and format handling.

## Sub-categories

| Category | Description | Skills |
|----------|-------------|--------|
| **alignment** | Pairwise and multiple sequence alignment | Biopython aligners, MSA methods |
| **alignment-files** | SAM/BAM/CRAM operations | Sorting, filtering, validation, indexing |
| **sequence-io** | FASTA/FASTQ/GenBank formats | Reading, writing, format conversion |
| **sequence-manipulation** | Sequence operations | Transcription, translation, motif search |
| **restriction-analysis** | Restriction enzymes | Site finding, mapping, enzyme selection |
| **primer-design** | PCR/qPCR primer design | Primer3 integration, probe design |

## Key Tools

- **Biopython** - Core sequence manipulation
- **pysam** - SAM/BAM file handling
- **samtools/bcftools** - CLI alignment operations
- **Primer3** - Primer design

## Usage

```python
from Bio import SeqIO
from Bio.Seq import Seq

# Read sequences
for record in SeqIO.parse("input.fasta", "fasta"):
    print(record.id, len(record))
```

---
*Source: bioSkills collection - integrated into BioMedical Skills Library*


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->