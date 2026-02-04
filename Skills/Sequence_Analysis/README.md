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
