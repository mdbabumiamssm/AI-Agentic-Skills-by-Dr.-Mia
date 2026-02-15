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

# alignment

## Overview

Sequence alignment using Biopython's Bio.Align and Bio.AlignIO modules for pairwise and multiple sequence alignments. Distinct from alignment-files which handles read-to-reference alignments (SAM/BAM).

**Tool type:** python | **Primary tools:** Bio.Align, Bio.AlignIO

## Skills

| Skill | Description |
|-------|-------------|
| pairwise-alignment | Global/local alignment using PairwiseAligner (Needleman-Wunsch, Smith-Waterman) |
| alignment-io | Read, write, convert MSA files (Clustal, PHYLIP, Stockholm, FASTA) |
| msa-parsing | Parse and analyze MSA content: gaps, conservation, filtering, consensus |
| msa-statistics | Calculate identity, conservation scores, entropy, substitution patterns |

## Example Prompts

- "Align these two DNA sequences and show the result"
- "Compare this protein to the reference using BLOSUM62"
- "Find the best matching region between these sequences"
- "What is the alignment score between seq1 and seq2?"
- "Read this Clustal alignment and show sequence IDs"
- "Convert my PHYLIP alignment to FASTA format"
- "Extract columns 100-200 from the alignment"
- "Save this alignment as Stockholm format"
- "Find conserved positions in this alignment"
- "Remove columns with more than 50% gaps"
- "Generate a consensus sequence"
- "Filter out sequences with too many gaps"
- "Calculate pairwise identity matrix"
- "Show conservation score at each position"
- "Calculate Shannon entropy for each column"

## Requirements

```bash
pip install biopython numpy
```

## Related Skills

- **alignment-files** - Process SAM/BAM/CRAM read alignments
- **phylogenetics** - Build phylogenetic trees from MSAs
- **sequence-io** - Read input sequences for alignment
- **sequence-manipulation** - Work with individual sequences


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->