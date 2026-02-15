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
name: bio-crispr-screens-base-editing-analysis
description: Analyzes base editing and prime editing outcomes including editing efficiency, bystander edits, and indel frequencies. Use when quantifying CRISPR base editor results, comparing ABE vs CBE efficiency, or assessing prime editing fidelity.
tool_type: python
primary_tool: CRISPResso2
measurable_outcome: Execute skill workflow successfully with valid output within 15 minutes.
allowed-tools:
  - read_file
  - run_shell_command
---

# Base Editing Analysis

## CRISPResso2 for Base Editing

```bash
# Analyze base editing with expected outcome
CRISPResso --fastq_r1 reads.fq.gz \
    --amplicon_seq ATGCGATCGATCGATCGATCGATCG \
    --guide_seq TCGATCGATCGATCGAT \
    --expected_hdr_amplicon_seq ATGCGATCGATCGTTCGATCGATCG \
    --base_editor_output \
    -o results/
```

## Key Metrics

| Metric | Description |
|--------|-------------|
| Editing efficiency | % reads with target base change |
| Bystander edits | Unintended edits in editing window |
| Indel frequency | Insertions/deletions (should be low) |
| Purity | Target edit without bystanders |

## Base Editor Types

### Cytosine Base Editors (CBE)

```bash
# C->T conversion (or G->A on opposite strand)
CRISPResso --fastq_r1 reads.fq.gz \
    --amplicon_seq $AMPLICON \
    --guide_seq $GUIDE \
    --base_editor_output \
    --conversion_nuc_from C \
    --conversion_nuc_to T
```

### Adenine Base Editors (ABE)

```bash
# A->G conversion (or T->C on opposite strand)
CRISPResso --fastq_r1 reads.fq.gz \
    --amplicon_seq $AMPLICON \
    --guide_seq $GUIDE \
    --base_editor_output \
    --conversion_nuc_from A \
    --conversion_nuc_to G
```

## Prime Editing Analysis

```bash
# Prime editing with pegRNA
CRISPResso --fastq_r1 reads.fq.gz \
    --amplicon_seq $AMPLICON \
    --guide_seq $SPACER \
    --expected_hdr_amplicon_seq $EDITED_AMPLICON \
    --prime_editing_pegRNA_extension_seq $EXTENSION \
    -o prime_edit_results/
```

## Editing Window Analysis

```python
import pandas as pd

# Load CRISPResso quantification
quant = pd.read_csv('CRISPResso_output/Quantification_window_nucleotide_percentage_table.txt',
                    sep='\t')

# Calculate per-position editing
editing_window = quant[(quant['Position'] >= -5) & (quant['Position'] <= 5)]
```

## Quality Thresholds

- Editing efficiency: >30% considered good for most applications
- Indel rate: <5% ideal for base editors
- Bystander rate: depends on application; <10% often acceptable

## Related Skills

- crispr-screens/crispresso-editing - General editing QC
- crispr-screens/library-design - Guide design considerations
- variant-calling/vcf-basics - Downstream variant analysis


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->