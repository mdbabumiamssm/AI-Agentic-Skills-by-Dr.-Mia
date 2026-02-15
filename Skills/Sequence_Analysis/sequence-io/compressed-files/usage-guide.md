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

# Compressed Files - Usage Guide

## Overview

This skill enables AI agents to help you read and write compressed sequence files (.gz, .bz2) using Biopython.

## Prerequisites

```bash
pip install biopython
```

No additional packages needed - gzip and bz2 are built into Python.

## Quick Start

Tell your AI agent what you want to do:

- "Read my gzipped FASTQ file"
- "Compress this FASTA file with gzip"
- "Count sequences in reads.fastq.gz"
- "Convert compressed GenBank to uncompressed FASTA"

## Example Prompts

### Reading
> "Parse the gzipped FASTA and show sequence lengths"

### Writing
> "Save these sequences to a gzip-compressed FASTA"

### Converting
> "Decompress reads.fastq.gz to reads.fastq"

## What the Agent Will Do
1. Detect compression format from file extension
2. Open file with appropriate decompression handler
3. Process sequences from the compressed stream
4. Close file handles properly

## Tips

- Always use text mode ('rt', 'wt') not binary mode with SeqIO
- Gzip is faster, bzip2 compresses better
- Cannot use SeqIO.index() on compressed files
- Process large compressed files as streams to save memory


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->