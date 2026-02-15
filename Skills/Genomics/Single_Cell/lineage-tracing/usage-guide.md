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

# Lineage Tracing Analysis - Usage Guide

## Overview

Reconstruct cell lineage trees from CRISPR barcodes, mitochondrial mutations, or other heritable markers.

## Prerequisites

```bash
pip install cassiopeia-lineage cospar
```

## Quick Start

- "Build lineage tree from CRISPR barcodes"
- "Infer clonal relationships between cells"
- "Track cell fate using lineage tracing"

## Example Prompts

### Tree Reconstruction

> "Reconstruct lineage tree from my barcode data"

> "Build phylogeny using Cassiopeia"

### Clonal Dynamics

> "Track how clones expand over time"

> "Infer fate biases for different lineages"

### Mitochondrial Tracing

> "Use mtDNA mutations for lineage tracing"

> "Find clonally related cells without barcodes"

## What the Agent Will Do

1. Process barcode sequences or variant calls
2. Build character matrix from mutations
3. Reconstruct lineage tree
4. Annotate tree with cell metadata
5. Analyze clonal dynamics

## Tips

- **Cassiopeia** - Gold standard for CRISPR barcode trees
- **CoSpar** - Combines lineage with transcriptome dynamics
- **Character matrix** - Mutation states (0=WT, 1,2,3...=mutations, -1=missing)
- **Tree solvers** - Greedy, neighbor-joining, or ILP-based
- **Validation** - Compare tree to known biology (e.g., tissue of origin)


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->