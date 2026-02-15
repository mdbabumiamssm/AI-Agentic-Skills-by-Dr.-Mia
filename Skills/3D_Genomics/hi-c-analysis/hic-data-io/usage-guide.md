<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Hi-C Data I/O - Usage Guide

## Overview

This skill covers loading, converting, and manipulating Hi-C contact matrices using the cooler format, the standard format for Hi-C data.

## Prerequisites

```bash
pip install cooler bioframe numpy pandas
# For .hic conversion:
pip install hic2cool
```

## Quick Start

Tell your AI agent what you want to do:
- "Load my Hi-C cooler file"
- "Convert .hic to cooler format"

## Example Prompts

### Loading Data
> "Load this cooler file"

> "Open the 10kb resolution from my mcool"

### Extracting Matrices
> "Get the contact matrix for chr1"

> "Extract contacts between these two regions"

### Conversion
> "Convert my .hic file to cooler"

> "Create a multi-resolution mcool"

## What the Agent Will Do

1. Open the cooler file
2. Provide basic statistics (bins, resolution, chromosomes)
3. Extract requested matrix or region
4. Return data as numpy array or DataFrame

## Tips

- **mcool vs cool** - mcool files contain multiple resolutions
- **Balance** - Use `balance=True` for normalized data
- **Sparse** - Use sparse matrices for large chromosomes
- **Resolution** - Choose resolution based on analysis needs


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->