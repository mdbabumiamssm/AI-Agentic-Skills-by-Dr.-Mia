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

# HDR Template Design - Usage Guide

## Overview

Design homology-directed repair donor templates for CRISPR knock-ins including ssODN, dsDNA, and plasmid donors with optimized homology arms.

## Prerequisites

```bash
pip install biopython primer3-py
```

## Quick Start

Tell your AI agent what you want to do:
- "Design an HDR template to add a GFP tag to MYC"
- "Create an ssODN to correct this point mutation"
- "Design a donor plasmid for a conditional knockout"

## Example Prompts

### Protein Tagging

> "Design an ssODN to add a FLAG tag to the C-terminus of my gene"

> "Create an HDR donor to insert GFP at the N-terminus"

### Point Mutations

> "Design a repair template to correct the sickle cell mutation"

> "Create an ssODN to change A to G at this position"

### Large Insertions

> "Design a dsDNA donor to insert a 2kb cassette"

> "Create PCR primers to amplify my HDR donor arms"

### Template Optimization

> "Add a silent PAM mutation to prevent re-cutting"

> "Design asymmetric arms to improve HDR efficiency"

## What the Agent Will Do

1. Identify the cut site and insertion/mutation location
2. Determine optimal template type (ssODN vs dsDNA)
3. Design homology arms of appropriate length
4. Incorporate the desired edit or insertion
5. Optionally add silent PAM mutations
6. Generate ordering sequences or cloning primers

## Tips

- **ssODN length** - Keep total length 100-200nt for efficient synthesis
- **Arm length** - ssODN: 30-60nt each; dsDNA: 200-800bp each
- **Asymmetric arms** - PAM-distal arm longer can improve efficiency
- **PAM disruption** - Mutate PAM silently to prevent re-cutting
- **Strand choice** - For ssODN, test both sense and antisense


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->