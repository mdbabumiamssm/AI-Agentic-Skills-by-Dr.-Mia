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

# Substructure Search - Usage Guide

## Overview
Search molecular libraries for specific structural patterns using SMARTS. Filter compounds by functional groups, pharmacophores, or scaffold features.

## Prerequisites
```bash
pip install rdkit
```

## Quick Start
Tell your AI agent what you want to do:
- "Find all compounds with a carboxylic acid group"
- "Filter out molecules containing reactive groups"
- "Search for benzimidazole scaffold in my library"
- "Highlight the matched atoms in my molecule"

## Example Prompts

### Pattern Searching
> "Find all molecules containing a primary amine."

> "Search for compounds with both hydroxyl and carbonyl groups."

### Library Filtering
> "Remove compounds with nitro groups or Michael acceptors."

> "Keep only molecules with aromatic rings."

### Visualization
> "Show me the carboxylic acid group highlighted in this molecule."

> "Mark the matched substructure in red."

## What the Agent Will Do
1. Parse SMARTS pattern
2. Search molecules for substructure matches
3. Return atom indices of matches
4. Filter library by inclusion/exclusion criteria
5. Generate visualizations with highlighting

## Tips
- SMARTS uses square brackets for atom properties: [OH], [NH2], [CX3]
- Aromatic atoms are lowercase: c for aromatic carbon, C for aliphatic
- Use recursive SMARTS for complex patterns: $(...) notation
- Common patterns: [OH] hydroxyl, [NH2] primary amine, c1ccccc1 benzene
- HasSubstructMatch is faster than GetSubstructMatches for presence check

## Related Skills
- molecular-io - Load molecules for searching
- similarity-searching - Fingerprint-based searching
- admet-prediction - Filter before ADMET analysis


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->