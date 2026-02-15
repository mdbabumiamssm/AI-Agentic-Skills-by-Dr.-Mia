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

# restriction-analysis

## Overview

Restriction enzyme analysis using Biopython Bio.Restriction. Find cut sites, create restriction maps, select enzymes for cloning, and predict fragment sizes. Includes data for 800+ enzymes from REBASE.

**Tool type:** python | **Primary tools:** Bio.Restriction

## Skills

| Skill | Description |
|-------|-------------|
| restriction-sites | Find where enzymes cut a sequence |
| restriction-mapping | Create restriction maps, visualize cut positions |
| enzyme-selection | Choose enzymes by criteria (cutters, overhangs, compatibility) |
| fragment-analysis | Predict fragment sizes, simulate gel electrophoresis |

## Example Prompts

- "Find all EcoRI sites in this sequence"
- "Where does BamHI cut in my plasmid?"
- "Show all restriction sites for common cloning enzymes"
- "Create a restriction map of this sequence"
- "Map EcoRI, BamHI, and HindIII sites"
- "Show distances between cut sites"
- "Find enzymes that cut this sequence exactly once"
- "Which enzymes don't cut my insert?"
- "Find enzymes with compatible sticky ends"
- "List all 6-cutter enzymes that cut my sequence"
- "Is my insert compatible with Golden Gate cloning?"
- "Find enzymes not affected by Dam methylation"
- "What fragments will EcoRI produce?"
- "Predict the gel pattern for this digest"
- "Calculate fragment sizes for a double digest"

## Requirements

```bash
pip install biopython
```

## Related Skills

- **sequence-io** - Read sequences for restriction analysis
- **sequence-manipulation** - Work with restriction fragments
- **primer-design** - Design primers around restriction sites


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->