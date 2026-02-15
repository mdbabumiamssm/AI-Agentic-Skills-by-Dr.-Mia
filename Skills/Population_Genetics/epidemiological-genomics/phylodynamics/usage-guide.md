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

# Phylodynamics - Usage Guide

## Overview

Build time-scaled phylogenies and infer evolutionary dynamics from pathogen sequences using TreeTime for outbreak analysis.

## Prerequisites

```bash
pip install phylo-treetime biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Build a time-scaled tree for this outbreak"
- "Estimate when this outbreak started"
- "Calculate the molecular clock rate for these sequences"

## Example Prompts

### Time Trees

> "Create a dated phylogeny from my sequences and collection dates"

> "When did the most recent common ancestor of these isolates exist?"

### Clock Analysis

> "What is the substitution rate for this virus?"

> "Is there clock-like evolution in my sequences?"

### Population Dynamics

> "Show how the outbreak population size changed over time"

> "Generate a skyline plot from my phylogeny"

### Geographic Spread

> "Reconstruct the geographic spread of this pathogen"

> "Where did this outbreak likely originate?"

## What the Agent Will Do

1. Load phylogenetic tree and sample dates
2. Optimize root position for molecular clock
3. Estimate clock rate and confidence intervals
4. Date internal nodes and root
5. Optionally reconstruct ancestral states
6. Generate time-scaled tree and reports

## Tips

- **Date format** - Use decimal years (2020.5) or ISO dates (2020-06-15)
- **Root optimization** - Let TreeTime find optimal root unless biologically known
- **Clock rate** - RNA viruses: 10^-3 to 10^-4; bacteria: 10^-6 to 10^-7
- **Temporal signal** - Need sufficient date range for clock estimation
- **Skyline prior** - Use for outbreak population dynamics analysis


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->