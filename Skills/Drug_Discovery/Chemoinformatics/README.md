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

# chemoinformatics

## Overview
Computational chemistry for drug discovery covering molecular representations, property prediction, similarity searching, virtual screening, and ADMET analysis.

**Tool type:** python | **Primary tools:** RDKit, DeepChem, AutoDock Vina

## Skills
| Skill | Description |
|-------|-------------|
| molecular-io | Read, write, convert molecular formats (SMILES, SDF, MOL2) |
| molecular-descriptors | Calculate fingerprints and physicochemical properties |
| similarity-searching | Find similar compounds using Tanimoto similarity |
| substructure-search | Filter compounds by SMARTS substructure patterns |
| admet-prediction | Predict ADMET properties and drug-likeness |
| virtual-screening | Dock compounds against protein targets |
| reaction-enumeration | Generate virtual libraries via reaction SMARTS |

## Example Prompts
- "Load my compound library from SDF and standardize structures"
- "Calculate ECFP4 fingerprints for my molecules"
- "Find compounds similar to my lead (Tanimoto > 0.7)"
- "Filter my library for Lipinski rule-of-5 compliant compounds"
- "Predict ADMET properties for my hit compounds"
- "Dock my library against this protein structure"

## Requirements
```bash
pip install rdkit deepchem vina openbabel-wheel
```

## Related Skills

- **structural-biology** - Protein structure handling for docking
- **machine-learning** - ML for activity prediction


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->