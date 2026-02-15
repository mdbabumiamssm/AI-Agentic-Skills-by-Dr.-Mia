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

# OpenMM Agent (Molecular Dynamics)

**ID:** `biomedical.drug_discovery.molecular_dynamics`
**Version:** 1.0.0
**Status:** Alpha
**Category:** Drug Discovery / Simulation

---

## Overview

The **OpenMM Agent** automates the setup and execution of Molecular Dynamics (MD) simulations. While AlphaFold provides static structures, this agent reveals how proteins *move* and change conformation over time, which is critical for understanding drug binding and allosteric effects.

## Key Capabilities

- **System Setup:** Cleans PDBs, adds hydrogens, solvates in water box, adds ions (neutralization).
- **Simulation Control:** Writes and executes Python scripts for Energy Minimization, Equilibration (NVT/NPT), and Production runs.
- **Analysis:**
    - **RMSD:** Structural stability over time.
    - **RMSF:** Flexibility per residue (identifying binding pockets).

## Integration

Can be chained with **AgentD** or **Docking Agents** to verify ligand stability in the binding pocket.

## References
- *OpenMM Python API*
- *GROMACS*


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->