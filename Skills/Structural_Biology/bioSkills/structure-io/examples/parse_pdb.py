# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Parse a PDB file and inspect structure contents'''

from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
structure = parser.get_structure('1abc', '1abc.pdb')

print(f'Structure: {structure.id}')
print(f'Models: {len(list(structure.get_models()))}')
print(f'Chains: {len(list(structure.get_chains()))}')
print(f'Residues: {len(list(structure.get_residues()))}')
print(f'Atoms: {len(list(structure.get_atoms()))}')

for model in structure:
    for chain in model:
        residues = list(chain.get_residues())
        print(f'  Chain {chain.id}: {len(residues)} residues')

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
