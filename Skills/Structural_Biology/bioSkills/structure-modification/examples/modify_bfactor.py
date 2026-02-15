# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Modify B-factors to encode custom data'''

from Bio.PDB import PDBParser, PDBIO

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

conservation_scores = {i: (i % 10) * 10 for i in range(1, 200)}

for residue in structure.get_residues():
    if residue.id[0] != ' ':
        continue
    resnum = residue.id[1]
    score = conservation_scores.get(resnum, 50.0)
    for atom in residue:
        atom.bfactor = score

io = PDBIO()
io.set_structure(structure)
io.save('colored_by_conservation.pdb')
print('B-factors set to conservation scores')
print('Open in PyMOL and color by B-factor: spectrum b')

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
