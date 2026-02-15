# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Transform structure coordinates'''

from Bio.PDB import PDBParser, PDBIO
import numpy as np

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

coords = np.array([a.coord for a in structure.get_atoms()])
center = coords.mean(axis=0)
print(f'Original center: {center}')

for atom in structure.get_atoms():
    atom.coord = atom.coord - center

coords = np.array([a.coord for a in structure.get_atoms()])
new_center = coords.mean(axis=0)
print(f'New center: {new_center}')

io = PDBIO()
io.set_structure(structure)
io.save('centered.pdb')
print('Saved centered structure to centered.pdb')

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
