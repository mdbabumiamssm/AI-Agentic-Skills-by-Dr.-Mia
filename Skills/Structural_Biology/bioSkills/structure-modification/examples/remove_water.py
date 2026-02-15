# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Remove water molecules from a structure'''

from Bio.PDB import PDBParser, PDBIO

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

water_count = 0
for model in structure:
    for chain in model:
        residues_to_remove = [r.id for r in chain if r.id[0] == 'W']
        water_count += len(residues_to_remove)
        for res_id in residues_to_remove:
            chain.detach_child(res_id)

print(f'Removed {water_count} water molecules')

io = PDBIO()
io.set_structure(structure)
io.save('no_water.pdb')
print('Saved to no_water.pdb')

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
