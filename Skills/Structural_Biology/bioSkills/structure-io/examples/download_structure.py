# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Download structures from RCSB PDB'''

from Bio.PDB import PDBList

pdbl = PDBList()

pdb_id = '4HHB'
file_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='mmCif')
print(f'Downloaded {pdb_id}: {file_path}')

# Download biological assembly
assembly_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='mmCif', assembly_num=1)
print(f'Downloaded assembly: {assembly_path}')

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
