# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Extract protein sequences from structure'''

from Bio.PDB import PDBParser, PPBuilder

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

ppb = PPBuilder()
for i, pp in enumerate(ppb.build_peptides(structure)):
    seq = pp.get_sequence()
    print(f'Polypeptide {i + 1}: length={len(seq)}')
    print(f'  {seq[:60]}...' if len(seq) > 60 else f'  {seq}')

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
