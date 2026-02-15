# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Read and inspect a multiple sequence alignment'''

from Bio import AlignIO

alignment = AlignIO.read('sample_alignment.aln', 'clustal')

print(f'Alignment length: {alignment.get_alignment_length()} columns')
print(f'Number of sequences: {len(alignment)}')
print(f'\nSequence IDs:')
for record in alignment:
    print(f'  {record.id}: {len(record.seq)} bp')

print(f'\nFirst 60 columns of first sequence:')
print(alignment[0].seq[:60])

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
