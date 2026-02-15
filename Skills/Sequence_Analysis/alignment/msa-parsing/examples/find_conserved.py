# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Find conserved positions in an alignment'''

from Bio import AlignIO
from collections import Counter

alignment = AlignIO.read('alignment.fasta', 'fasta')

def find_conserved_positions(alignment, threshold=0.8):
    conserved = []
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx]
        counts = Counter(column)
        if '-' in counts:
            del counts['-']
        if not counts:
            continue
        most_common_char, most_common_count = counts.most_common(1)[0]
        conservation = most_common_count / len(alignment)
        if conservation >= threshold:
            conserved.append((col_idx, most_common_char, conservation))
    return conserved

print('Fully conserved positions (100%):')
for pos, char, cons in find_conserved_positions(alignment, threshold=1.0):
    print(f'  Position {pos}: {char}')

print('\nHighly conserved positions (80%+):')
for pos, char, cons in find_conserved_positions(alignment, threshold=0.8):
    print(f'  Position {pos}: {char} ({cons*100:.0f}%)')

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
