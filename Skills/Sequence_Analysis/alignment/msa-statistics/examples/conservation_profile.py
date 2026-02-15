# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Calculate conservation score at each position'''

from Bio import AlignIO
from collections import Counter

def column_conservation(alignment, col_idx, ignore_gaps=True):
    column = alignment[:, col_idx]
    if ignore_gaps:
        column = column.replace('-', '')
    if not column:
        return 0.0
    counts = Counter(column)
    most_common_count = counts.most_common(1)[0][1]
    return most_common_count / len(column)

alignment = AlignIO.read('alignment.fasta', 'fasta')
print(f'Alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} columns\n')

conservation_scores = []
for i in range(alignment.get_alignment_length()):
    score = column_conservation(alignment, i)
    conservation_scores.append(score)

print('Conservation profile:')
for i, score in enumerate(conservation_scores):
    bar = '#' * int(score * 20)
    print(f'  {i:3d}: {score*100:5.1f}% {bar}')

avg = sum(conservation_scores) / len(conservation_scores)
print(f'\nAverage conservation: {avg*100:.1f}%')

fully_conserved = sum(1 for s in conservation_scores if s == 1.0)
print(f'Fully conserved columns: {fully_conserved} ({fully_conserved/len(conservation_scores)*100:.1f}%)')

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
