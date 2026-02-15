# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Generate consensus sequence from alignment'''

from Bio import AlignIO
from collections import Counter

def consensus_sequence(alignment, threshold=0.5, gap_char='-', ambiguous='N'):
    consensus = []
    for col_idx in range(alignment.get_alignment_length()):
        column = alignment[:, col_idx]
        counts = Counter(column)
        most_common_char, most_common_count = counts.most_common(1)[0]
        if most_common_char == gap_char:
            counts.pop(gap_char, None)
            if counts:
                most_common_char, most_common_count = counts.most_common(1)[0]
            else:
                most_common_char = gap_char

        if most_common_count / len(alignment) >= threshold:
            consensus.append(most_common_char)
        else:
            consensus.append(ambiguous)
    return ''.join(consensus)

alignment = AlignIO.read('alignment.fasta', 'fasta')
print(f'Alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} columns\n')

consensus_50 = consensus_sequence(alignment, threshold=0.5, ambiguous='N')
print(f'Consensus (50% threshold):\n{consensus_50}\n')

consensus_70 = consensus_sequence(alignment, threshold=0.7, ambiguous='N')
print(f'Consensus (70% threshold):\n{consensus_70}\n')

consensus_100 = consensus_sequence(alignment, threshold=1.0, ambiguous='N')
print(f'Consensus (100% threshold):\n{consensus_100}')

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
