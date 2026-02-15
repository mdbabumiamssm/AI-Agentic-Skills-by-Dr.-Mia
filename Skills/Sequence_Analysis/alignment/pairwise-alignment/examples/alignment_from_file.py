# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Align sequences from a FASTA file'''

from Bio import SeqIO
from Bio.Align import PairwiseAligner

fasta_file = 'sequences.fasta'
records = list(SeqIO.parse(fasta_file, 'fasta'))

if len(records) < 2:
    print('Need at least 2 sequences to align')
else:
    seq1, seq2 = records[0].seq, records[1].seq

    aligner = PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-0.5)
    alignments = aligner.align(seq1, seq2)

    print(f'Aligning {records[0].id} vs {records[1].id}')
    print(f'Score: {alignments[0].score}\n')
    print(alignments[0])

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
