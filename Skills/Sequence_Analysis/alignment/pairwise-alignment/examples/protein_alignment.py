# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Protein alignment using BLOSUM62 substitution matrix'''

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

protein1 = Seq('MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIAAALEHHHHHH')
protein2 = Seq('MKTAYIAKQRQISFVKSHFSRQLEERLDLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIA')

blosum62 = substitution_matrices.load('BLOSUM62')
aligner = PairwiseAligner(mode='global', substitution_matrix=blosum62, open_gap_score=-11, extend_gap_score=-1)

alignments = aligner.align(protein1, protein2)
print(f'Score: {alignments[0].score}\n')
print(alignments[0])

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
