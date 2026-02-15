# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Slice and subset alignments'''

from Bio import AlignIO

alignment = AlignIO.read('alignment.aln', 'clustal')
print(f'Original: {len(alignment)} sequences, {alignment.get_alignment_length()} columns')

subset_seqs = alignment[0:5]
print(f'First 5 sequences: {len(subset_seqs)} sequences')

trimmed = alignment[:, 50:150]
print(f'Columns 50-150: {trimmed.get_alignment_length()} columns')

region = alignment[0:5, 50:150]
print(f'Combined slice: {len(region)} sequences, {region.get_alignment_length()} columns')

AlignIO.write(region, 'trimmed_subset.fasta', 'fasta')
print('Wrote trimmed subset to trimmed_subset.fasta')

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
