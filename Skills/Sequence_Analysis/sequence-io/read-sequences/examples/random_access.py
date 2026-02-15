# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Random access to sequences using to_dict() and index()'''
from Bio import SeqIO

# Method 1: to_dict() - loads all into memory
print('=== Using to_dict() ===')
records_dict = SeqIO.to_dict(SeqIO.parse('sample.fasta', 'fasta'))
print(f'Available IDs: {list(records_dict.keys())}')
print(f'Accessing seq2: {records_dict["seq2"].seq[:30]}...')

# Method 2: index() - disk-based, memory efficient for large files
print('\n=== Using index() ===')
records_index = SeqIO.index('sample.fasta', 'fasta')
print(f'Available IDs: {list(records_index.keys())}')
print(f'Accessing seq2: {records_index["seq2"].seq[:30]}...')
records_index.close()

# Practical example: extract specific sequences by ID
print('\n=== Extracting Specific Sequences ===')
wanted_ids = {'seq1', 'seq3'}
records_dict = SeqIO.to_dict(SeqIO.parse('sample.fasta', 'fasta'))
for seq_id in wanted_ids:
    if seq_id in records_dict:
        print(f'{seq_id}: {len(records_dict[seq_id].seq)} bp')

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
