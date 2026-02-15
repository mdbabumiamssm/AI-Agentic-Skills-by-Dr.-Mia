# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Basic sequence writing examples'''
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Create simple records
records = [
    SeqRecord(Seq('ATGCGATCGATCGATCGATCG'), id='seq1', description='First sequence'),
    SeqRecord(Seq('GCTAGCTAGCTAGCTAGCTA'), id='seq2', description='Second sequence'),
    SeqRecord(Seq('TTAATTAATTAATTAATTAA'), id='seq3', description='Third sequence')
]

# Write to FASTA
count = SeqIO.write(records, 'output.fasta', 'fasta')
print(f'Wrote {count} records to output.fasta')

# Get formatted string without writing
for record in records:
    print(record.format('fasta'))

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
