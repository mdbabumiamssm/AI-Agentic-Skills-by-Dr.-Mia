# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Convert alignment between different formats'''

from Bio import AlignIO

input_file = 'alignment.aln'
input_format = 'clustal'

conversions = [
    ('output.fasta', 'fasta'),
    ('output.phy', 'phylip-relaxed'),
    ('output.nex', 'nexus'),
]

alignment = AlignIO.read(input_file, input_format)
print(f'Read alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} columns')

for output_file, output_format in conversions:
    AlignIO.write(alignment, output_file, output_format)
    print(f'Wrote: {output_file} ({output_format})')

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
