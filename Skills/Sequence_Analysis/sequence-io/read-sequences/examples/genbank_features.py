# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Reading GenBank files and extracting features'''
from Bio import SeqIO

print('=== Parsing GenBank File ===')
for record in SeqIO.parse('sample.gb', 'genbank'):
    print(f'ID: {record.id}')
    print(f'Description: {record.description}')
    print(f'Sequence length: {len(record.seq)} bp')
    print(f'Number of features: {len(record.features)}')

    print('\n=== Features ===')
    for feature in record.features:
        print(f'  Type: {feature.type}')
        print(f'  Location: {feature.location}')
        if feature.type == 'CDS':
            product = feature.qualifiers.get('product', ['Unknown'])[0]
            print(f'  Product: {product}')
        print()

    print('=== Annotations ===')
    for key, value in record.annotations.items():
        print(f'  {key}: {value}')

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
