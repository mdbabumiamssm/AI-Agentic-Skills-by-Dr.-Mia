# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

#!/usr/bin/env python3
'''Read and write compressed sequence files'''

import gzip
from Bio import SeqIO
from pathlib import Path

def read_gzipped(filepath, format_type='fasta'):
    '''Read gzipped FASTA/FASTQ'''
    with gzip.open(filepath, 'rt') as handle:
        for record in SeqIO.parse(handle, format_type):
            yield record

def write_gzipped(records, filepath, format_type='fasta'):
    '''Write sequences to gzipped file'''
    with gzip.open(filepath, 'wt') as handle:
        count = SeqIO.write(records, handle, format_type)
    return count

def auto_detect_format(filepath):
    '''Detect format from file extension'''
    path = Path(filepath)
    suffixes = path.suffixes

    is_gzipped = '.gz' in suffixes
    if '.fq' in suffixes or '.fastq' in suffixes:
        return 'fastq', is_gzipped
    else:
        return 'fasta', is_gzipped

def smart_open(filepath):
    '''Open file, handling gzip automatically'''
    if str(filepath).endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')

def count_sequences(filepath):
    '''Count sequences in file'''
    format_type, is_gzipped = auto_detect_format(filepath)

    count = 0
    with smart_open(filepath) as handle:
        for _ in SeqIO.parse(handle, format_type):
            count += 1
    return count

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        filepath = sys.argv[1]
        count = count_sequences(filepath)
        print(f'{filepath}: {count} sequences')

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
