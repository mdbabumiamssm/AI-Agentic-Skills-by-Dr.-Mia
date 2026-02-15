# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Balance a Hi-C matrix'''

import cooler
import cooltools
import numpy as np

clr = cooler.Cooler('matrix.cool')

print('Before balancing:')
print(f'  Has weights: {"weight" in clr.bins().columns}')

print('\nBalancing (ICE)...')
cooler.balance_cooler('matrix.cool', store=True, cis_only=True)

clr = cooler.Cooler('matrix.cool')
weights = clr.bins()['weight'][:]
print(f'\nAfter balancing:')
print(f'  Weight range: {weights.min():.4f} - {weights.max():.4f}')
print(f'  NaN weights: {np.isnan(weights).sum()}')

print('\nComputing expected values...')
expected = cooltools.expected_cis(clr, ignore_diags=2)
print(expected.head(10))

balanced = clr.matrix(balance=True).fetch('chr1')
print(f'\nBalanced chr1 matrix:')
print(f'  Shape: {balanced.shape}')
print(f'  Range: {np.nanmin(balanced):.6f} - {np.nanmax(balanced):.6f}')

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
