# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

'''Plot Hi-C contact matrix'''

import cooler
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

clr = cooler.Cooler('matrix.mcool::resolutions/10000')
print(f'Loaded at {clr.binsize}bp resolution')

region = 'chr1:50000000-60000000'
matrix = clr.matrix(balance=True).fetch(region)
print(f'Matrix shape: {matrix.shape}')

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

im1 = axes[0].imshow(matrix, cmap='Reds', norm=LogNorm(vmin=0.001, vmax=0.1))
axes[0].set_title(f'{region} (balanced)')
plt.colorbar(im1, ax=axes[0], label='Contacts')

log_matrix = np.log2(matrix + 1e-10)
im2 = axes[1].imshow(log_matrix, cmap='Reds', vmin=-10, vmax=-3)
axes[1].set_title(f'{region} (log2)')
plt.colorbar(im2, ax=axes[1], label='log2(contacts)')

plt.tight_layout()
plt.savefig('hic_matrix.png', dpi=150)
print('Saved to hic_matrix.png')

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
