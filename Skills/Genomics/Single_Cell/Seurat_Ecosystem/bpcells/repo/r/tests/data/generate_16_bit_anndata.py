# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Writes a uint16_t anndata file, with versioning based on currently installed anndata version.

import anndata as ad
import numpy as np
import scipy.sparse

if __name__ == "__main__":
    x = np.array([[1., 2., 0., 3., 0.],
                  [0., 0., 2., 2., 1.],
                  [0., 0., 2., 0., 2.]],
                 dtype = np.uint16)

    adata = ad.AnnData(
        X = scipy.sparse.csr_matrix(x),
        layers = {"transpose": scipy.sparse.csc_matrix(x), "dense": x},
        obsm = {"obs_mat": scipy.sparse.csc_matrix(x[:,:2])},
        varm = {"var_mat": x[:2,:].T},
        )

    adata.write_h5ad(f"mini_mat_uint_16.anndata-v{ad.__version__}.h5ad")
__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
