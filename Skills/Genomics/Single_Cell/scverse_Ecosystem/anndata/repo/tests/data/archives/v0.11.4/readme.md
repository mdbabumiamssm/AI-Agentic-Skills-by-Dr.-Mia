<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

These files were written with

```bash
uvx '--with=anndata==0.11.4' '--with=zarr<3' python -c '
import zarr
from anndata import AnnData
adata = AnnData(shape=(10, 20))
adata.write_zarr(zarr.ZipStore("tests/data/archives/v0.11.4/adata.zarr.zip"))
adata.write_h5ad("tests/data/archives/v0.11.4/adata.h5ad")'
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->