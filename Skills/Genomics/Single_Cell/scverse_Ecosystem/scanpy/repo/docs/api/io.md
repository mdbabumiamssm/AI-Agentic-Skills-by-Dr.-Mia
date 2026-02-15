<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal AI Agentic Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

(reading)=
(reading-and-writing)=

## Reading and Writing

```{eval-rst}
.. currentmodule:: scanpy
```

Write {class}`~anndata.AnnData` objects using its {doc}`writing <anndata:api>` methods

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   write
```

```{note}
For reading annotation use {ref}`pandas.read_… <pandas:io>`
and add it to your {class}`~anndata.AnnData` object. The following read functions are
intended for the numeric data in the data matrix `X`.
```

Read common file formats using

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   read
```

Read 10x formatted hdf5 files and directories containing `.mtx` files using

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   read_10x_h5
   read_10x_mtx
   read_visium
```

Read other formats using functions borrowed from {mod}`anndata`

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   read_h5ad
   read_csv
   read_excel
   read_hdf
   read_loom
   read_mtx
   read_text
   read_umi_tools

```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->