# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from __future__ import annotations

from ._core.sparse_dataset import sparse_dataset
from ._io.h5ad import read_h5ad, write_h5ad
from ._io.read import (
    read_csv,
    read_excel,
    read_hdf,
    read_loom,  # noqa: F401
    read_mtx,
    read_text,
    read_umi_tools,
)
from ._io.specs import read_elem, write_elem
from ._io.write import write_csvs, write_loom  # noqa: F401
from ._io.zarr import read_zarr, write_zarr

__all__ = [
    "read_csv",
    "read_elem",
    "read_excel",
    "read_h5ad",
    "read_hdf",
    "read_mtx",
    "read_text",
    "read_umi_tools",
    "read_zarr",
    "sparse_dataset",
    "write_csvs",
    "write_elem",
    "write_h5ad",
    "write_zarr",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
