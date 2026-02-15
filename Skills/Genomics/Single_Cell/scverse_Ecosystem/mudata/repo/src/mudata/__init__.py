# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Multimodal datasets"""

try:  # See https://github.com/maresb/hatch-vcs-footgun-example
    from setuptools_scm import get_version

    __version__ = get_version(root="../..", relative_to=__file__)
except (ImportError, LookupError):
    try:
        from ._version import __version__
    except ModuleNotFoundError:
        raise RuntimeError("mudata is not correctly installed. Please install it, e.g. with pip.")

from anndata import AnnData

from ._core import utils
from ._core.config import set_options
from ._core.io import (
    read,
    read_anndata,
    read_h5ad,
    read_h5mu,
    read_zarr,
    write,
    write_anndata,
    write_h5ad,
    write_h5mu,
    write_zarr,
)
from ._core.merge import concat
from ._core.mudata import MuData
from ._core.to_ import to_anndata, to_mudata

__anndataversion__ = "0.1.0"
__mudataversion__ = "0.1.0"

__all__ = [
    "__version__",
    "MuData",
    "AnnData",
    "utils",
    "set_options",
    "to_anndata",
    "to_mudata",
    "concat",
    "read",
    "read_h5ad",
    "read_anndata",
    "read_h5mu",
    "read_zarr",
    "write",
    "write_h5ad",
    "write_anndata",
    "write_h5mu",
    "write_zarr",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
