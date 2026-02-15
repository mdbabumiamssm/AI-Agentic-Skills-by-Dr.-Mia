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

from squidpy.datasets._datasets import (
    # Type aliases for dataset names
    AnnDataDatasets,
    ImageDatasets,
    SpatialDataDatasets,
    VisiumDatasets,
    # AnnData datasets
    four_i,
    imc,
    merfish,
    mibitof,
    sc_mouse_cortex,
    seqfish,
    slideseqv2,
    # 10x Genomics Visium
    visium,
    visium_fluo_adata,
    visium_fluo_adata_crop,
    # Image datasets
    visium_fluo_image_crop,
    visium_hne_adata,
    visium_hne_adata_crop,
    visium_hne_image,
    visium_hne_image_crop,
    visium_hne_sdata,
)

__all__ = [
    # Type aliases
    "VisiumDatasets",
    "AnnDataDatasets",
    "ImageDatasets",
    "SpatialDataDatasets",
    # Datasets by format:
    # AnnData
    "four_i",
    "imc",
    "seqfish",
    "visium_hne_adata",
    "visium_hne_adata_crop",
    "visium_fluo_adata",
    "visium_fluo_adata_crop",
    "sc_mouse_cortex",
    "mibitof",
    "merfish",
    "slideseqv2",
    # AnnData with Image
    "visium",
    # Image
    "visium_fluo_image_crop",
    "visium_hne_image_crop",
    "visium_hne_image",
    # SpatialData
    "visium_hne_sdata",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
