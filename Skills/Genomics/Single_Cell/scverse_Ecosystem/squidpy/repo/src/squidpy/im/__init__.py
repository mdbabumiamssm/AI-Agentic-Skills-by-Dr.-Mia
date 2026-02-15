# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""The image module."""

from __future__ import annotations

from squidpy.im._container import ImageContainer
from squidpy.im._feature import calculate_image_features
from squidpy.im._process import process
from squidpy.im._segment import (
    SegmentationCustom,
    SegmentationModel,
    SegmentationWatershed,
    segment,
)

__all__ = [
    "ImageContainer",
    "calculate_image_features",
    "process",
    "SegmentationCustom",
    "SegmentationModel",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
