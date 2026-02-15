# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Experimental preprocessing functions."""

from __future__ import annotations

from scanpy.experimental.pp._highly_variable_genes import highly_variable_genes
from scanpy.experimental.pp._normalization import (
    normalize_pearson_residuals,
    normalize_pearson_residuals_pca,
)
from scanpy.experimental.pp._recipes import recipe_pearson_residuals

__all__ = [
    "highly_variable_genes",
    "normalize_pearson_residuals",
    "normalize_pearson_residuals_pca",
    "recipe_pearson_residuals",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
