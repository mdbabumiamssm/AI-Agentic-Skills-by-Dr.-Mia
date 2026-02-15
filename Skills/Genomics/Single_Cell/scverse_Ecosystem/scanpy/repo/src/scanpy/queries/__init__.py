# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Biomart queries."""

from __future__ import annotations

from ._queries import (
    biomart_annotations,
    enrich,  # gprofiler queries
    gene_coordinates,
    mitochondrial_genes,
)

__all__ = [
    "biomart_annotations",
    "enrich",
    "gene_coordinates",
    "mitochondrial_genes",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
