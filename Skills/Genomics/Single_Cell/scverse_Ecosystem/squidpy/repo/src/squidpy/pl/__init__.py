# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""The plotting module."""

from __future__ import annotations

from squidpy.pl._graph import (
    centrality_scores,
    co_occurrence,
    interaction_matrix,
    nhood_enrichment,
    ripley,
)

# from squidpy.pl._interactive import Interactive  # type: ignore[attr-defined] # deprecated
from squidpy.pl._ligrec import ligrec
from squidpy.pl._spatial import spatial_scatter, spatial_segment
from squidpy.pl._utils import extract
from squidpy.pl._var_by_distance import var_by_distance

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
