# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Get data from AnnData."""

from __future__ import annotations

from ._aggregated import aggregate
from .get import (
    _check_mask,
    _get_obs_rep,
    _ObsRep,
    _set_obs_rep,
    obs_df,
    rank_genes_groups_df,
    var_df,
)

__all__ = [
    "_ObsRep",
    "_check_mask",
    "_get_obs_rep",
    "_set_obs_rep",
    "aggregate",
    "obs_df",
    "rank_genes_groups_df",
    "var_df",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
