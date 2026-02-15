# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import NamedTuple


class _SCVIVA_REGISTRY_KEYS_NT(NamedTuple):
    SAMPLE_KEY: str = "sample"
    Z1_MEAN_KEY: str = "latent_mean"
    NICHE_INDEXES_KEY: str = "niche_indexes"
    NICHE_DISTANCES_KEY: str = "niche_distances"
    NICHE_COMPOSITION_KEY: str = "niche_composition"
    Z1_MEAN_CT_KEY: str = "latent_mean_ct_key"
    CELL_COORDINATES_KEY: str = "spatial"


SCVIVA_REGISTRY_KEYS = _SCVIVA_REGISTRY_KEYS_NT()


class _SCVIVA_MODULE_KEYS(NamedTuple):
    # generative model
    NICHE_MEAN: str = "niche_mean"
    NICHE_VARIANCE: str = "niche_variance"
    P_NICHE_COMPOSITION: str = "niche_composition"
    P_NICHE_EXPRESSION: str = "niche_expression"
    NICHE_ATTENTION: str = "niche_attention"
    # loss
    NLL_NICHE_COMPOSITION_KEY: str = "niche_compo"
    NLL_NICHE_EXPRESSION_KEY: str = "niche_reconst"
    SPATIAL_WEIGHT_KEY: str = "spatial_weight"


SCVIVA_MODULE_KEYS = _SCVIVA_MODULE_KEYS()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
