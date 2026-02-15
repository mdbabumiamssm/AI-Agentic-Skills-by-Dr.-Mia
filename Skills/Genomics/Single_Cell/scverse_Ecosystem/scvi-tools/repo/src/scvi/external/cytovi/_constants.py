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


class _CYTOVI_REGISTRY_KEYS_NT(NamedTuple):
    X_KEY: str = "X"
    BATCH_KEY: str = "batch"
    LABELS_KEY: str = "labels"
    SAMPLE_KEY: str = "sample_id"
    CAT_COVS_KEY: str = "extra_categorical_covs"
    CONT_COVS_KEY: str = "extra_continuous_covs"
    INDICES_KEY: str = "ind_x"
    SIZE_FACTOR_KEY: str = "size_factor"
    MINIFY_TYPE_KEY: str = "minify_type"
    LATENT_QZM_KEY: str = "latent_qzm"
    LATENT_QZV_KEY: str = "latent_qzv"
    PROTEIN_NAN_MASK: str = "nan_layer"


CYTOVI_REGISTRY_KEYS = _CYTOVI_REGISTRY_KEYS_NT()

CYTOVI_DEFAULT_REP = "X_CytoVI"
CYTOVI_SCATTER_FEATS = ("FSC", "Fsc", "fsc", "SSC", "Ssc", "ssc")

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
