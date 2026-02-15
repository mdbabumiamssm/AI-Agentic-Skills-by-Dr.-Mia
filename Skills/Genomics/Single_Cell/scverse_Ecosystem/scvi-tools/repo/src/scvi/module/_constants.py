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


class _MODULE_KEYS(NamedTuple):
    X_KEY: str = "x"
    # inference
    Z_KEY: str = "z"
    QZ_KEY: str = "qz"
    QZM_KEY: str = "qzm"
    QZV_KEY: str = "qzv"
    LIBRARY_KEY: str = "library"
    QL_KEY: str = "ql"
    BATCH_INDEX_KEY: str = "batch_index"
    Y_KEY: str = "y"
    CONT_COVS_KEY: str = "cont_covs"
    CAT_COVS_KEY: str = "cat_covs"
    SIZE_FACTOR_KEY: str = "size_factor"
    # generative
    PX_KEY: str = "px"
    PL_KEY: str = "pl"
    PZ_KEY: str = "pz"
    # loss
    KL_L_KEY: str = "kl_divergence_l"
    KL_Z_KEY: str = "kl_divergence_z"


MODULE_KEYS = _MODULE_KEYS()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
