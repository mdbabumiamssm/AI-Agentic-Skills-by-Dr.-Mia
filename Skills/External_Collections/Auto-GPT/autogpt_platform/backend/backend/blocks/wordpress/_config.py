# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from backend.sdk import BlockCostType, ProviderBuilder

from ._oauth import WordPressOAuthHandler, WordPressScope

wordpress = (
    ProviderBuilder("wordpress")
    .with_base_cost(1, BlockCostType.RUN)
    .with_oauth(
        WordPressOAuthHandler,
        scopes=[
            v.value
            for v in [
                WordPressScope.POSTS,
            ]
        ],
        client_id_env_var="WORDPRESS_CLIENT_ID",
        client_secret_env_var="WORDPRESS_CLIENT_SECRET",
    )
    .build()
)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
