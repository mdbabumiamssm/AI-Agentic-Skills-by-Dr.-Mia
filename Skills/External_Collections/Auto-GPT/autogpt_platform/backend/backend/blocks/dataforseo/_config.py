# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""
Configuration for all DataForSEO blocks using the new SDK pattern.
"""

from backend.sdk import BlockCostType, ProviderBuilder

# Build the DataForSEO provider with username/password authentication
dataforseo = (
    ProviderBuilder("dataforseo")
    .with_user_password(
        username_env_var="DATAFORSEO_USERNAME",
        password_env_var="DATAFORSEO_PASSWORD",
        title="DataForSEO Credentials",
    )
    .with_base_cost(1, BlockCostType.RUN)
    .build()
)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
