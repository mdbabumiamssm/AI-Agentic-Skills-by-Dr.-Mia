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
Shared configuration for SDK test providers using the SDK pattern.
"""

from backend.sdk import BlockCostType, ProviderBuilder

# Configure test providers
test_api = (
    ProviderBuilder("test_api")
    .with_api_key("TEST_API_KEY", "Test API Key")
    .with_base_cost(5, BlockCostType.RUN)
    .build()
)

test_service = (
    ProviderBuilder("test_service")
    .with_api_key("TEST_SERVICE_API_KEY", "Test Service API Key")
    .with_base_cost(10, BlockCostType.RUN)
    .build()
)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
