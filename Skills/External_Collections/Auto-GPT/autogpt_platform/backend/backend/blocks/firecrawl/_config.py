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

firecrawl = (
    ProviderBuilder("firecrawl")
    .with_api_key("FIRECRAWL_API_KEY", "Firecrawl API Key")
    .with_base_cost(1, BlockCostType.RUN)
    .build()
)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
