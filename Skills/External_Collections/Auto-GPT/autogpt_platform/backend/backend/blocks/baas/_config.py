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
Shared configuration for all Meeting BaaS blocks using the SDK pattern.
"""

from backend.sdk import BlockCostType, ProviderBuilder

# Configure the Meeting BaaS provider with API key authentication
baas = (
    ProviderBuilder("baas")
    .with_api_key("MEETING_BAAS_API_KEY", "Meeting BaaS API Key")
    .with_base_cost(5, BlockCostType.RUN)  # Higher cost for meeting recording service
    .build()
)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
