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
Shared configuration for all Airtable blocks using the SDK pattern.
"""

from backend.sdk import BlockCostType, ProviderBuilder

from ._oauth import AirtableOAuthHandler, AirtableScope
from ._webhook import AirtableWebhookManager

# Configure the Airtable provider with API key authentication
airtable = (
    ProviderBuilder("airtable")
    .with_api_key("AIRTABLE_API_KEY", "Airtable Personal Access Token")
    .with_webhook_manager(AirtableWebhookManager)
    .with_base_cost(1, BlockCostType.RUN)
    .with_oauth(
        AirtableOAuthHandler,
        scopes=[
            v.value
            for v in [
                AirtableScope.DATA_RECORDS_READ,
                AirtableScope.DATA_RECORDS_WRITE,
                AirtableScope.SCHEMA_BASES_READ,
                AirtableScope.SCHEMA_BASES_WRITE,
                AirtableScope.WEBHOOK_MANAGE,
            ]
        ],
        client_id_env_var="AIRTABLE_CLIENT_ID",
        client_secret_env_var="AIRTABLE_CLIENT_SECRET",
    )
    .build()
)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
