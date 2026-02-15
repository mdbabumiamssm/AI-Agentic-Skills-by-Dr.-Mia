# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from typing import Literal

from pydantic import SecretStr

from backend.data.model import CredentialsField, CredentialsMetaInput, OAuth2Credentials
from backend.integrations.providers import ProviderName
from backend.util.settings import Secrets

secrets = Secrets()
NOTION_OAUTH_IS_CONFIGURED = bool(
    secrets.notion_client_id and secrets.notion_client_secret
)

NotionCredentials = OAuth2Credentials
NotionCredentialsInput = CredentialsMetaInput[
    Literal[ProviderName.NOTION], Literal["oauth2"]
]


def NotionCredentialsField() -> NotionCredentialsInput:
    """Creates a Notion OAuth2 credentials field."""
    return CredentialsField(
        description="Connect your Notion account. Ensure the pages/databases are shared with the integration."
    )


# Test credentials for Notion OAuth2
TEST_CREDENTIALS = OAuth2Credentials(
    id="01234567-89ab-cdef-0123-456789abcdef",
    provider="notion",
    access_token=SecretStr("test_access_token"),
    title="Mock Notion OAuth",
    scopes=["read_content", "insert_content", "update_content"],
    username="testuser",
)

TEST_CREDENTIALS_INPUT = {
    "provider": TEST_CREDENTIALS.provider,
    "id": TEST_CREDENTIALS.id,
    "type": TEST_CREDENTIALS.type,
    "title": TEST_CREDENTIALS.title,
}

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
