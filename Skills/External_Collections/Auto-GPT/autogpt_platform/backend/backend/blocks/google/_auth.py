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

# --8<-- [start:GoogleOAuthIsConfigured]
secrets = Secrets()
GOOGLE_OAUTH_IS_CONFIGURED = bool(
    secrets.google_client_id and secrets.google_client_secret
)
# --8<-- [end:GoogleOAuthIsConfigured]
GoogleCredentials = OAuth2Credentials
GoogleCredentialsInput = CredentialsMetaInput[
    Literal[ProviderName.GOOGLE], Literal["oauth2"]
]


def GoogleCredentialsField(scopes: list[str]) -> GoogleCredentialsInput:
    """
    Creates a Google credentials input on a block.

    Params:
        scopes: The authorization scopes needed for the block to work.
    """
    return CredentialsField(
        required_scopes=set(scopes),
        description="The Google integration requires OAuth2 authentication.",
    )


TEST_CREDENTIALS = OAuth2Credentials(
    id="01234567-89ab-cdef-0123-456789abcdef",
    provider="google",
    access_token=SecretStr("mock-google-access-token"),
    refresh_token=SecretStr("mock-google-refresh-token"),
    access_token_expires_at=1234567890,
    scopes=[
        "https://www.googleapis.com/auth/gmail.readonly",
        "https://www.googleapis.com/auth/gmail.send",
    ],
    title="Mock Google OAuth2 Credentials",
    username="mock-google-username",
    refresh_token_expires_at=1234567890,
)

TEST_CREDENTIALS_INPUT = {
    "provider": TEST_CREDENTIALS.provider,
    "id": TEST_CREDENTIALS.id,
    "type": TEST_CREDENTIALS.type,
    "title": TEST_CREDENTIALS.title,
}

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
