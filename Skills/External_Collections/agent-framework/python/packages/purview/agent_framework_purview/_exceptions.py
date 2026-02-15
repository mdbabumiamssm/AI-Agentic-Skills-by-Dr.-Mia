# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Copyright (c) Microsoft. All rights reserved.
"""Purview specific exceptions (minimal error shaping)."""

from agent_framework.exceptions import ServiceResponseException

__all__ = [
    "PurviewAuthenticationError",
    "PurviewPaymentRequiredError",
    "PurviewRateLimitError",
    "PurviewRequestError",
    "PurviewServiceError",
]


class PurviewServiceError(ServiceResponseException):
    """Base exception for Purview errors."""


class PurviewAuthenticationError(PurviewServiceError):
    """Authentication / authorization failure (401/403)."""


class PurviewPaymentRequiredError(PurviewServiceError):
    """Payment required (402)."""


class PurviewRateLimitError(PurviewServiceError):
    """Rate limiting or throttling (429)."""


class PurviewRequestError(PurviewServiceError):
    """Other non-success HTTP errors."""

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
