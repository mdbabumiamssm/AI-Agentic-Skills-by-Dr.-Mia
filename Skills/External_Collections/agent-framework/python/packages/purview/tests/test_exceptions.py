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

"""Tests for Purview exceptions."""

from agent_framework_purview import (
    PurviewAuthenticationError,
    PurviewPaymentRequiredError,
    PurviewRateLimitError,
    PurviewRequestError,
    PurviewServiceError,
)


class TestPurviewExceptions:
    """Test custom Purview exception classes."""

    def test_purview_service_error(self) -> None:
        """Test PurviewServiceError base exception."""
        error = PurviewServiceError("Service error occurred")
        assert str(error) == "Service error occurred"
        assert isinstance(error, Exception)

    def test_purview_authentication_error(self) -> None:
        """Test PurviewAuthenticationError exception."""
        error = PurviewAuthenticationError("Authentication failed")
        assert str(error) == "Authentication failed"
        assert isinstance(error, PurviewServiceError)

    def test_purview_payment_required_error(self) -> None:
        """Test PurviewPaymentRequiredError exception."""
        error = PurviewPaymentRequiredError("Payment required")
        assert str(error) == "Payment required"
        assert isinstance(error, PurviewServiceError)

    def test_purview_rate_limit_error(self) -> None:
        """Test PurviewRateLimitError exception."""
        error = PurviewRateLimitError("Rate limit exceeded")
        assert str(error) == "Rate limit exceeded"
        assert isinstance(error, PurviewServiceError)

    def test_purview_request_error(self) -> None:
        """Test PurviewRequestError exception."""
        error = PurviewRequestError("Request failed")
        assert str(error) == "Request failed"
        assert isinstance(error, PurviewServiceError)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
