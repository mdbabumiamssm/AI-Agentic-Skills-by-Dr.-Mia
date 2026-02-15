# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import pytest
import pytest_socket
import requests


def test_socket_disabled() -> None:
    """This test should fail."""
    with pytest.raises(pytest_socket.SocketBlockedError):
        # Ignore S113 since we don't need a timeout here as the request
        # should fail immediately
        requests.get("https://www.example.com")  # noqa: S113

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
