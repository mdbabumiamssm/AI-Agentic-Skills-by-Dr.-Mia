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

from urllib.parse import urlparse

from semantic_kernel.utils.feature_stage_decorator import experimental


@experimental
class Uri:
    """The Uri class that represents the URI."""

    def __init__(self, uri):
        """Initialize the Uri."""
        self.uri = uri

    def get_left_part(self):
        """Get the left part of the URI."""
        parsed_uri = urlparse(self.uri)
        return f"{parsed_uri.scheme}://{parsed_uri.netloc}"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
