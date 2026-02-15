# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Base abstraction and in-memory implementation of rate limiters.

These rate limiters can be used to limit the rate of requests to an API.

The rate limiters can be used together with `BaseChatModel`.
"""

from langchain_core.rate_limiters import BaseRateLimiter, InMemoryRateLimiter

__all__ = [
    "BaseRateLimiter",
    "InMemoryRateLimiter",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
