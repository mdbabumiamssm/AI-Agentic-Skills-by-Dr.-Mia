# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""In memory store that is not thread safe and has no eviction policy.

This is a simple implementation of the BaseStore using a dictionary that is useful
primarily for unit testing purposes.
"""

from langchain_core.stores import InMemoryBaseStore, InMemoryByteStore, InMemoryStore

__all__ = [
    "InMemoryBaseStore",
    "InMemoryByteStore",
    "InMemoryStore",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
