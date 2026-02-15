# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Embeddings models.

!!! warning "Modules moved"
    With the release of `langchain 1.0.0`, several embeddings modules were moved to
    `langchain-classic`, such as `CacheBackedEmbeddings` and all community
    embeddings. See [list](https://github.com/langchain-ai/langchain/blob/bdf1cd383ce36dc18381a3bf3fb0a579337a32b5/libs/langchain/langchain/embeddings/__init__.py)
    of moved modules to inform your migration.
"""

from langchain_core.embeddings import Embeddings

from langchain.embeddings.base import init_embeddings

__all__ = [
    "Embeddings",
    "init_embeddings",
]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
