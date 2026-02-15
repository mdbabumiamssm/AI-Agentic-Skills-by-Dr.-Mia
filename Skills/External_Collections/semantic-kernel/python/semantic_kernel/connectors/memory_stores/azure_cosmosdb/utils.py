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

from enum import Enum

from semantic_kernel.utils.feature_stage_decorator import experimental


@experimental
class CosmosDBSimilarityType(str, Enum):
    """Cosmos DB Similarity Type as enumerator."""

    COS = "COS"
    """CosineSimilarity"""
    IP = "IP"
    """inner - product"""
    L2 = "L2"
    """Euclidean distance"""


@experimental
class CosmosDBVectorSearchType(str, Enum):
    """Cosmos DB Vector Search Type as enumerator."""

    VECTOR_IVF = "vector-ivf"
    """IVF vector index"""
    VECTOR_HNSW = "vector-hnsw"
    """HNSW vector index"""

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
