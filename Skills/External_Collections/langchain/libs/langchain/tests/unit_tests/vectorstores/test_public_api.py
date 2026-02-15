# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test the public API of the tools package."""

from langchain_classic.vectorstores import __all__ as public_api

_EXPECTED = [
    "AlibabaCloudOpenSearch",
    "AlibabaCloudOpenSearchSettings",
    "AnalyticDB",
    "Annoy",
    "AstraDB",
    "AtlasDB",
    "AwaDB",
    "AzureCosmosDBVectorSearch",
    "AzureSearch",
    "Bagel",
    "Cassandra",
    "Chroma",
    "Clarifai",
    "Clickhouse",
    "ClickhouseSettings",
    "DashVector",
    "DatabricksVectorSearch",
    "DeepLake",
    "Dingo",
    "DocArrayHnswSearch",
    "DocArrayInMemorySearch",
    "DuckDB",
    "EcloudESVectorStore",
    "ElasticKnnSearch",
    "ElasticsearchStore",
    "ElasticVectorSearch",
    "Epsilla",
    "FAISS",
    "Hologres",
    "LanceDB",
    "LLMRails",
    "Marqo",
    "MatchingEngine",
    "Meilisearch",
    "Milvus",
    "MomentoVectorIndex",
    "MongoDBAtlasVectorSearch",
    "MyScale",
    "MyScaleSettings",
    "Neo4jVector",
    "NeuralDBClientVectorStore",
    "NeuralDBVectorStore",
    "OpenSearchVectorSearch",
    "PGEmbedding",
    "PGVector",
    "Pinecone",
    "Qdrant",
    "Redis",
    "Rockset",
    "ScaNN",
    "SemaDB",
    "SingleStoreDB",
    "SKLearnVectorStore",
    "SQLiteVSS",
    "StarRocks",
    "SupabaseVectorStore",
    "Tair",
    "TencentVectorDB",
    "TileDB",
    "TimescaleVector",
    "Typesense",
    "USearch",
    "Vald",
    "Vearch",
    "Vectara",
    "VectorStore",
    "VespaStore",
    "Weaviate",
    "Yellowbrick",
    "ZepVectorStore",
    "Zilliz",
]


def test_public_api() -> None:
    """Test for regressions or changes in the public API."""
    # Check that the public API is as expected
    assert set(public_api) == set(_EXPECTED)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
