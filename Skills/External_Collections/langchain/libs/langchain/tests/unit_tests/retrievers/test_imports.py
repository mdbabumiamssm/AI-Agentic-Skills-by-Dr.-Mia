# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_classic import retrievers

EXPECTED_ALL = [
    "AmazonKendraRetriever",
    "AmazonKnowledgeBasesRetriever",
    "ArceeRetriever",
    "ArxivRetriever",
    "AzureAISearchRetriever",
    "AzureCognitiveSearchRetriever",
    "BM25Retriever",
    "ChaindeskRetriever",
    "ChatGPTPluginRetriever",
    "CohereRagRetriever",
    "ContextualCompressionRetriever",
    "DocArrayRetriever",
    "DriaRetriever",
    "ElasticSearchBM25Retriever",
    "EmbedchainRetriever",
    "EnsembleRetriever",
    "GoogleCloudEnterpriseSearchRetriever",
    "GoogleDocumentAIWarehouseRetriever",
    "GoogleVertexAIMultiTurnSearchRetriever",
    "GoogleVertexAISearchRetriever",
    "KayAiRetriever",
    "KNNRetriever",
    "LlamaIndexGraphRetriever",
    "LlamaIndexRetriever",
    "MergerRetriever",
    "MetalRetriever",
    "MilvusRetriever",
    "MultiQueryRetriever",
    "MultiVectorRetriever",
    "NeuralDBRetriever",
    "OutlineRetriever",
    "ParentDocumentRetriever",
    "PineconeHybridSearchRetriever",
    "PubMedRetriever",
    "RemoteLangChainRetriever",
    "RePhraseQueryRetriever",
    "SelfQueryRetriever",
    "SVMRetriever",
    "TavilySearchAPIRetriever",
    "TFIDFRetriever",
    "TimeWeightedVectorStoreRetriever",
    "VespaRetriever",
    "WeaviateHybridSearchRetriever",
    "WebResearchRetriever",
    "WikipediaRetriever",
    "ZepRetriever",
    "ZillizRetriever",
]


def test_imports() -> None:
    assert sorted(retrievers.__all__) == sorted(EXPECTED_ALL)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
