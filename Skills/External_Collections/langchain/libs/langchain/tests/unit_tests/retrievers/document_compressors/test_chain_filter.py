# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.documents import Document
from langchain_core.language_models import FakeListChatModel

from langchain_classic.retrievers.document_compressors import LLMChainFilter


def test_llm_chain_filter() -> None:
    documents = [
        Document(
            page_content="Candlepin bowling is popular in New England.",
            metadata={"a": 1},
        ),
        Document(
            page_content="Candlepin bowling balls are smaller.",
            metadata={"b": 2},
        ),
        Document(page_content="The moon is round.", metadata={"c": 3}),
    ]
    llm = FakeListChatModel(responses=["YES", "YES", "NO"])
    doc_compressor = LLMChainFilter.from_llm(llm)
    output = doc_compressor.compress_documents(
        documents,
        "Tell me about Candlepin bowling.",
    )
    expected = documents[:2]
    assert output == expected


async def test_llm_chain_extractor_async() -> None:
    documents = [
        Document(
            page_content="Candlepin bowling is popular in New England.",
            metadata={"a": 1},
        ),
        Document(
            page_content="Candlepin bowling balls are smaller.",
            metadata={"b": 2},
        ),
        Document(page_content="The moon is round.", metadata={"c": 3}),
    ]
    llm = FakeListChatModel(responses=["YES", "YES", "NO"])
    doc_compressor = LLMChainFilter.from_llm(llm)
    output = await doc_compressor.acompress_documents(
        documents,
        "Tell me about Candlepin bowling.",
    )
    expected = documents[:2]
    assert output == expected

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
