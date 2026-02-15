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

from langchain_classic.retrievers.document_compressors.listwise_rerank import (
    LLMListwiseRerank,
)


def test_list_rerank() -> None:
    from langchain_openai import ChatOpenAI

    documents = [
        Document("Sally is my friend from school"),
        Document("Steve is my friend from home"),
        Document("I didn't always like yogurt"),
        Document("I wonder why it's called football"),
        Document("Where's waldo"),
    ]

    reranker = LLMListwiseRerank.from_llm(
        llm=ChatOpenAI(model="gpt-3.5-turbo"),
        top_n=3,
    )
    compressed_docs = reranker.compress_documents(documents, "Who is steve")
    assert len(compressed_docs) == 3
    assert "Steve" in compressed_docs[0].page_content

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
