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
from langchain_core.language_models import FakeListLLM
from langchain_core.prompts import PromptTemplate

from langchain_classic.chains import create_history_aware_retriever
from tests.unit_tests.retrievers.parrot_retriever import FakeParrotRetriever


def test_create() -> None:
    answer = "I know the answer!"
    llm = FakeListLLM(responses=[answer])
    retriever = FakeParrotRetriever()
    question_gen_prompt = PromptTemplate.from_template("hi! {input} {chat_history}")
    chain = create_history_aware_retriever(llm, retriever, question_gen_prompt)
    expected_output = [Document(page_content="What is the answer?")]
    output = chain.invoke({"input": "What is the answer?", "chat_history": []})
    assert output == expected_output

    output = chain.invoke({"input": "What is the answer?"})
    assert output == expected_output

    expected_output = [Document(page_content="I know the answer!")]
    output = chain.invoke(
        {
            "input": "What is the answer?",
            "chat_history": ["hi", "hi"],
        },
    )
    assert output == expected_output

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
