# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test LLMSummarization functionality."""

import pytest

from langchain_classic.chains.llm_summarization_checker.base import (
    ARE_ALL_TRUE_PROMPT,
    CHECK_ASSERTIONS_PROMPT,
    CREATE_ASSERTIONS_PROMPT,
    REVISED_SUMMARY_PROMPT,
    LLMSummarizationCheckerChain,
)
from tests.unit_tests.llms.fake_llm import FakeLLM


def test_input_variables() -> None:
    assert CREATE_ASSERTIONS_PROMPT.input_variables == ["summary"]
    assert CHECK_ASSERTIONS_PROMPT.input_variables == ["assertions"]
    assert REVISED_SUMMARY_PROMPT.input_variables == ["checked_assertions", "summary"]
    assert ARE_ALL_TRUE_PROMPT.input_variables == ["checked_assertions"]


@pytest.fixture
def fake_llm_summarization_checker_chain() -> LLMSummarizationCheckerChain:
    """Fake LLMCheckerChain for testing."""
    queries = {
        CREATE_ASSERTIONS_PROMPT.format(
            summary="a",
        ): "b",
        CHECK_ASSERTIONS_PROMPT.format(
            assertions="b",
        ): "- b - True",
        REVISED_SUMMARY_PROMPT.format(
            checked_assertions="- b - True", summary="a"
        ): "b",
        ARE_ALL_TRUE_PROMPT.format(
            checked_assertions="- b - True",
        ): "True",
    }
    fake_llm = FakeLLM(queries=queries)
    return LLMSummarizationCheckerChain.from_llm(
        fake_llm, input_key="q", output_key="a"
    )


def test_simple_text(
    fake_llm_summarization_checker_chain: LLMSummarizationCheckerChain,
) -> None:
    """Test simple question that should not need python."""
    question = "a"
    output = fake_llm_summarization_checker_chain.run(question)
    assert output == "b"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
