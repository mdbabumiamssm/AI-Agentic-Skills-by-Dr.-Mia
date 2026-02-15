# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Test LLM Math functionality."""

import pytest

from langchain_classic.chains.llm_math.base import LLMMathChain
from langchain_classic.chains.llm_math.prompt import _PROMPT_TEMPLATE
from tests.unit_tests.llms.fake_llm import FakeLLM


@pytest.fixture
def fake_llm_math_chain() -> LLMMathChain:
    """Fake LLM Math chain for testing."""
    complex_question = _PROMPT_TEMPLATE.format(question="What is the square root of 2?")
    queries = {
        _PROMPT_TEMPLATE.format(question="What is 1 plus 1?"): "Answer: 2",
        complex_question: "```text\n2**.5\n```",
        _PROMPT_TEMPLATE.format(question="foo"): "foo",
    }
    fake_llm = FakeLLM(queries=queries)
    return LLMMathChain.from_llm(fake_llm, input_key="q", output_key="a")


@pytest.mark.requires("numexpr")
def test_simple_question(fake_llm_math_chain: LLMMathChain) -> None:
    """Test simple question that should not need python."""
    question = "What is 1 plus 1?"
    output = fake_llm_math_chain.run(question)
    assert output == "Answer: 2"


@pytest.mark.requires("numexpr")
def test_complex_question(fake_llm_math_chain: LLMMathChain) -> None:
    """Test complex question that should need python."""
    question = "What is the square root of 2?"
    output = fake_llm_math_chain.run(question)
    assert output == f"Answer: {2**0.5}"


@pytest.mark.requires("numexpr")
def test_error(fake_llm_math_chain: LLMMathChain) -> None:
    """Test question that raises error."""
    with pytest.raises(ValueError, match="unknown format from LLM: foo"):
        fake_llm_math_chain.run("foo")

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
