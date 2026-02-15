# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""LLM Chain for generating examples for question answering."""

from __future__ import annotations

from typing import Any

from langchain_core.language_models import BaseLanguageModel
from langchain_core.output_parsers import BaseLLMOutputParser
from pydantic import Field
from typing_extensions import override

from langchain_classic.chains.llm import LLMChain
from langchain_classic.evaluation.qa.generate_prompt import PROMPT
from langchain_classic.output_parsers.regex import RegexParser

_QA_OUTPUT_PARSER = RegexParser(
    regex=r"QUESTION: (.*?)\n+ANSWER: (.*)",
    output_keys=["query", "answer"],
)


class QAGenerateChain(LLMChain):
    """LLM Chain for generating examples for question answering."""

    output_parser: BaseLLMOutputParser = Field(default=_QA_OUTPUT_PARSER)
    output_key: str = "qa_pairs"

    @classmethod
    @override
    def is_lc_serializable(cls) -> bool:
        return False

    @classmethod
    def from_llm(cls, llm: BaseLanguageModel, **kwargs: Any) -> QAGenerateChain:
        """Load QA Generate Chain from LLM."""
        return cls(llm=llm, prompt=PROMPT, **kwargs)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
