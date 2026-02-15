# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.output_parsers import BaseOutputParser
from langchain_core.prompts import PromptTemplate
from typing_extensions import override


class FinishedOutputParser(BaseOutputParser[tuple[str, bool]]):
    """Output parser that checks if the output is finished."""

    finished_value: str = "FINISHED"
    """Value that indicates the output is finished."""

    @override
    def parse(self, text: str) -> tuple[str, bool]:
        cleaned = text.strip()
        finished = self.finished_value in cleaned
        return cleaned.replace(self.finished_value, ""), finished


PROMPT_TEMPLATE = """\
Respond to the user message using any relevant context. \
If context is provided, you should ground your answer in that context. \
Once you're done responding return FINISHED.

>>> CONTEXT: {context}
>>> USER INPUT: {user_input}
>>> RESPONSE: {response}\
"""

PROMPT = PromptTemplate(
    template=PROMPT_TEMPLATE,
    input_variables=["user_input", "context", "response"],
)


QUESTION_GENERATOR_PROMPT_TEMPLATE = """\
Given a user input and an existing partial response as context, \
ask a question to which the answer is the given term/entity/phrase:

>>> USER INPUT: {user_input}
>>> EXISTING PARTIAL RESPONSE: {current_response}

The question to which the answer is the term/entity/phrase "{uncertain_span}" is:"""
QUESTION_GENERATOR_PROMPT = PromptTemplate(
    template=QUESTION_GENERATOR_PROMPT_TEMPLATE,
    input_variables=["user_input", "current_response", "uncertain_span"],
)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
