# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.prompts.prompt import PromptTemplate

_CREATE_DRAFT_ANSWER_TEMPLATE = """{question}\n\n"""
CREATE_DRAFT_ANSWER_PROMPT = PromptTemplate(
    input_variables=["question"], template=_CREATE_DRAFT_ANSWER_TEMPLATE
)

_LIST_ASSERTIONS_TEMPLATE = """Here is a statement:
{statement}
Make a bullet point list of the assumptions you made when producing the above statement.\n\n"""  # noqa: E501
LIST_ASSERTIONS_PROMPT = PromptTemplate(
    input_variables=["statement"], template=_LIST_ASSERTIONS_TEMPLATE
)

_CHECK_ASSERTIONS_TEMPLATE = """Here is a bullet point list of assertions:
{assertions}
For each assertion, determine whether it is true or false. If it is false, explain why.\n\n"""  # noqa: E501
CHECK_ASSERTIONS_PROMPT = PromptTemplate(
    input_variables=["assertions"], template=_CHECK_ASSERTIONS_TEMPLATE
)

_REVISED_ANSWER_TEMPLATE = """{checked_assertions}

Question: In light of the above assertions and checks, how would you answer the question '{question}'?

Answer:"""  # noqa: E501
REVISED_ANSWER_PROMPT = PromptTemplate(
    input_variables=["checked_assertions", "question"],
    template=_REVISED_ANSWER_TEMPLATE,
)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
