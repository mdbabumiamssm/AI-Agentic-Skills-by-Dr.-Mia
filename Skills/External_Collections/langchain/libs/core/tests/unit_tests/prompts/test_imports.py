# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.prompts import __all__

EXPECTED_ALL = [
    "AIMessagePromptTemplate",
    "BaseChatPromptTemplate",
    "BasePromptTemplate",
    "ChatMessagePromptTemplate",
    "ChatPromptTemplate",
    "DictPromptTemplate",
    "FewShotPromptTemplate",
    "FewShotPromptWithTemplates",
    "FewShotChatMessagePromptTemplate",
    "format_document",
    "aformat_document",
    "HumanMessagePromptTemplate",
    "MessagesPlaceholder",
    "PromptTemplate",
    "StringPromptTemplate",
    "SystemMessagePromptTemplate",
    "load_prompt",
    "check_valid_template",
    "get_template_variables",
    "jinja2_formatter",
    "validate_jinja2",
]


def test_all_imports() -> None:
    assert set(__all__) == set(EXPECTED_ALL)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
