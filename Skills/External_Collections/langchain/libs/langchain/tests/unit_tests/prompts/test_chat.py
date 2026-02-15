# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_classic.prompts.chat import __all__

EXPECTED_ALL = [
    "MessageLike",
    "MessageLikeRepresentation",
    "MessagePromptTemplateT",
    "AIMessagePromptTemplate",
    "BaseChatPromptTemplate",
    "BaseMessagePromptTemplate",
    "BaseStringMessagePromptTemplate",
    "ChatMessagePromptTemplate",
    "ChatPromptTemplate",
    "ChatPromptValue",
    "ChatPromptValueConcrete",
    "HumanMessagePromptTemplate",
    "MessagesPlaceholder",
    "SystemMessagePromptTemplate",
    "_convert_to_message",
    "_create_template_from_message_type",
]


def test_all_imports() -> None:
    assert set(__all__) == set(EXPECTED_ALL)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
