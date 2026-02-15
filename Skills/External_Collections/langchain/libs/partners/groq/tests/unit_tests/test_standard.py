# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

"""Standard LangChain interface tests."""

from langchain_core.language_models import BaseChatModel
from langchain_tests.unit_tests.chat_models import (
    ChatModelUnitTests,
)

from langchain_groq import ChatGroq


class TestGroqStandard(ChatModelUnitTests):
    """Run ChatGroq on LangChain standard tests."""

    @property
    def chat_model_class(self) -> type[BaseChatModel]:
        return ChatGroq

    @property
    def chat_model_params(self) -> dict:
        return {"model": "llama-3.1-8b-instant"}

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
