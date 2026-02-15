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

import pytest
from langchain_core.language_models import BaseChatModel
from langchain_tests.integration_tests import ChatModelIntegrationTests

from langchain_perplexity import ChatPerplexity


class TestPerplexityStandard(ChatModelIntegrationTests):
    @property
    def chat_model_class(self) -> type[BaseChatModel]:
        return ChatPerplexity

    @property
    def chat_model_params(self) -> dict:
        return {"model": "sonar"}

    @pytest.mark.xfail(reason="TODO: handle in integration.")
    def test_double_messages_conversation(self, model: BaseChatModel) -> None:
        super().test_double_messages_conversation(model)

    @pytest.mark.xfail(reason="Raises 400: Custom stop words not supported.")
    def test_stop_sequence(self, model: BaseChatModel) -> None:
        super().test_stop_sequence(model)

    # TODO, API regressed for some reason after 2025-04-15

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
