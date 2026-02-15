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
from langchain_tests.unit_tests import (  # type: ignore[import-not-found]
    ChatModelUnitTests,  # type: ignore[import-not-found]
)

from langchain_fireworks import ChatFireworks


class TestFireworksStandard(ChatModelUnitTests):
    @property
    def chat_model_class(self) -> type[BaseChatModel]:
        return ChatFireworks

    @property
    def chat_model_params(self) -> dict:
        return {
            "model": "accounts/fireworks/models/llama-v3p1-70b-instruct",
            "api_key": "test_api_key",
        }

    @property
    def init_from_env_params(self) -> tuple[dict, dict, dict]:
        return (
            {
                "FIREWORKS_API_KEY": "api_key",
                "FIREWORKS_API_BASE": "https://base.com",
            },
            {
                "model": "accounts/fireworks/models/llama-v3p1-70b-instruct",
            },
            {
                "fireworks_api_key": "api_key",
                "fireworks_api_base": "https://base.com",
            },
        )


def test_profile() -> None:
    """Test that model profile is loaded correctly."""
    model = ChatFireworks(
        model="accounts/fireworks/models/gpt-oss-20b",
        api_key="test_key",  # type: ignore[arg-type]
    )
    assert model.profile

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
