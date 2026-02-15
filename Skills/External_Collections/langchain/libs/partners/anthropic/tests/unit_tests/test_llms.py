# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import os

from langchain_anthropic import AnthropicLLM

os.environ["ANTHROPIC_API_KEY"] = "foo"


def test_anthropic_model_params() -> None:
    # Test standard tracing params
    llm = AnthropicLLM(model="foo")  # type: ignore[call-arg]

    ls_params = llm._get_ls_params()
    assert ls_params == {
        "ls_provider": "anthropic",
        "ls_model_type": "llm",
        "ls_model_name": "foo",
        "ls_max_tokens": 1024,
    }

    llm = AnthropicLLM(model="foo", temperature=0.1)  # type: ignore[call-arg]

    ls_params = llm._get_ls_params()
    assert ls_params == {
        "ls_provider": "anthropic",
        "ls_model_type": "llm",
        "ls_model_name": "foo",
        "ls_max_tokens": 1024,
        "ls_temperature": 0.1,
    }

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
