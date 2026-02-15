# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import json

import pytest

from langchain_classic.chains import OpenAIModerationChain
from langchain_classic.chains.openai_functions.openapi import get_openapi_chain

api_spec = {
    "openapi": "3.0.0",
    "info": {"title": "JSONPlaceholder API", "version": "1.0.0"},
    "servers": [{"url": "https://jsonplaceholder.typicode.com"}],
    "paths": {
        "/posts": {
            "get": {
                "summary": "Get posts",
                "parameters": [
                    {
                        "name": "_limit",
                        "in": "query",
                        "required": False,
                        "schema": {"type": "integer", "example": 2},
                        "description": "Limit the number of results",
                    },
                ],
            },
        },
    },
}


@pytest.mark.requires("openapi_pydantic")
@pytest.mark.requires("langchain_openai")
def test_openai_openapi_chain() -> None:
    from langchain_openai import ChatOpenAI

    llm = ChatOpenAI(model="gpt-4o-mini", temperature=0)
    chain = get_openapi_chain(json.dumps(api_spec), llm)
    output = chain.invoke({"query": "Fetch the top two posts."})
    assert len(output["response"]) == 2


@pytest.mark.requires("openai")
def test_openai_moderation_chain_instantiation() -> None:
    """Test OpenAIModerationChain."""
    api_key = "foo"

    moderation = OpenAIModerationChain(openai_api_key=api_key)

    assert isinstance(moderation, OpenAIModerationChain)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
