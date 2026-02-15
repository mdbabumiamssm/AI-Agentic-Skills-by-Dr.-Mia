# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

from langchain_core.utils.function_calling import convert_to_openai_function
from pydantic import BaseModel, Field


def test_convert_pydantic_to_openai_function() -> None:
    class Data(BaseModel):
        """The data to return."""

        key: str = Field(..., description="API key")
        days: int = Field(default=0, description="Number of days to forecast")

    actual = convert_to_openai_function(Data)
    expected = {
        "name": "Data",
        "description": "The data to return.",
        "parameters": {
            "type": "object",
            "properties": {
                "key": {"description": "API key", "type": "string"},
                "days": {
                    "description": "Number of days to forecast",
                    "default": 0,
                    "type": "integer",
                },
            },
            "required": ["key"],
        },
    }
    assert actual == expected


def test_convert_pydantic_to_openai_function_nested() -> None:
    class Data(BaseModel):
        """The data to return."""

        key: str = Field(..., description="API key")
        days: int = Field(default=0, description="Number of days to forecast")

    class Model(BaseModel):
        """The model to return."""

        data: Data

    actual = convert_to_openai_function(Model)
    expected = {
        "name": "Model",
        "description": "The model to return.",
        "parameters": {
            "type": "object",
            "properties": {
                "data": {
                    "description": "The data to return.",
                    "type": "object",
                    "properties": {
                        "key": {
                            "description": "API key",
                            "type": "string",
                        },
                        "days": {
                            "description": "Number of days to forecast",
                            "default": 0,
                            "type": "integer",
                        },
                    },
                    "required": ["key"],
                },
            },
            "required": ["data"],
        },
    }
    assert actual == expected

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
