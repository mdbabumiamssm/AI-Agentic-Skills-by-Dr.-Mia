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

from langchain_core.load import dump, loads
from langchain_core.prompts import ChatPromptTemplate


def test_image_prompt_template_deserializable() -> None:
    """Test that the image prompt template is serializable."""
    loads(
        dump.dumps(
            ChatPromptTemplate.from_messages(
                [("system", [{"type": "image", "image_url": "{img}"}])]
            )
        ),
    )


def test_image_prompt_template_deserializable_old() -> None:
    """Test that the image prompt template is serializable."""
    loads(
        json.dumps(
            {
                "lc": 1,
                "type": "constructor",
                "id": ["langchain", "prompts", "chat", "ChatPromptTemplate"],
                "kwargs": {
                    "messages": [
                        {
                            "lc": 1,
                            "type": "constructor",
                            "id": [
                                "langchain",
                                "prompts",
                                "chat",
                                "SystemMessagePromptTemplate",
                            ],
                            "kwargs": {
                                "prompt": [
                                    {
                                        "lc": 1,
                                        "type": "constructor",
                                        "id": [
                                            "langchain",
                                            "prompts",
                                            "prompt",
                                            "PromptTemplate",
                                        ],
                                        "kwargs": {
                                            "template": "Foo",
                                            "input_variables": [],
                                            "template_format": "f-string",
                                            "partial_variables": {},
                                        },
                                    }
                                ]
                            },
                        },
                        {
                            "lc": 1,
                            "type": "constructor",
                            "id": [
                                "langchain",
                                "prompts",
                                "chat",
                                "HumanMessagePromptTemplate",
                            ],
                            "kwargs": {
                                "prompt": [
                                    {
                                        "lc": 1,
                                        "type": "constructor",
                                        "id": [
                                            "langchain",
                                            "prompts",
                                            "image",
                                            "ImagePromptTemplate",
                                        ],
                                        "kwargs": {
                                            "template": {
                                                "url": "data:image/png;base64,{img}"
                                            },
                                            "input_variables": ["img"],
                                        },
                                    },
                                    {
                                        "lc": 1,
                                        "type": "constructor",
                                        "id": [
                                            "langchain",
                                            "prompts",
                                            "prompt",
                                            "PromptTemplate",
                                        ],
                                        "kwargs": {
                                            "template": "{input}",
                                            "input_variables": ["input"],
                                            "template_format": "f-string",
                                            "partial_variables": {},
                                        },
                                    },
                                ]
                            },
                        },
                    ],
                    "input_variables": ["img", "input"],
                },
            }
        ),
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
