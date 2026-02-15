# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import pytest
from packaging import version

from langchain_core.prompts.string import mustache_schema
from langchain_core.utils.pydantic import PYDANTIC_VERSION

PYDANTIC_VERSION_AT_LEAST_29 = version.parse("2.9") <= PYDANTIC_VERSION


@pytest.mark.skipif(
    not PYDANTIC_VERSION_AT_LEAST_29,
    reason=(
        "Only test with most recent version of pydantic. "
        "Pydantic introduced small fixes to generated JSONSchema on minor versions."
    ),
)
def test_mustache_schema_parent_child() -> None:
    template = "{{x.y}} {{x}}"
    expected = {
        "$defs": {
            "x": {
                "properties": {"y": {"default": None, "title": "Y", "type": "string"}},
                "title": "x",
                "type": "object",
            }
        },
        "properties": {"x": {"$ref": "#/$defs/x", "default": None}},
        "title": "PromptInput",
        "type": "object",
    }
    actual = mustache_schema(template).model_json_schema()
    assert expected == actual

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
