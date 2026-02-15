# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

# Copyright (c) Microsoft. All rights reserved.

from pydantic import ValidationError
from pytest import raises

from semantic_kernel.template_engine.blocks.block import Block


def test_init():
    block = Block(content="test content")
    assert block.content == "test content"


def test_content_strip():
    block = Block(content=" test content ")
    assert block.content == "test content"


def test_no_content():
    with raises(ValidationError):
        Block()

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
