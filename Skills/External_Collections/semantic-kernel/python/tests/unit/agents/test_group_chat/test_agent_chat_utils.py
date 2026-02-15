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

import base64
from hashlib import sha256

import pytest

from semantic_kernel.agents.group_chat.agent_chat_utils import KeyEncoder
from semantic_kernel.exceptions.agent_exceptions import AgentExecutionException


def test_generate_hash_valid_keys():
    keys = ["key1", "key2", "key3"]
    expected_joined_keys = ":".join(keys).encode("utf-8")
    expected_hash = sha256(expected_joined_keys).digest()
    expected_base64 = base64.b64encode(expected_hash).decode("utf-8")

    result = KeyEncoder.generate_hash(keys)

    assert result == expected_base64


def test_generate_hash_empty_keys():
    with pytest.raises(
        AgentExecutionException, match="Channel Keys must not be empty. Unable to generate channel hash."
    ):
        KeyEncoder.generate_hash([])

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
