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

import pytest

from semantic_kernel.contents.chat_history import ChatHistory


@pytest.fixture()
def service_id() -> str:
    return "test_service_id"


@pytest.fixture()
def chat_history() -> ChatHistory:
    chat_history = ChatHistory()
    chat_history.add_system_message("system_prompt")
    chat_history.add_user_message("test_prompt")
    return chat_history


@pytest.fixture()
def prompt() -> str:
    return "test_prompt"

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
