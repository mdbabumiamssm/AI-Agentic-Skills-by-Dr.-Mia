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

import copy
import os
from collections.abc import Awaitable, Callable
from typing import Any

import pytest
from pytest import MonkeyPatch, mark, param

from samples.getting_started.threads.custom_chat_message_store_thread import main as threads_custom_store
from samples.getting_started.threads.suspend_resume_thread import main as threads_suspend_resume

# Environment variable for controlling sample tests
RUN_SAMPLES_TESTS = "RUN_SAMPLES_TESTS"

# All thread samples
thread_samples = [
    param(
        threads_custom_store,
        [],  # Non-interactive sample
        id="threads_custom_store",
        marks=[
            pytest.mark.openai,
            pytest.mark.skipif(os.getenv(RUN_SAMPLES_TESTS, None) is None, reason="Not running sample tests."),
        ],
    ),
    param(
        threads_suspend_resume,
        [],  # Non-interactive sample
        id="threads_suspend_resume",
        marks=[
            pytest.mark.openai,
            pytest.mark.skipif(os.getenv(RUN_SAMPLES_TESTS, None) is None, reason="Not running sample tests."),
        ],
    ),
]


@mark.parametrize("sample, responses", thread_samples)
async def test_thread_samples(sample: Callable[..., Awaitable[Any]], responses: list[str], monkeypatch: MonkeyPatch):
    """Test thread samples with input mocking and retry logic."""
    saved_responses = copy.deepcopy(responses)

    def reset():
        responses.clear()
        responses.extend(saved_responses)

    def mock_input(prompt: str = "") -> str:
        return responses.pop(0) if responses else "exit"

    monkeypatch.setattr("builtins.input", mock_input)
    await sample

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
