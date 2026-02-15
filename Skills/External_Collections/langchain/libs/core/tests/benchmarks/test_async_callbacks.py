# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

import asyncio
from itertools import cycle
from typing import Any
from uuid import UUID

import pytest
from pytest_benchmark.fixture import BenchmarkFixture
from typing_extensions import override

from langchain_core.callbacks.base import AsyncCallbackHandler
from langchain_core.language_models import GenericFakeChatModel
from langchain_core.messages import AIMessage, BaseMessage
from langchain_core.outputs import ChatGenerationChunk, GenerationChunk


class MyCustomAsyncHandler(AsyncCallbackHandler):
    @override
    async def on_chat_model_start(
        self,
        serialized: dict[str, Any],
        messages: list[list[BaseMessage]],
        *,
        run_id: UUID,
        parent_run_id: UUID | None = None,
        tags: list[str] | None = None,
        metadata: dict[str, Any] | None = None,
        **kwargs: Any,
    ) -> Any:
        # Do nothing
        # Required to implement since this is an abstract method
        pass

    @override
    async def on_llm_new_token(
        self,
        token: str,
        *,
        chunk: GenerationChunk | ChatGenerationChunk | None = None,
        run_id: UUID,
        parent_run_id: UUID | None = None,
        tags: list[str] | None = None,
        **kwargs: Any,
    ) -> None:
        await asyncio.sleep(0)


@pytest.mark.benchmark
async def test_async_callbacks_in_sync(benchmark: BenchmarkFixture) -> None:
    infinite_cycle = cycle([AIMessage(content=" ".join(["hello", "goodbye"] * 5))])
    model = GenericFakeChatModel(messages=infinite_cycle)

    @benchmark  # type: ignore[misc]
    def sync_callbacks() -> None:
        for _ in range(5):
            for _ in model.stream("meow", {"callbacks": [MyCustomAsyncHandler()]}):
                pass

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
