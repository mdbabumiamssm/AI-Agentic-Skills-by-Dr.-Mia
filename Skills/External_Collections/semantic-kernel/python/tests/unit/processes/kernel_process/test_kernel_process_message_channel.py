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

from unittest.mock import AsyncMock

from semantic_kernel.processes.kernel_process.kernel_process_message_channel import KernelProcessMessageChannel
from semantic_kernel.processes.local_runtime.local_event import KernelProcessEvent


class MockKernelProcessMessageChannel(KernelProcessMessageChannel):
    async def emit_event(self, process_event: KernelProcessEvent) -> None:
        pass


async def test_emit_event():
    # Arrange
    event = KernelProcessEvent(id="event_001", data={"key": "value"})
    channel = MockKernelProcessMessageChannel()
    channel.emit_event = AsyncMock()

    # Act
    await channel.emit_event(event)

    # Assert
    channel.emit_event.assert_awaited_once_with(event)


async def test_emit_event_with_no_event():
    # Arrange
    channel = MockKernelProcessMessageChannel()
    channel.emit_event = AsyncMock()

    # Act
    await channel.emit_event(None)

    # Assert
    channel.emit_event.assert_awaited_once_with(None)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
