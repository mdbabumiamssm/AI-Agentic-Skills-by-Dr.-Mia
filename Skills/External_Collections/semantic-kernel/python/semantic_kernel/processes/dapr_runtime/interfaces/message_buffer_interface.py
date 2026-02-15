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

from abc import ABC, abstractmethod

from dapr.actor import ActorInterface, actormethod

from semantic_kernel.utils.feature_stage_decorator import experimental


@experimental
class MessageBufferInterface(ActorInterface, ABC):
    """Abstract base class for a message event buffer that follows the ActorInterface."""

    @abstractmethod
    @actormethod(name="enqueue")
    async def enqueue(self, message: str) -> None:
        """Enqueues a message event into the buffer.

        Args:
            message: The message event to enqueue.
        """
        ...

    @abstractmethod
    @actormethod(name="dequeue_all")
    async def dequeue_all(self) -> list[str]:
        """Dequeues all process events from the buffer.

        Returns:
            The dequeued message event as a list of string
            representing a ProcessEvent.
        """
        ...

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
