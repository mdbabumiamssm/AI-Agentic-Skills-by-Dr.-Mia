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
class ExternalEventBufferInterface(ActorInterface, ABC):
    """Abstract base class for an external event buffer that follows the ActorInterface."""

    @abstractmethod
    @actormethod(name="enqueue")
    async def enqueue(self, external_event: str) -> None:
        """Enqueues an external event into the buffer.

        Args:
            external_event: The external event to enqueue.
        """
        ...

    @abstractmethod
    @actormethod(name="dequeue_all")
    async def dequeue_all(self) -> list[str]:
        """Dequeues all external events from the buffer as.

        The list is of string representations of KernelProcessEvent.

        Returns:
            The dequeued external event.
        """
        ...

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
