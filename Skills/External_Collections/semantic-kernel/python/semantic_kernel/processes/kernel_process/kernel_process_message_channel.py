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

from semantic_kernel.processes.local_runtime.local_event import KernelProcessEvent
from semantic_kernel.utils.feature_stage_decorator import experimental


@experimental
class KernelProcessMessageChannel(ABC):
    """Abstract base class for emitting events from a step."""

    @abstractmethod
    async def emit_event(self, process_event: "KernelProcessEvent") -> None:
        """Emits the specified event from the step."""
        pass

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
