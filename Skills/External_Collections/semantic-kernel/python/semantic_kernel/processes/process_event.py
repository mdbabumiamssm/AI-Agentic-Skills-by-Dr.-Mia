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

from typing import Any

from semantic_kernel.kernel_pydantic import KernelBaseModel
from semantic_kernel.processes.kernel_process.kernel_process_event import (
    KernelProcessEvent,
    KernelProcessEventVisibility,
)
from semantic_kernel.utils.feature_stage_decorator import experimental


@experimental
class ProcessEvent(KernelBaseModel):
    """A wrapper around KernelProcessEvent that helps to manage the namespace of the event."""

    namespace: str | None = None
    inner_event: KernelProcessEvent

    @property
    def id(self) -> str:
        """The Id of the event."""
        return f"{self.namespace}.{self.inner_event.id}"

    @property
    def data(self) -> Any | None:
        """The data of the event."""
        return self.inner_event.data

    @property
    def visibility(self) -> "KernelProcessEventVisibility":
        """The visibility of the event."""
        return self.inner_event.visibility

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
