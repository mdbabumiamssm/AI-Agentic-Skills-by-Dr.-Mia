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

from semantic_kernel.functions import kernel_function
from semantic_kernel.processes.kernel_process import (
    KernelProcessEventVisibility,
    KernelProcessStep,
    KernelProcessStepContext,
)


class ExternalStep(KernelProcessStep):
    external_event_name: str

    def __init__(self, external_event_name: str):
        super().__init__(external_event_name=external_event_name)

    @kernel_function()
    async def emit_external_event(self, context: KernelProcessStepContext, data: Any):
        await context.emit_event(
            process_event=self.external_event_name, data=data, visibility=KernelProcessEventVisibility.Public
        )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
