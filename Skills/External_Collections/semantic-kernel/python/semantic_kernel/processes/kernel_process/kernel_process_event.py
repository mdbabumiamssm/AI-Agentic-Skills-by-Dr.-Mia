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

from enum import Enum
from typing import Any

from pydantic import ConfigDict

from semantic_kernel.kernel_pydantic import KernelBaseModel
from semantic_kernel.utils.feature_stage_decorator import experimental


@experimental
class KernelProcessEventVisibility(Enum):
    """Visibility of a kernel process event."""

    # The event is visible inside the process as well as outside the process. This is useful
    # when the event is intended to be consumed by other processes or external systems.
    Public = "Public"

    # The event is only visible to steps within the same process.
    Internal = "Internal"


@experimental
class KernelProcessEvent(KernelBaseModel):
    """A kernel process event."""

    id: str
    data: Any | None = None
    visibility: KernelProcessEventVisibility = KernelProcessEventVisibility.Internal

    model_config = ConfigDict(use_enum_values=False)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
