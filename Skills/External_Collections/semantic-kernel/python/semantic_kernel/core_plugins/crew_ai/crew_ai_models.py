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

from semantic_kernel.kernel_pydantic import KernelBaseModel


class CrewAIEnterpriseKickoffState(str, Enum):
    """The Crew.AI Enterprise kickoff state."""

    Pending = "PENDING"
    Started = "STARTED"
    Running = "RUNNING"
    Success = "SUCCESS"
    Failed = "FAILED"
    Failure = "FAILURE"
    Not_Found = "NOT FOUND"


class CrewAIStatusResponse(KernelBaseModel):
    """Represents the status response from Crew AI."""

    state: CrewAIEnterpriseKickoffState
    result: str | None = None
    last_step: dict[str, Any] | None = None


class CrewAIKickoffResponse(KernelBaseModel):
    """Represents the kickoff response from Crew AI."""

    kickoff_id: str


class CrewAIRequiredInputs(KernelBaseModel):
    """Represents the required inputs for Crew AI."""

    inputs: dict[str, str]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
