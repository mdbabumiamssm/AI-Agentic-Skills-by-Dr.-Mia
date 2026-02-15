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

from abc import ABC
from typing import TYPE_CHECKING, Generic, TypeVar

from semantic_kernel.kernel_pydantic import KernelBaseModel
from semantic_kernel.utils.feature_stage_decorator import experimental

if TYPE_CHECKING:
    from semantic_kernel.processes.kernel_process.kernel_process_step_state import KernelProcessStepState

TState = TypeVar("TState")


@experimental
class KernelProcessStep(ABC, KernelBaseModel, Generic[TState]):
    """A KernelProcessStep Base class for process steps."""

    state: TState | None = None

    async def activate(self, state: "KernelProcessStepState[TState]"):
        """Activates the step and sets the state."""
        pass  # pragma: no cover

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
