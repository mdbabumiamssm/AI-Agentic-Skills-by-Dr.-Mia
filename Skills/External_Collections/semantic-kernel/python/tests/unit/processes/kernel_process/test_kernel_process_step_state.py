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

import pytest
from pydantic import ValidationError

from semantic_kernel.processes.kernel_process.kernel_process_step_state import KernelProcessStepState


def test_initialization_with_name_and_id():
    # Arrange
    name = "step_state"
    step_id = "5678"

    # Act
    step_state = KernelProcessStepState(name=name, id=step_id, version="1.0")

    # Assert
    assert step_state.name == name
    assert step_state.id == step_id
    assert step_state.state is None


def test_initialization_with_name_only():
    # Arrange
    name = "step_state_without_id"

    # Act
    step_state = KernelProcessStepState(name=name, version="1.0")

    # Assert
    assert step_state.name == name
    assert step_state.id is None
    assert step_state.state is None


def test_setting_step_state_value():
    # Arrange
    name = "step_state"
    state_value = {"status": "in_progress"}

    # Act
    step_state = KernelProcessStepState(name=name, version="1.0")
    step_state.state = state_value

    # Assert
    assert step_state.state == state_value


def test_initialization_with_invalid_name():
    # Arrange
    name = 12345  # Invalid type for name

    # Act & Assert
    with pytest.raises(ValidationError):
        KernelProcessStepState(name=name)

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
