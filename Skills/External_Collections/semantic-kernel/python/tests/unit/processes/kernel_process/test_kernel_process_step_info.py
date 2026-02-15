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

from unittest.mock import MagicMock

import pytest
from pydantic import ValidationError

from semantic_kernel.processes.kernel_process.kernel_process_edge import KernelProcessEdge
from semantic_kernel.processes.kernel_process.kernel_process_step_info import KernelProcessStepInfo
from semantic_kernel.processes.kernel_process.kernel_process_step_state import KernelProcessStepState


class MockStepType:
    """A mock class to simulate a step type."""

    pass


def test_initialization():
    # Arrange
    inner_step_type = MockStepType
    state = MagicMock(spec=KernelProcessStepState)
    edge = MagicMock(spec=KernelProcessEdge)
    edges = {"output": [edge]}

    # Act
    step_info = KernelProcessStepInfo(inner_step_type=inner_step_type, state=state, output_edges=edges)

    # Assert
    assert step_info.inner_step_type == inner_step_type
    assert step_info.state == state
    assert step_info.output_edges == edges


def test_initialization_with_missing_arguments():
    # Arrange
    inner_step_type = MockStepType
    state = MagicMock(spec=KernelProcessStepState)
    edge = MagicMock(spec=KernelProcessEdge)
    edges = {"output": [edge]}

    # Act & Assert
    with pytest.raises(ValidationError):
        KernelProcessStepInfo(inner_step_type=None, state=state, output_edges=edges)

    with pytest.raises(ValidationError):
        KernelProcessStepInfo(inner_step_type=inner_step_type, state=None, output_edges=edges)

    with pytest.raises(ValidationError):
        KernelProcessStepInfo(inner_step_type=inner_step_type, state=state, output_edges=None)


def test_edges_property():
    # Arrange
    inner_step_type = MockStepType
    state = MagicMock(spec=KernelProcessStepState)
    edge_1 = MagicMock(spec=KernelProcessEdge)
    edge_2 = MagicMock(spec=KernelProcessEdge)
    edges = {"output": [edge_1, edge_2]}

    # Act
    step_info = KernelProcessStepInfo(inner_step_type=inner_step_type, state=state, output_edges=edges)
    copied_edges = step_info.edges

    # Assert
    assert copied_edges == edges
    assert copied_edges is not edges
    assert copied_edges["output"] is not edges["output"]

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
