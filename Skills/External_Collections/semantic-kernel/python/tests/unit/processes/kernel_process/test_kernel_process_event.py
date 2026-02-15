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

from semantic_kernel.processes.kernel_process.kernel_process_event import (
    KernelProcessEvent,
    KernelProcessEventVisibility,
)


def test_initialization():
    # Arrange
    event_id = "event_001"
    event_data = {"key": "value"}

    # Act
    event = KernelProcessEvent(id=event_id, data=event_data)

    # Assert
    assert event.id == event_id
    assert event.data == event_data
    assert event.visibility == KernelProcessEventVisibility.Internal


def test_initialization_with_visibility():
    # Arrange
    event_id = "event_002"
    event_data = [1, 2, 3]
    visibility = KernelProcessEventVisibility.Public

    # Act
    event = KernelProcessEvent(id=event_id, data=event_data, visibility=visibility)

    # Assert
    assert event.id == event_id
    assert event.data == event_data
    assert event.visibility == KernelProcessEventVisibility.Public


def test_invalid_visibility():
    # Arrange
    event_id = "event_003"
    event_data = "sample data"
    invalid_visibility = "Hidden"

    # Act & Assert
    with pytest.raises(ValueError):
        KernelProcessEvent(id=event_id, data=event_data, visibility=invalid_visibility)


def test_none_data_initialization():
    # Arrange
    event_id = "event_004"
    event_data = None

    # Act
    event = KernelProcessEvent(id=event_id, data=event_data)

    # Assert
    assert event.id == event_id
    assert event.data is None
    assert event.visibility == KernelProcessEventVisibility.Internal

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
