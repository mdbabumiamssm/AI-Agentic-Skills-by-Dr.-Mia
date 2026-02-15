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

import logging
from dataclasses import dataclass

from ._checkpoint import WorkflowCheckpoint
from ._const import EXECUTOR_STATE_KEY
from ._events import RequestInfoEvent

logger = logging.getLogger(__name__)


@dataclass
class WorkflowCheckpointSummary:
    """Human-readable summary of a workflow checkpoint."""

    checkpoint_id: str
    timestamp: str
    iteration_count: int
    targets: list[str]
    executor_ids: list[str]
    status: str
    pending_request_info_events: list[RequestInfoEvent]


def get_checkpoint_summary(checkpoint: WorkflowCheckpoint) -> WorkflowCheckpointSummary:
    targets = sorted(checkpoint.messages.keys())
    executor_ids = sorted(checkpoint.shared_state.get(EXECUTOR_STATE_KEY, {}).keys())
    pending_request_info_events = [
        RequestInfoEvent.from_dict(request) for request in checkpoint.pending_request_info_events.values()
    ]

    status = "idle"
    if pending_request_info_events:
        status = "awaiting request response"
    elif not checkpoint.messages and "finalise" in executor_ids:
        status = "completed"
    elif checkpoint.messages:
        status = "awaiting next superstep"

    return WorkflowCheckpointSummary(
        checkpoint_id=checkpoint.checkpoint_id,
        timestamp=checkpoint.timestamp,
        iteration_count=checkpoint.iteration_count,
        targets=targets,
        executor_ids=executor_ids,
        status=status,
        pending_request_info_events=pending_request_info_events,
    )

__AUTHOR_SIGNATURE__ = "9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE"
